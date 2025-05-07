#include <fstream>

#include "pzlog.h"
#include <iostream>
#include <string>

#include <cmath>
#include <set>
#include "TPZRefPatternDataBase.h"
#include "TPZRefPattern.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"

#include "TPZMultiphysicsCompMesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzgeoelbc.h"
#include "TPZApproxCreator.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZMixedSymElasticityND.h"
#include "TPZInterfaceSymTensor.h"

#include "pzelementgroup.h"
#include "pzcondensedcompel.h"

#include "pzfstrmatrix.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"  //symmetric sparse matrix storage
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZAnalyticSolution.h"

#ifdef PZ_LOG
static TPZLogger logger("testmhm");
#endif

/// uniformly refine the mesh
void UniformRefine(TPZGeoMesh *gmesh, int nref);

/// divide the volumetric elements with a specific refinement pattern
void CreateJohnsonMercier(TPZGeoMesh *gmesh);

/// Add geometric boundary elements
void AddWrapElements(TPZGeoMesh *gmesh);

/// create a discontinuous mesh with continuity of center nodes
TPZCompMesh *CreateTensorSpace(TPZGeoMesh *gmesh);

/// create the displacement space and lagrange multipliers
TPZCompMesh *CreateDisplacementSpace(TPZGeoMesh *gmesh);

TPZMultiphysicsCompMesh *GenerateMultiphysicsMesh(TPZGeoMesh *gmesh);

void AddInterfaceElements(TPZMultiphysicsCompMesh *mfmesh);

std::set<int64_t> InactiveInterfaceConnects(TPZMultiphysicsCompMesh *mfmesh);

void GroupElements(TPZMultiphysicsCompMesh *mfmesh);

void UniformStretch(TPZMultiphysicsCompMesh *mfmesh, TPZFMatrix<STATE> &sol);

void CheckMatrixConsistency(TPZMultiphysicsCompMesh *mfmesh);

void SolveProblem(TPZMultiphysicsCompMesh *mfmesh);

enum matids {matID = 1,matBCD = 2, matBCN = 3, matWrap, matIntFacePos, matIntFaceNeg, matLagrDisp};

TElasticity2DAnalytic gAnalytic;

int main(int argc, char *argv[]) {
    //    TPZMaterial::gBigNumber = 1.e16;
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    TPZGmshReader Square;
    Square.GetDimNamePhysical()[2]["domain"] = matID;
    Square.GetDimNamePhysical()[1]["dirichlet"] = matBCD;
    Square.GetDimNamePhysical()[1]["neuman"] = matBCN;
    gAnalytic.fProblemType = TElasticity2DAnalytic::Etest1;
#ifdef MACOSX
    TPZGeoMesh *gmesh = Square.GeometricGmshMesh("../TriangleBC.msh");
#else
    TPZGeoMesh *gmesh = Square.GeometricGmshMesh("Triangle.msh");
#endif
    UniformRefine(gmesh, 2);
    CreateJohnsonMercier(gmesh);
    AddWrapElements(gmesh);
    {
        std::ofstream out("Mercier.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    }
    auto mfmesh = GenerateMultiphysicsMesh(gmesh);
    mfmesh->ComputeNodElCon();
    {
        std::ofstream out("MFMesh.txt");
        mfmesh->Print(out);
    }
    GroupElements(mfmesh);
    {
        std::ofstream out("MFMeshCondensed.txt");
        mfmesh->Print(out);
    }
    SolveProblem(mfmesh);
    TPZManVector<STATE> errors(5,0.);
    mfmesh->ElementSolution().Redim(mfmesh->NElements(), 5);
    mfmesh->EvaluateError(1, errors);
    std::cout << "errors " << errors << std::endl;
    
    return 0;
}

TPZRefPattern *TriangleRef() {
    char buf[] =
    "4     4  "
    "100       TensTri    "
    "0.     0.     0. "
    "1.     0.     0. "
    "0.     1.     0. "
    "0.3333    0.3333     0. "
    "2     3     0     1     2 "
    "2     3     0     1     3 "
    "2     3     1     2     3 "
    "2     3     2     0     3 ";
    std::istringstream str(buf);
    TPZRefPattern *refpat = new TPZRefPattern(str);
    return refpat;
}

TPZRefPattern *QuadRef() {
    char buf[] =
    "5     5  "
    "110       TensQuad    "
    "-1.    -1.     0. "
    " 1.    -1.     0. "
    " 1.     1.     0. "
    "-1.     1.     0. "
    " 0.     0.     0. "
    " 3    4     0     1     2     3 "
    "2     3     0     1     4 "
    "2     3     1     2     4 "
    "2     3     2     3     4 "
    "2     3     3     0     4 ";
    std::istringstream str(buf);
    TPZRefPattern *refpat = new TPZRefPattern(str);
    return refpat;
}

/// divide the volumetric elements with a specific refinement pattern
void CreateJohnsonMercier(TPZGeoMesh *gmesh) {
    
    TPZAutoPointer<TPZRefPattern> Trirefpat = TriangleRef();
    TPZAutoPointer<TPZRefPattern> Quadrefpat = QuadRef();
    {
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->HasSubElement()) continue;
            TPZManVector<TPZGeoEl *> subels;
            if(gel->Type() == EQuadrilateral) {
                gel->SetRefPattern(Quadrefpat);
                gel->Divide(subels);
            } else if (gel->Type() == ETriangle) {
                gel->SetRefPattern(Trirefpat);
                gel->Divide(subels);
            }
        }
    }
}

/// create a discontinuous mesh with continuity of center nodes
TPZCompMesh *CreateTensorSpace(TPZGeoMesh *gmesh) {
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = gmesh->Dimension();
    cmesh->SetDefaultOrder(1);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    TPZNullMaterial<STATE> *null = new TPZNullMaterial<STATE>(matID,gmesh->Dimension(),3);
    cmesh->InsertMaterialObject(null);
    TPZNullMaterial<STATE> *materialwrap = new TPZNullMaterial<STATE>(matWrap,gmesh->Dimension()-1,3);
    
    cmesh->AutoBuild();
    
    cmesh->InsertMaterialObject(materialwrap);
    
    std::map<int64_t,int64_t> topconnectindex;

    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() != gmesh->Dimension()) continue;
        int ncorner = gel->NCornerNodes();
        int64_t cindex = cel->ConnectIndex(ncorner-1);
        int64_t nodeindex = gel->NodeIndex(ncorner-1);
        if(topconnectindex.find(nodeindex) == topconnectindex.end()) {
            topconnectindex[nodeindex] = cindex;
        } else {
            cindex = topconnectindex[nodeindex];
        }
        cel->SetConnectIndex(ncorner-1, cindex);
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    
    cmesh->Reference()->ResetReference();
    nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() != gmesh->Dimension()) continue;
        if(gel->HasSubElement()) DebugStop();
        int firstside = gel->FirstSide(dim-1);
        int nsides = gel->NSides();
        for (int side = firstside; side < nsides-1 ; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour(gelside.Neighbour());
            if(neighbour.Element()->MaterialId() != matWrap) DebugStop();
            cel->LoadElementReference();
            TPZCompEl *wrap =cmesh->ApproxSpace().CreateCompEl(neighbour.Element(), *cmesh);
            gel->ResetReference();
            wrap->Reference()->ResetReference();
        }
    }
    {
        std::ofstream out("tensormesh.txt");
        cmesh->Print(out);
    }
    return cmesh;
}

/// Add geometric boundary elements
void AddWrapElements(TPZGeoMesh *gmesh) {
    int64_t nel = gmesh->NElements();
    int dim = gmesh->Dimension();
    for (int64_t el = 0; el<nel ; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->HasSubElement()) continue;
        if(gel->Dimension() != dim) continue;
        int firstside = gel->FirstSide(dim-1);
        int nsides = gel->NSides();
        for (int side = firstside; side < nsides-1; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElBC gbc(gelside,matWrap);
            TPZGeoElSide wrapside = gbc.CreatedElement();
            TPZGeoElSide dispsside = wrapside.HasNeighbour(matLagrDisp);
            if(dispsside) {
                TPZGeoElBC gbc(wrapside,matIntFaceNeg);
            } else {
                TPZGeoElBC gbc(wrapside,matIntFacePos);
                TPZGeoElSide gelsideintpos = gbc.CreatedElement();
                TPZGeoElSide boundside = gelsideintpos.HasNeighbour({matBCD,matBCN});
                if(!boundside) {
                    TPZGeoElBC gbc2(gelsideintpos,matLagrDisp);
                }
            }
        }
    }
}

/// create the displacement space and lagrange multipliers
TPZCompMesh *CreateDisplacementSpace(TPZGeoMesh *gmesh) {
    int dim = gmesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(0);
    cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    
    TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial<STATE>(matID,dim,2);
    cmesh->InsertMaterialObject(nullmat);
    
    std::set<int> matids= {matID};
    cmesh->AutoBuild(matids);
    
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    gmesh->ResetReference();
    
    matids.clear();
    matids = {matBCD,matBCN,matLagrDisp};
    
    TPZNullMaterial<STATE> *nullBCD = new TPZNullMaterial<STATE>(matBCD,dim-1,2);
    cmesh->InsertMaterialObject(nullBCD);
    TPZNullMaterial<STATE> *nullBCN = new TPZNullMaterial<STATE>(matBCN,dim-1,2);
    cmesh->InsertMaterialObject(nullBCN);
    TPZNullMaterial<> *nullLagr = new TPZNullMaterial<>(matLagrDisp,dim-1,2);
    cmesh->InsertMaterialObject(nullLagr);
    cmesh->AutoBuild(matids);
    cmesh->ExpandSolution();
    {
        std::ofstream out("DisplacementMesh.txt");
        cmesh->Print(out);
    }
    return cmesh;
}

TPZMultiphysicsCompMesh *GenerateMultiphysicsMesh(TPZGeoMesh *gmesh)
{
    int dim = gmesh->Dimension();
    TPZMultiphysicsCompMesh *mfmesh = new TPZMultiphysicsCompMesh(gmesh);
    TPZManVector<TPZCompMesh *> meshvec(2,0);
    meshvec[0] = CreateTensorSpace(gmesh);
    meshvec[1] = CreateDisplacementSpace(gmesh);
    
    TPZNullMaterialCS<STATE> *dispinterface = new TPZNullMaterialCS<STATE> (matLagrDisp,dim-1, dim);
    mfmesh->InsertMaterialObject(dispinterface);
    gAnalytic.gE = 1.;
    gAnalytic.gPoisson = 0.;
    TPZMixedSymElasticityND *mixed = new TPZMixedSymElasticityND(matID,dim);
    mixed->SetElasticity(gAnalytic.gE, gAnalytic.gPoisson);
    mixed->SetForcingFunction(gAnalytic.ForceFunc(), 4);
    mixed->SetExactSol(gAnalytic.ExactSolution(), 4);
    mfmesh->InsertMaterialObject(mixed);
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<STATE> val2(2,0.);
    auto *bcD = mixed->CreateBC(mixed, matBCD, 0, val1, val2);
    bcD->SetForcingFunctionBC(gAnalytic.ExactSolution(), 3);
    mfmesh->InsertMaterialObject(bcD);
    auto *bcN = mixed->CreateBC(mixed, matBCN, 1, val1, val2);
    bcN->SetForcingFunctionBC(gAnalytic.ExactSolution(), 3);
    mfmesh->InsertMaterialObject(bcN);
    TPZNullMaterialCS<STATE> *materialwrap = new TPZNullMaterialCS<STATE>(matWrap,dim-1,dim);
    mfmesh->InsertMaterialObject(materialwrap);
    mfmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    mfmesh->BuildMultiphysicsSpace(meshvec);
    AddInterfaceElements(mfmesh);
    return mfmesh;
}

void AddInterfaceElements(TPZMultiphysicsCompMesh *mfmesh) {
    TPZGeoMesh *gmesh = mfmesh->Reference();
    int dim = gmesh->Dimension();
    TPZInterfaceSymTensor *sym1 = new TPZInterfaceSymTensor(matIntFacePos, dim-1);
    sym1->SetMultiplier(1);
    mfmesh->InsertMaterialObject(sym1);
    TPZInterfaceSymTensor *sym2 = new TPZInterfaceSymTensor(matIntFaceNeg, dim-1);
    sym2->SetMultiplier(-1);
    mfmesh->InsertMaterialObject(sym2);
    
    mfmesh->LoadReferences();
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel ; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() == matWrap) {
            TPZGeoElSide gelside(gel);
            TPZCompElSide celsidewrap = gelside.Reference();
            TPZGeoElSide neighint = gelside.Neighbour();
            TPZGeoEl *neighgel = neighint.Element();
            if(neighgel->MaterialId() != matIntFaceNeg && neighgel->MaterialId() != matIntFacePos) DebugStop();
            TPZGeoElSide neighdisp = gelside.HasNeighbour(matLagrDisp);
            if(neighdisp) {
                auto intface = new TPZMultiphysicsInterfaceElement(*mfmesh,neighgel,celsidewrap,neighdisp.Reference());
                continue;
            }
            neighdisp = gelside.HasNeighbour({matBCD,matBCN});
            if(!neighdisp) DebugStop();
            auto intface = new TPZMultiphysicsInterfaceElement(*mfmesh,neighgel,celsidewrap,neighdisp.Reference());
        }
    }
    
}

std::set<int64_t> InactiveInterfaceConnects(TPZMultiphysicsCompMesh *mfmesh) {
    std::set<int64_t> inactive;
    TPZGeoMesh *gmesh = mfmesh->Reference();
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->HasSubElement()) continue;
        if(gel->MaterialId() != matID) continue;
        int nnodes = gel->NNodes();
        int side = nnodes-1;
        TPZGeoElSide gelside(gel,side);
        for(auto neigh = gelside.Neighbour(); neigh!=gelside; neigh++) {
            TPZGeoEl *neighgel = neigh.Element();
            if(neighgel->MaterialId() == matLagrDisp) {
                TPZCompEl *cel = neighgel->Reference();
                int neighside = neigh.Side();
                inactive.insert(cel->ConnectIndex(neighside));
            }
        }
    }
    mfmesh->CleanUpUnconnectedNodes();
    return inactive;
}

void GroupElements(TPZMultiphysicsCompMesh *mfmesh) {
    /// map from geometric node to set of geometric elements
    std::map<int64_t,std::set<int64_t> > nodetoel;
    auto inactive = InactiveInterfaceConnects(mfmesh);

    
    TPZGeoMesh *gmesh = mfmesh->Reference();
    int dim = gmesh->Dimension();
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->HasSubElement()) continue;
        if(gel->Dimension() != dim) continue;
        int nnodes = gel->NCornerNodes();
        int side = nnodes-1;
        int64_t nodeindex = gel->NodeIndex(side);
        if(nodetoel.find(nodeindex) != nodetoel.end()) continue;
        nodetoel[nodeindex].insert(el);
//        std::cout << "Including " << el << " matid " << gel->MaterialId() << std::endl;
        TPZGeoElSide gelside(gel,side);
        for (TPZGeoElSide neigh = gelside.Neighbour(); neigh != gelside; neigh++) {
            TPZGeoEl *gelneigh = neigh.Element();
            if(gelneigh->HasSubElement()) continue;
            int64_t neighindex = gelneigh->Index();
//            std::cout << "Including " << gelneigh->Index() << " matid " << gelneigh->MaterialId() << std::endl;
            nodetoel[nodeindex].insert(neighindex);
        }
    }
    // for the volumetric elements, include the wrap and interface element
    for (auto &it : nodetoel) {
        int64_t nodeindex = it.first;
        std::set<int64_t> elset = it.second;
        for (auto elindex : elset) {
            TPZGeoEl *gel = gmesh->Element(elindex);
            if(gel->Dimension() != dim) continue;
            int firstside = gel->FirstSide(dim-1);
            int lastside = gel->NSides()-1;
            for(int side = firstside; side < lastside; side++) {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neigh1 = gelside.Neighbour();
                TPZGeoElSide neigh2 = neigh1.Neighbour();
                it.second.insert(neigh1.Element()->Index());
                it.second.insert(neigh2.Element()->Index());
            }
        }
    }
    std::set<TPZElementGroup *> groupset;
    for (auto it : nodetoel) {
        int64_t nodeindex = it.first;
        std::set<int64_t> &elset = it.second;
        TPZElementGroup *elgroup = new TPZElementGroup(*mfmesh);
        groupset.insert(elgroup);
        for(auto gelindex : elset) {
            TPZGeoEl *gel = gmesh->Element(gelindex);
            TPZCompEl *cel = gel->Reference();
            if(!cel) DebugStop();
            elgroup->AddElement(cel);
        }
    }
    mfmesh->ComputeNodElCon();
    for(auto it : inactive) {
        TPZConnect &c = mfmesh->ConnectVec()[it];
        c.IncrementElConnected();
    }
    for(auto it : groupset) {
        TPZCondensedCompElT<STATE> *compel = new TPZCondensedCompElT<STATE>(it);
    }
    for(auto it : inactive) {
        TPZConnect &c = mfmesh->ConnectVec()[it];
        c.SetCondensed(true);
    }
    mfmesh->CleanUpUnconnectedNodes();
}

void UniformStretch(TPZMultiphysicsCompMesh *mfmesh, TPZFMatrix<STATE> &sol) {
    auto &meshvec = mfmesh->MeshVector();
    auto stressmesh = meshvec[0];
    auto dispmesh = meshvec[1];
    TPZFMatrix<STATE> &solstress = stressmesh->Solution();
    int64_t nstr = solstress.Rows();
    for(int64_t i=0; i<nstr; i+=3) solstress(i,0) = 1.;
    int64_t neldisp = dispmesh->NElements();
    TPZFMatrix<STATE> &dispsol = dispmesh->Solution();
    for(int64_t el = 0; el<neldisp; el++) {
        TPZCompEl *cel = dispmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() == 1) {
            for(int n = 0; n<2; n++) {
                double x = gel->NodePtr(n)->Coord(0);
                TPZConnect &c = cel->Connect(n);
                int64_t seq = c.SequenceNumber();
                int64_t eq = dispmesh->Block().Index(seq, 0);
                dispsol(eq,0) = x;
            }
        } else if(gel->Dimension()==2) {
            TPZManVector<REAL,3> xi(2,0.),x(3,0.);
            gel->CenterPoint(gel->NSides()-1, xi);
            gel->X(xi, x);
            int nc = cel->NConnects();
            TPZConnect &c = cel->Connect(nc-1);
            int64_t seq = c.SequenceNumber();
            int64_t eq = dispmesh->Block().Index(seq, 0);
            dispsol(eq,0) = x[0];

        }
    }
    mfmesh->LoadSolutionFromMeshes();
    sol = mfmesh->Solution();
}
void CheckMatrixConsistency(TPZMultiphysicsCompMesh *mfmesh) {
    TPZFMatrix<STATE> solglob;
    UniformStretch(mfmesh, solglob);
    TPZFStructMatrix<> strmat(mfmesh);
    TPZFMatrix<STATE> rhs;
    auto glob = strmat.CreateAssemble(rhs);
    TPZFMatrix<STATE> *fglob = dynamic_cast<TPZFMatrix<STATE>*>(glob);
    {
        std::ofstream globmat("glob.txt");
        fglob->Print("GK = ",globmat,EMathematicaInput);
        solglob.Print("Sol = ",globmat,EMathematicaInput);
    }

}

#include "TPZVTKGenerator.h"

void SolveProblem(TPZMultiphysicsCompMesh *mfmesh) {
    TPZLinearAnalysis analysis(mfmesh);
    TPZSkylineStructMatrix<STATE> matskl(mfmesh);
    matskl.SetNumThreads(0);
    analysis.SetStructuralMatrix(matskl);

    /// Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt); // ELU //ECholesky // ELDLt
    analysis.SetSolver(step);
    
    analysis.Run();

    TPZManVector<std::string> fields = {"Displacement", "SigmaX", "SigmaY", "TauXY"};
    TPZVTKGenerator vtk(mfmesh, fields, "Solution.vtk", 0);
    vtk.SetNThreads(0);
    vtk.Do();
    
}

#include "pzcheckgeom.h"
/// uniformly refine the mesh
void UniformRefine(TPZGeoMesh *gmesh, int nref) {
    TPZCheckGeom check(gmesh);
    check.UniformRefine(nref);
}

