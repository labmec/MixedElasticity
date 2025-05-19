//
//  MeshConditioning.cpp
//  SymTensor
//
//  Created by Philippe Devloo on 13/05/25.
//

#include "MeshConditioning.h"
#include "TPZLagrangeMultiplierCS.h"

TElasticity2DAnalytic gAnalytic;

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
TPZCompMesh *CreateTensorSpace(TPZGeoMesh *gmesh, int porder, int porderlow) {
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = gmesh->Dimension();
    cmesh->SetDefaultOrder(porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    TPZNullMaterial<STATE> *null = new TPZNullMaterial<STATE>(matID,gmesh->Dimension(),3);
    cmesh->InsertMaterialObject(null);
    TPZNullMaterial<STATE> *materialwrap = new TPZNullMaterial<STATE>(matWrap,gmesh->Dimension()-1,3);
    
    cmesh->AutoBuild();
    
    cmesh->InsertMaterialObject(materialwrap);
    
    std::map<int64_t,int64_t> topconnectindex;
    std::set<int> matids = {matHEXAGON};
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() != gmesh->Dimension()) continue;
        int ncorner = gel->NCornerNodes();
        int ntop = 0;
        for(int side = 0; side<ncorner; side++) {
            TPZGeoElSide gelside(gel,side);
            if(gelside.HasNeighbour(matids)) continue;
            ntop++;
            int64_t cindex = cel->ConnectIndex(side);
            int64_t nodeindex = gel->NodeIndex(side);
            if(topconnectindex.find(nodeindex) == topconnectindex.end()) {
                topconnectindex[nodeindex] = cindex;
            } else {
                cindex = topconnectindex[nodeindex];
            }
            cel->SetConnectIndex(side, cindex);

        }
        if(ntop != 1) DebugStop();
//        int64_t cindex = cel->ConnectIndex(ncorner-1);
//        int64_t nodeindex = gel->NodeIndex(ncorner-1);
//        if(topconnectindex.find(nodeindex) == topconnectindex.end()) {
//            topconnectindex[nodeindex] = cindex;
//        } else {
//            cindex = topconnectindex[nodeindex];
//        }
//        cel->SetConnectIndex(ncorner-1, cindex);
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    
    cmesh->Reference()->ResetReference();
    nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() == gmesh->Dimension()) {
            cmesh->SetDefaultOrder(porder);
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
    }
    if(porderlow >= 0) {
        TPZNullMaterial<STATE> *materialtens = new TPZNullMaterial<STATE>(matTensionLoworder,gmesh->Dimension()-1,2);
        cmesh->InsertMaterialObject(materialtens);

        cmesh->SetDefaultOrder(porderlow);
        // if the order is zero, discontinuous elements need be created
        if(porderlow == 0) {
            cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
            cmesh->SetDimModel(dim-1);
        }
        nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->MaterialId() == matTensionLoworder) {
                cmesh->ApproxSpace().CreateCompEl(gel, *cmesh);
                gel->ResetReference();
            }
        }
        cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        cmesh->SetDimModel(dim);
    }
    cmesh->ExpandSolution();
    if(0)
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
    if(dim != 2) DebugStop();
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
            } else { // there is no displacement boundary element
                TPZGeoElBC gbc(wrapside,matIntFacePos);
                TPZGeoElSide gelsideintpos = gbc.CreatedElement();
                TPZGeoElBC gbc2(gelsideintpos,matLagrDisp);
            }
        }
    }
    // testing the integrity of the mesh
    std::set<int> matids = {matLagrDisp, matBCD, matBCN};
    for (int64_t el = 0; el<nel ; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->HasSubElement()) continue;
        if(gel->Dimension() != dim) continue;
        int firstside = gel->FirstSide(dim-1);
        int nsides = gel->NSides();
        for (int side = firstside; side < nsides-1; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neigh = gelside.Neighbour();
            if(neigh.Element()->MaterialId() != matWrap) DebugStop();
            neigh++;
            int neighmat = neigh.Element()->MaterialId();
            if(neighmat != matIntFaceNeg && neighmat != matIntFacePos) DebugStop();
            if(neighmat == matIntFacePos) {
                if(!neigh.HasNeighbour(matIntFacePos)) DebugStop();
            } else if (neighmat == matIntFaceNeg) {
                if(!neigh.HasNeighbour(matIntFacePos)) DebugStop();
            }
            if(!neigh.HasNeighbour(matids)) DebugStop();
        }
    }
}

/// Add geometric boundary elements to represent triple hybridization
void AddWrapElementsTripleHybrid(TPZGeoMesh *gmesh)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() != matWrap) continue;
        TPZGeoElSide gelwrap(gel);
        // we create double hybrid elements only along the boundary of the polygons
        if(!gelwrap.HasNeighbour(matHEXAGON)) continue;
        TPZGeoElSide intface = gelwrap.Neighbour();
        int matintface = intface.Element()->MaterialId();
        if(matintface == matIntFacePos) {
            TPZGeoElSide lagrHigh = intface.Neighbour();
            if(lagrHigh.Element()->MaterialId() != matLagrDisp) DebugStop();
            TPZGeoElBC gbcIntface(lagrHigh,matIntFaceNegLoworder);
            TPZGeoElBC gbctenslow(gbcIntface,matTensionLoworder);
            TPZGeoElBC gbcIntPos(gbctenslow,matIntFacePosLoworder);
            TPZGeoElBC gbcDispLow(gbcIntPos,matLagrDispLoworder);
        } else {
            TPZGeoElBC gbcDisp(intface,matLagrDisp);
            TPZGeoElBC gbcIntface(gbcDisp,matIntFacePosLoworder);
            TPZGeoElBC gbctenslow(gbcIntface,matTensionLoworder);
            TPZGeoElBC gbcIntNeg(gbctenslow,matIntFaceNegLoworder);
        }
    }
    
}

/// create the displacement space and lagrange multipliers
TPZCompMesh *CreateDisplacementSpace(TPZGeoMesh *gmesh, int porder, int porderlow) {
    int dim = gmesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(porder-1);
    if(porder == 1) {
        cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    } else {
        cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    
    TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial<STATE>(matID,dim,2);
    cmesh->InsertMaterialObject(nullmat);
    
    std::set<int> matids= {matID};
    cmesh->AutoBuild(matids);
    
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->SetDefaultOrder(porder);
    gmesh->ResetReference();
    
    matids.clear();
    
    TPZNullMaterial<STATE> *nullBCD = new TPZNullMaterial<STATE>(matBCD,dim-1,2);
    cmesh->InsertMaterialObject(nullBCD);
    TPZNullMaterial<STATE> *nullBCN = new TPZNullMaterial<STATE>(matBCN,dim-1,2);
    cmesh->InsertMaterialObject(nullBCN);
    TPZNullMaterial<> *nullLagr = new TPZNullMaterial<>(matLagrDisp,dim-1,2);
    cmesh->InsertMaterialObject(nullLagr);
    matids = {matLagrDisp};
    cmesh->AutoBuild(matids);
    if(porderlow >= 0) {
        TPZNullMaterial<> *nullLagr = new TPZNullMaterial<>(matLagrDispLoworder,dim-1,2);
        cmesh->InsertMaterialObject(nullLagr);
        if(porderlow == 0) {
            cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
            cmesh->SetDimModel(dim-1); // otherwise dim-1 elements will not be created
        } else {
            cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        }
        cmesh->SetDefaultOrder(porderlow);
        matids = {matLagrDispLoworder};
        cmesh->AutoBuild(matids);

    }
    matids.clear();
    matids = {matBCD,matBCN};
    // create boundary condition elements with the same connects as the neighbouring displacement element
    int64_t nel = cmesh->NElements();
    for(int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->MaterialId() != matLagrDisp && porderlow < 0) continue;
        if(gel->MaterialId() != matLagrDispLoworder && porderlow >= 0) continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neigh = gelside.HasNeighbour(matids);
        if(neigh) {
            cel->LoadElementReference();
            TPZCompEl *celbc = cmesh->ApproxSpace().CreateCompEl(neigh.Element(), *cmesh);
            // if the order is zero, the discontinuous element will not adopt the connect index of the neighbour
            if(porderlow == 0 && celbc->NConnects() != 1) DebugStop();
            if(porderlow == 0) {
                celbc->SetConnectIndex(0, cel->ConnectIndex(0));
            }
            gel->ResetReference();
            neigh.Element()->ResetReference();
        }
    }
    if(porderlow == 0) {
        cmesh->ComputeNodElCon();
        cmesh->CleanUpUnconnectedNodes();
    }
    cmesh->ExpandSolution();
    if(0)
    {
        std::ofstream out("DisplacementMesh.txt");
        cmesh->Print(out);
    }
    return cmesh;
}

TPZMultiphysicsCompMesh *GenerateMultiphysicsMesh(TPZGeoMesh *gmesh, int porder, int porderlow)
{
    int dim = gmesh->Dimension();
    TPZMultiphysicsCompMesh *mfmesh = new TPZMultiphysicsCompMesh(gmesh);
    TPZManVector<TPZCompMesh *> meshvec(2,0);
    meshvec[0] = CreateTensorSpace(gmesh, porder, porderlow);
    meshvec[1] = CreateDisplacementSpace(gmesh, porder, porderlow);
    
    TPZNullMaterialCS<STATE> *dispinterface = new TPZNullMaterialCS<STATE> (matLagrDisp,dim-1, dim);
    mfmesh->InsertMaterialObject(dispinterface);
    TPZMixedSymElasticityND *mixed = new TPZMixedSymElasticityND(matID,dim);
    mixed->SetPlaneStress();
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
    if(porderlow >= 0) {
        TPZNullMaterialCS<STATE> *materialtenslow = new TPZNullMaterialCS<STATE>(matTensionLoworder,dim-1,dim);
        mfmesh->InsertMaterialObject(materialtenslow);
        TPZNullMaterialCS<STATE> *materialdisplow = new TPZNullMaterialCS<STATE>(matLagrDispLoworder,dim-1,dim);
        mfmesh->InsertMaterialObject(materialdisplow);
    }
    mfmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    mfmesh->BuildMultiphysicsSpace(meshvec);
    AddInterfaceElements(mfmesh,porderlow);
    return mfmesh;
}

void AddInterfaceElements(TPZMultiphysicsCompMesh *mfmesh, int porderlow) {
    TPZGeoMesh *gmesh = mfmesh->Reference();
    int dim = gmesh->Dimension();
    TPZInterfaceSymTensor *sym1 = new TPZInterfaceSymTensor(matIntFacePos, dim-1);
    sym1->SetMultiplier(1);
    mfmesh->InsertMaterialObject(sym1);
    TPZInterfaceSymTensor *sym2 = new TPZInterfaceSymTensor(matIntFaceNeg, dim-1);
    sym2->SetMultiplier(1);
    mfmesh->InsertMaterialObject(sym2);
    
    // indicates that we will use triple hybridization
    if(porderlow >= 0) {
        int nstate = 2;
        TPZLagrangeMultiplierCS<STATE> *matintpos = new TPZLagrangeMultiplierCS<STATE>(matIntFacePosLoworder, dim-1, nstate);
        matintpos->SetMultiplier(1);
        mfmesh->InsertMaterialObject(matintpos);
        TPZLagrangeMultiplierCS<STATE> *matintneg = new TPZLagrangeMultiplierCS<STATE>(matIntFaceNegLoworder, dim-1, nstate);
        matintneg->SetMultiplier(-1);
        mfmesh->InsertMaterialObject(matintneg);
    }
    
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
            if(!neighdisp) DebugStop();
            auto intface = new TPZMultiphysicsInterfaceElement(*mfmesh,neighgel,celsidewrap,neighdisp.Reference());
            if(porderlow >= 0 && gelside.HasNeighbour(matHEXAGON)) {
                
                TPZGeoElSide gelintfaceloworder1 = neighdisp.Neighbour();
                int intfacematid = gelintfaceloworder1.Element()->MaterialId();
                if(intfacematid != matIntFaceNegLoworder && intfacematid != matIntFacePosLoworder) DebugStop();
                TPZGeoElSide geltractionloworder = gelintfaceloworder1.Neighbour();
                if(geltractionloworder.Element()->MaterialId() != matTensionLoworder) DebugStop();
                new TPZMultiphysicsInterfaceElement(*mfmesh,gelintfaceloworder1.Element(),neighdisp.Reference(),geltractionloworder.Reference());
                TPZGeoElSide gelintfaceloworder2 = geltractionloworder.Neighbour();
                intfacematid = gelintfaceloworder2.Element()->MaterialId();
                if(intfacematid != matIntFaceNegLoworder && intfacematid != matIntFacePosLoworder) DebugStop();
                TPZGeoElSide geldisploworder = gelintfaceloworder2.HasNeighbour(matLagrDispLoworder);
                if(!geldisploworder) DebugStop();
                new TPZMultiphysicsInterfaceElement(*mfmesh,gelintfaceloworder2.Element(),geltractionloworder.Reference(),geldisploworder.Reference());
            }
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
        int firstside = 0;
        int lastside = gel->NCornerNodes();
        for (int side = firstside; side <lastside; side++) {
            TPZGeoElSide gelside(gel,side);
            if(gelside.HasNeighbour(matHEXAGON)) continue;
            for(auto neigh = gelside.Neighbour(); neigh!=gelside; neigh++) {
                TPZGeoEl *neighgel = neigh.Element();
                if(neighgel->MaterialId() == matLagrDisp) {
                    TPZCompEl *cel = neighgel->Reference();
                    int neighside = neigh.Side();
                    inactive.insert(cel->ConnectIndex(neighside));
                }
            }
        }
    }
    mfmesh->CleanUpUnconnectedNodes();
    return inactive;
}

void GroupElements(TPZMultiphysicsCompMesh *mfmesh, int ploworder) {
    /// map from geometric node to set of geometric elements
    std::map<int64_t,std::set<int64_t> > nodetoel;
    auto inactive = InactiveInterfaceConnects(mfmesh);

    // build a data structure from the geometric node index to the computational element
    TPZGeoMesh *gmesh = mfmesh->Reference();
    if(0)
    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);
    }
    int dim = gmesh->Dimension();
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->HasSubElement()) continue;
        if(gel->Dimension() != dim) continue;
        int nnodes = gel->NCornerNodes();
        int side = -1;
        // look for a corner node that is not connected to the boundary of the polygon
        for(side = 0 ; side<nnodes; side++) {
            TPZGeoElSide gelside(gel,side);
            auto neigh = gelside.HasNeighbour(matHEXAGON);
            if(!neigh) break;
        }
        if(side == nnodes) DebugStop();
        int64_t nodeindex = gel->NodeIndex(side);
        // the element set has already been formed
        if(nodetoel.find(nodeindex) != nodetoel.end()) continue;
        nodetoel[nodeindex].insert(el);
//        std::cout << "Including " << el << " matid " << gel->MaterialId() << std::endl;
        // include ALL the neighbouring elements
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
                if(neigh1.Element()->MaterialId() != matWrap) DebugStop();
                TPZGeoElSide neigh2 = neigh1.Neighbour();
                int neigh2matid = neigh2.Element()->MaterialId();
                if(neigh2matid != matIntFaceNeg && neigh2matid != matIntFacePos) DebugStop();
                it.second.insert(neigh1.Element()->Index());
                it.second.insert(neigh2.Element()->Index());
                if(ploworder >= 0 && gelside.HasNeighbour(matHEXAGON)) {
                    TPZGeoElSide neighdisp = neigh2.Neighbour();
                    int dispmat = neighdisp.Element()->MaterialId();
                    if(dispmat != matLagrDisp) DebugStop();
                    it.second.insert(neighdisp.Element()->Index());
                    TPZGeoElSide neighIntLow1 = neighdisp.Neighbour();
                    int matint1 = neighIntLow1.Element()->MaterialId();
                    if(matint1 != matIntFaceNegLoworder && matint1 != matIntFacePosLoworder) DebugStop();
                    it.second.insert(neighIntLow1.Element()->Index());
                    TPZGeoElSide neighTensionLow = neighIntLow1.Neighbour();
                    int matstress = neighTensionLow.Element()->MaterialId();
                    if(matstress != matTensionLoworder) DebugStop();
                    it.second.insert(neighTensionLow.Element()->Index());
                    TPZGeoElSide neighIntLow2 = neighTensionLow.Neighbour();
                    int matint2 = neighIntLow2.Element()->MaterialId();
                    if(matint2 != matIntFaceNegLoworder && matint2 != matIntFacePosLoworder) DebugStop();
                    it.second.insert(neighIntLow2.Element()->Index());
                    // matDisp, matIntfacePosloworder, matTraction, matIntfacPosLoworder
                }
            }
        }
    }
    std::set<TPZElementGroup *> groupset;
    for (auto it : nodetoel) {
        int64_t nodeindex = it.first;
        std::set<int64_t> &elset = it.second;
        TPZElementGroup *elgroup = new TPZElementGroup(*mfmesh);
        groupset.insert(elgroup);
        std::set<int64_t> connectindexes;
        for(auto gelindex : elset) {
            TPZGeoEl *gel = gmesh->Element(gelindex);
            TPZCompEl *cel = gel->Reference();
            {
                int nc = cel->NConnects();
                for (int ic = 0; ic<nc; ic++) {
                    int64_t cindex = cel->ConnectIndex(ic);
                    connectindexes.insert(cindex);
                }
            }
            if(!cel) DebugStop();
            elgroup->AddElement(cel);
        }
        int64_t ngroupeq = 0;
        for(auto it : connectindexes) {
            if(inactive.find(it) != inactive.end()) continue;
            TPZConnect &c = mfmesh->ConnectVec()[it];
            ngroupeq += c.NState()*c.NShape();
        }
//        std::cout << "Number of equations associated with group index " << elgroup->Index() << " is " << ngroupeq << std::endl;
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

#include "pzcheckgeom.h"
/// uniformly refine the mesh
void UniformRefine(TPZGeoMesh *gmesh, int nref) {
    TPZCheckGeom check(gmesh);
    check.UniformRefine(nref);
}

void AddBCHoneycomb(TPZGeoMesh *gmesh, int bcidD, int bcidN) {
    TPZManVector<REAL,3> minimum(3),maximum(3);
    gmesh->NodeVec()[0].GetCoordinates(minimum);
    maximum = minimum;
    int64_t nnodes = gmesh->NNodes();
    for (int64_t n = 0; n<nnodes; n++) {
        for(int i = 0; i<3; i++) {
            REAL x = gmesh->NodeVec()[n].Coord(i);
            minimum[i] = minimum[i] > x ? x : minimum[i];
            maximum[i] = maximum[i] > x ? maximum[i] : x;
        }
    }
    TPZManVector<REAL,3> delta(3,0.);
    for(int i=0; i<3; i++) delta[i] = maximum[i]-minimum[i];
    for (int64_t n = 0; n<nnodes; n++) {
        TPZManVector<REAL,3> newco(3);
        for(int i = 0; i<2; i++) {
            REAL x = gmesh->NodeVec()[n].Coord(i);
            newco[i] = (x-minimum[i])/delta[i];
        }
        gmesh->NodeVec()[n].SetCoord(newco);
    }
    int64_t nel = gmesh->NElements();
    for(int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Dimension() != 2) continue;
        TPZGeoElSide gelside(gel);
        if(gelside.Neighbour() != gelside) {
            gel->RemoveConnectivities();
            delete gel;
        }
    }
    
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) continue;
        if(gel->MaterialId() != matHEXAGON) continue;
        int numd2 = 0;
        TPZGeoElSide gelside(gel);
        for (TPZGeoElSide neigh = gelside.Neighbour(); neigh != gelside; neigh++) {
            if(neigh.Element()->MaterialId() == matID) numd2++;
        }
        if(numd2 == 1) {
            REAL xco = gel->NodePtr(0)->Coord(0);
            if(xco < 0.5)
            {
                TPZGeoElBC(gelside,bcidD);
            }
            else {
                TPZGeoElBC(gelside,bcidN);
            }
        }
    }
}

/// Add matHexagon elements indicating the contour of the Johnson Mercier elements
void AddMatHexagon(TPZGeoMesh *gmesh) {
    int64_t nel = gmesh->NElements();
    std::set<int> bcids = {matHEXAGON};
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->Dimension() != 2) continue;
        int firstside = gel->FirstSide(1);
        int lastside = gel->NSides()-1;
        for (int side = firstside; side<lastside; side++) {
            TPZGeoElSide gelside(gel,side);
            if(gelside.HasNeighbour(bcids)) continue;
            TPZGeoElBC(gelside,matHEXAGON);
        }
    }
}
