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
#include "MeshConditioning.h"

#ifdef PZ_LOG
static TPZLogger logger("testmhm");
#endif


void UniformStretch(TPZMultiphysicsCompMesh *mfmesh, TPZFMatrix<STATE> &sol);

void CheckMatrixConsistency(TPZMultiphysicsCompMesh *mfmesh);

void SolveProblem(TPZMultiphysicsCompMesh *mfmesh);

TPZGeoMesh *OnePolygon(int nfaces);

void CheckOneElement(int nfaces, int porder, int porderlow);


int main(int argc, char *argv[]) {
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
//    CheckOneElement(4, 1, -1);
//    return 0;
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gAnalytic.fProblemType = TElasticity2DAnalytic::Etest2;
    gAnalytic.gE = 1.;
    gAnalytic.gPoisson = 0.3;
    int porder = 2;
    int porderlow = -1;
    //    TPZMaterial::gBigNumber = 1.e16;
    TPZGmshReader Square;
    Square.GetDimNamePhysical()[2]["domain"] = matID;
    Square.GetDimNamePhysical()[1]["dirichlet"] = matBCD;
    Square.GetDimNamePhysical()[1]["neuman"] = matBCN;
    std::string meshfile = "TriangleBC.msh";
//    std::string meshfile = "Quad.msh";
#ifdef MACOSX
    TPZGeoMesh *gmeshorig = Square.GeometricGmshMesh("../"+meshfile);
//    TPZGeoMesh *gmeshorig = Square.GeometricGmshMesh("../Quad.msh");
//    TPZGeoMesh *gmeshorig = OneTriangle();
#else
    TPZGeoMesh *gmeshorig = Square.GeometricGmshMesh(meshfile);
#endif
    {
        std::ofstream out("gmesh.txt");
        gmeshorig->Print(out);
    }
    {
        std::ofstream out("original.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshorig, out, true);
    }
    std::cout << meshfile << " Internal order " << porder << " Boundary order " << (porderlow == -1 ? porder : porderlow) << std::endl;
    UniformRefine(gmeshorig, 1);
    for(int iref = 0; iref<4; iref++) {
        UniformRefine(gmeshorig, 1);
        AddMatHexagon(gmeshorig);
        TPZGeoMesh *gmesh = new TPZGeoMesh(*gmeshorig);
        TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
        CreateJohnsonMercier(gmesh);
        AddWrapElements(gmesh);
        if(porderlow >= 0) {
            AddWrapElementsTripleHybrid(gmesh);
        }
        if(1)
        {
            std::ofstream out("Mercier.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
//            std::ofstream out3("gmesh.txt");
//            gmesh->Print(out3);
        }
        auto mfmesh = GenerateMultiphysicsMesh(gmesh, porder, porderlow);
//        mfmesh->ComputeNodElCon();
        if(0)
        {
            std::ofstream out("MFMesh.txt");
            mfmesh->Print(out);
        }
        GroupElements(mfmesh,porderlow);
        SolveProblem(mfmesh);
        if(0)
        {
            std::ofstream out("MFMeshCondensed.txt");
            mfmesh->Print(out);
        }
        TPZManVector<STATE> errors(5,0.);
        mfmesh->ElementSolution().Redim(mfmesh->NElements(), 5);
        mfmesh->EvaluateError(1, errors);
        auto *mat = mfmesh->FindMaterial(1);
        auto matmix = dynamic_cast<TPZMixedSymElasticityND *>(mat);
        auto names = matmix->ErrorNames();
        int64_t neq_cond = mfmesh->NEquations();
        int64_t neq = mfmesh->Solution().Rows();
        std::cout << "Global system neq " << neq_cond << " Total equations " << neq << std::endl;
        for(int i=0; i<names.size(); i++) std::cout << names[i] << " ";
        std::cout << std::endl;
        std::cout << "errors " << errors << std::endl;
    }
    
    return 0;
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
//    analysis.Assemble();
//    return;
    analysis.Run();

    TPZManVector<std::string> fields = {"Displacement", "SigmaX", "SigmaY", "TauXY"};
    TPZVTKGenerator vtk(mfmesh, fields, "Solution.vtk", 0);
    vtk.SetNThreads(0);
    vtk.Do();
    
}

TPZGeoMesh *OnePolygon(int nfaces) {
    int dim = 2;
    REAL angle = 2.*M_PI/nfaces;
    REAL radius = 0.5/sin(angle/2.);
    REAL height = radius*cos(angle/2.);
    REAL firstalpha = -M_PI_2-angle/2.;
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    TPZManVector<REAL,3> center = {0.,height,0.};
    {
        int64_t in = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[in].Initialize(center, *gmesh);
    }
    for(int i = 0; i<nfaces; i++) {
        TPZManVector<REAL,3> co(3,0.);
        REAL alpha = firstalpha+i*angle;
        co[0] = center[0]+radius*cos(alpha);
        co[1] = center[1]+radius*sin(alpha);
        int64_t in = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[in].Initialize(co, *gmesh);
    }
    for(int i = 0; i<nfaces; i++) {
        TPZManVector<int64_t,3> nodeindexes = {i+1,i+2,0};
        if(i==nfaces-1) nodeindexes[1] = 1;
        int64_t index;
        auto geoel = gmesh->CreateGeoElement(ETriangle, nodeindexes, matID, index);
    }
    gmesh->BuildConnectivity();
    for(int el = 0; el<nfaces; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int side = 3;
        TPZGeoElBC gbc1(gel,side,matBCN);
        TPZGeoElBC gbc2(gel,side,matHEXAGON);
//        std::cout << gbc.CreatedElement()->Dimension() << std::endl;
    }
    return gmesh;
}

void GetRigidBodyModes(TPZCompMesh *cmesh, TPZFMatrix<STATE> &rgb) {
    int64_t neq = cmesh->NEquations();
    TPZBlock &block = cmesh->Block();
    rgb.Redim(neq, 3);
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        int nc = gel->NCornerNodes();
        for(int corner = 0; corner<nc; corner++) {
            TPZManVector<REAL,3> co(3,0.);
            gel->NodePtr(corner)->GetCoordinates(co);
            TPZConnect &c = cel->Connect(corner);
            int64_t seq = c.SequenceNumber();
            int64_t pos = block.Position(seq);
            // x translation
            rgb(pos,0) = 1.;
            // y translation
            rgb(pos+1,1) = 1.;
            // rotation
            rgb(pos,2) = -co[1];
            rgb(pos+1,2) = co[0];
        }
    }
}

void CheckOneElement(int nfaces, int porder, int porderlow) {
    gAnalytic.fProblemType = TElasticity2DAnalytic::EDispx;
    gAnalytic.gE = 1.;
    gAnalytic.gPoisson = 0.3;

    TPZGeoMesh *gmeshorig = OnePolygon(nfaces);
    {
        std::ofstream out("gmesh.txt");
        gmeshorig->Print(out);
    }
    {
        std::ofstream out("original.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshorig, out, true);
    }
    TPZGeoMesh *gmesh = gmeshorig;
    AddWrapElements(gmesh);
    if(porderlow >= 0) {
        AddWrapElementsTripleHybrid(gmesh);
    }
    {
        std::ofstream out("Mercier.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
        std::ofstream out3("gmesh2.txt");
        gmesh->Print(out3);
    }

    auto mfmesh = GenerateMultiphysicsMesh(gmesh, porder, porderlow);
    GroupElements(mfmesh,porderlow);
    if(1)
    {
        std::ofstream out("MFMeshCondensed.txt");
        mfmesh->Print(out);
    }
    TPZLinearAnalysis analysis(mfmesh);
    TPZFStructMatrix<STATE> matskl(mfmesh);
    matskl.SetNumThreads(0);
    analysis.SetStructuralMatrix(matskl);

    /// Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt); // ELU //ECholesky // ELDLt
    analysis.SetSolver(step);
    analysis.Assemble();
    TPZMatrixSolver<STATE> *matsolver = dynamic_cast<TPZMatrixSolver<STATE> *>(analysis.Solver());
    TPZFMatrix<STATE> *mat = dynamic_cast<TPZFMatrix<STATE> *>(matsolver->Matrix().operator->());
    TPZFMatrix<STATE> matcopy(*mat),rgb;
    GetRigidBodyModes(mfmesh, rgb);
    TPZFMatrix<STATE> force(rgb.Rows(),rgb.Cols(),0.);
    mat->Multiply(rgb, force);
    force.Print(std::cout);
    mat->AddContribution(0, 0, rgb, 0, rgb, 1, 1.);
    matcopy = *mat;
    TPZManVector<std::complex<STATE> > eigenvalues(mat->Rows(),0.);
    TPZFNMatrix<100,std::complex<STATE>> eigenvectors(mat->Rows(), mat->Rows(),0.);
    matcopy.SolveEigenProblem(eigenvalues, eigenvectors);
    std::cout << eigenvalues << std::endl;
    if(1)
    {
        int64_t r = mat->Rows();
        TPZManVector<std::string> fields = {"Displacement", "SigmaX", "SigmaY", "TauXY"};
        TPZVTKGenerator vtk(mfmesh, fields, "Solution.vtk", 0);
        vtk.SetNThreads(0);
        for(int i = 0; i<r; i++) {
            if(norm(eigenvalues[i]) >= 1.e-10) continue;
            TPZFNMatrix<100,STATE> vec(r,1,0.), lambdavec(r,1,0.);
            for(int j=0; j<r; j++) vec(j,0) = eigenvectors(j,i).real();
//            vec.Print(std::cout);
            mfmesh->LoadSolution(vec);
            mfmesh->TransferMultiphysicsSolution();
//            mfmesh->Solution().Print("sol");
            vtk.Do();

        }
    }
    return;
}
