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



int main(int argc, char *argv[]) {
    gAnalytic.fProblemType = TElasticity2DAnalytic::Etest2;
    gAnalytic.gE = 1.;
    gAnalytic.gPoisson = 0.3;
    int porder = 3;
    int porderlow = 1;
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    //    TPZMaterial::gBigNumber = 1.e16;
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
#ifdef MACOSX
    std::vector<std::string> filenames = {"../honeycomb.1.msh","../honeycomb.2.msh","../honeycomb.4.msh","../honeycomb.8.msh"};
#else
    TPZGeoMesh *gmeshorig = Square.GeometricGmshMesh("Triangle.msh");
#endif
    std::cout << filenames[0] << " Internal order " << porder << " Boundary order " << (porderlow == -1 ? porder : porderlow) << std::endl;
    for(int iref = 0; iref<filenames.size(); iref++) {
        TPZGmshReader Square;
        Square.SetVerbose(0);
        Square.GetDimNamePhysical()[2]["domain"] = matID;
        Square.GetDimNamePhysical()[1]["hexagon"] = matHEXAGON;
        TPZGeoMesh *gmeshorig = Square.GeometricGmshMesh(filenames[iref]);
        AddBCHoneycomb(gmeshorig,matBCD, matBCN);
        if(0)
        {
            std::ofstream out("gmesh.txt");
            gmeshorig->Print(out);
        }
        {
            std::stringstream sout;
            sout << "hexagon." << iref+1 << ".vtk";
            std::ofstream out(sout.str());
            TPZVTKGeoMesh::PrintGMeshVTK(gmeshorig, out, true);
        }

        TPZGeoMesh *gmesh = gmeshorig;
        TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
        AddWrapElements(gmesh);
        if(porderlow >= 0) {
            AddWrapElementsTripleHybrid(gmesh);
        }
        if(0)
        {
            std::ofstream out("Mercier.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
            std::ofstream out3("gmesh.txt");
            gmesh->Print(out3);
        }
        auto mfmesh = GenerateMultiphysicsMesh(gmesh, porder, porderlow);
        
        mfmesh->ComputeNodElCon();
        if(0)
        {
            std::ofstream out("MFMesh.txt");
            mfmesh->Print(out);
        }
        GroupElements(mfmesh, porderlow);
        if(0)
        {
            std::ofstream out("MFMeshCondensed.txt");
            mfmesh->Print(out);
        }
        SolveProblem(mfmesh);
        int64_t neq_cond = mfmesh->NEquations();
        int64_t neq = mfmesh->Solution().Rows();
        std::cout << "Global system neq " << neq_cond << " Total equations " << neq << std::endl;

        TPZManVector<STATE> errors(5,0.);
        mfmesh->ElementSolution().Redim(mfmesh->NElements(), 5);
        mfmesh->EvaluateError(1, errors);
        auto *mat = mfmesh->FindMaterial(1);
        auto matmix = dynamic_cast<TPZMixedSymElasticityND *>(mat);
        auto names = matmix->ErrorNames();
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


