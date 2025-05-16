#include <fstream>
#include "pzlog.h"
#include <iostream>
#include <string>
#include <TPZHDivApproxCreator.h>
#include <pzgmesh.h>
#include <TPZMultiphysicsCompMesh.h>
#include <Material/Elasticity/TPZMixedElasticityND.h>
#include <TPZVTKGenerator.h>
#include <TPZAnalyticSolution.h>
#include <TPZGmshReader.h>
#include "ProblemData.h"
#include "dirs_config.h"
#include <TPZLinearAnalysis.h>
#include <pzstrmatrixot.h>
#include <pzskylstrmatrix.h>
#ifdef PZ_USING_MKL
#include <TPZSSpStructMatrix.h>
#endif
#include <pzstepsolver.h>

const std::unordered_map<std::string, TElasticity2DAnalytic::EDefState> ProblemsMap2D = {
    {"none", TElasticity2DAnalytic::ENone},
    {"dispx", TElasticity2DAnalytic::EDispx},
    {"dispy", TElasticity2DAnalytic::EDispy},
    {"stretchx", TElasticity2DAnalytic::EStretchx},
    {"stretchy", TElasticity2DAnalytic::EStretchy},
    {"shear", TElasticity2DAnalytic::EShear},
    {"bend", TElasticity2DAnalytic::EBend},
    {"loadedbeam", TElasticity2DAnalytic::ELoadedBeam},
    {"test1", TElasticity2DAnalytic::Etest1},
    {"test2", TElasticity2DAnalytic::Etest2}
};

const std::unordered_map<std::string, TElasticity3DAnalytic::EDefState> ProblemsMap3D = {
    {"none", TElasticity3DAnalytic::ENone}
};

const int numthreads = 0;

TPZGeoMesh *ReadMeshFromGmsh(std::string file, ProblemData &simData);

int main(int argc, char *argv[]) {

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    ProblemData simData;
    simData.ReadJson("shear-beam.json");

    const int dim = simData.Dimension;
    REAL E = simData.Domains[0].E;
    REAL nu = simData.Domains[0].nu;
    int PlaneStress = simData.PlaneStress;

    TPZAnalyticSolution *analyticSol = 0;

    switch (dim) {
        case 2: {
            TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
            elas->gE = E;
            elas->gPoisson = nu;
            TElasticity2DAnalytic::EDefState problemType = ProblemsMap2D.at(simData.AnalyticSolution);
            elas->fProblemType = problemType;
            elas->fPlaneStress = PlaneStress;
            analyticSol = elas;
            break;
        }
        case 3: {
            TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
            elas->fE = E;
            elas->fPoisson = nu;
            TElasticity3DAnalytic::EDefState problemType = ProblemsMap3D.at(simData.AnalyticSolution);
            elas->fProblemType = problemType;
            analyticSol = elas;
            break;
        }
        default:
            DebugStop();
    }
#ifdef PZDEBUG
    if (!analyticSol) {
        std::cout << "Error: Analytic solution not found." << std::endl;
        DebugStop();
    }
#endif

    TPZGeoMesh *gmesh = ReadMeshFromGmsh(simData.MeshName, simData);

    TPZHDivApproxCreator approx(gmesh);
    approx.ProbType() = ProblemType::EElastic;
    approx.SetDefaultOrder(simData.StressOrder);
    approx.SetExtraInternalOrder(0);
    approx.HdivFamily() = HDivFamily::EHDivStandard;
    approx.SetShouldCondense(simData.ShouldCondense);
    approx.IsRigidBodySpaces() = simData.FiveSpaces;
    approx.HybridType() = HybridizationType::ENone;

    for (auto &domain : simData.Domains) {
        TPZMixedElasticityND* matdomain = new TPZMixedElasticityND(domain.matID, domain.E, domain.nu, 0., 0., PlaneStress, dim);
        if (analyticSol) {
            matdomain->SetExactSol(analyticSol->ExactSolution(),4);
            matdomain->SetForcingFunction(analyticSol->ForceFunc(), 4);
        }
        approx.InsertMaterialObject(matdomain);

        TPZFMatrix<STATE> val1(3, 3, 0.);
        TPZManVector<STATE> val2(3, 0.);
        for (auto &bc : simData.Boundaries) {
            if (bc.domainID != domain.matID) continue;
            val2 = bc.value;
            TPZBndCond *matbc = matdomain->CreateBC(matdomain, bc.matID, bc.type, val1, val2);
            auto matbcT = dynamic_cast<TPZBndCondT<STATE> *>(matbc);
            if (analyticSol) matbcT->SetForcingFunctionBC(analyticSol->ExactSolution(), 4);
            approx.InsertMaterialObject(matbc);
        }
    }

    TPZCompMesh *cmesh = approx.CreateApproximationSpace();

    TPZLinearAnalysis an(cmesh, RenumType::EMetis);

#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE,TPZStructMatrixOT<STATE>> matstr(cmesh);
#else
    TPZSkylineStructMatrix<STATE> matstr(cmesh);
#endif
    matstr.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matstr);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);

    std::cout << "Assemble matrix with NDoF = " << cmesh->NEquations() << "." << std::endl;
    an.Assemble();
    std::cout << "Assemble finished." << std::endl;

    an.Solve();
    std::cout << "Solve finished." << std::endl;

    std::string vtkfile = simData.AnalyticSolution + ".vtk";

    TPZStack<std::string> fieldnames;
    fieldnames.Push("SigmaX");
    fieldnames.Push("SigmaY");
    fieldnames.Push("TauXY");
    if (dim == 3)
    {
        fieldnames.Push("SigmaZ");
        fieldnames.Push("TauXZ");
        fieldnames.Push("TauYZ");
    }
    fieldnames.Push("Displacement");

    TPZVTKGenerator vtk(cmesh, fieldnames, vtkfile, simData.VTKRes, dim);
    vtk.SetNThreads(numthreads);
    vtk.Do();

    return 0;
}

TPZGeoMesh *ReadMeshFromGmsh(std::string file, ProblemData &simData)
{
    // read mesh from gmsh
    std::string path = std::string(INPUTDIR) + "/" + file;
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    TPZGmshReader reader;
    reader.GeometricGmshMesh(path, gmesh);
    gmesh->BuildConnectivity();
    return gmesh;
}