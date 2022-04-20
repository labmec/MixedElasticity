#include "TPZMixedElasticityUtils.h"
#include "pzstepsolver.h"
#include "TPZStructMatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPZSSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzcmesh.h"
#include "TPZMatErrorCombinedSpaces.h"


void TPZMixedElasticityUtils::SolveProblem(TPZLinearAnalysis &an, int fDomainMatID){

    int numthreads = 0;

    TPZCompMesh * cmesh = an.Mesh();
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>>  matskl(cmesh);

    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);

    std::cout << "Assemble matrix with NDoF = " << cmesh->NEquations() << "." << std::endl;
    an.Assemble(); //Assembles the global stiffness matrix (and load vector)
    std::cout << "Assemble finished." << std::endl;

    // TPZMatErrorCombinedSpaces<STATE> *materr = dynamic_cast<TPZMatErrorCombinedSpaces<STATE> *>(cmesh_m_HDiv->FindMaterial(fDomainMatID));
    // TPZManVector<REAL, 6> Errors(materr->NEvalErrors());
    // TElasticity2DAnalytic example;
    // an.SetExact(example.Exact(),3);
    //            an.PostProcessError(Errors,std::cout);

#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    if (false) {
        std::ofstream filestiff("stiffness.nb");
        an.MatrixSolver<STATE>().Matrix()->Print("K1 = ", filestiff, EMathematicaInput);

        std::ofstream filerhs("rhs.nb");
        an.Rhs().Print("R = ", filerhs, EMathematicaInput);
    }
#endif

    std::cout << "Solving." << std::endl;
    an.Solve();
    std::cout << "Solved." << std::endl;


    {
        TPZStepSolver<STATE> solver;
        an.SetSolver(solver);
    }
#ifdef PZDEBUG
    if (0) {
        std::ofstream file("file.txt");
        an.Solution().Print("sol=", file, EMathematicaInput);

    }
#endif


}


void TPZMixedElasticityUtils::PrintResults(TPZLinearAnalysis &an){

    std::string plotfile;
    {
        std::stringstream sout;
        sout << "results.vtk";
        plotfile = sout.str();
    }
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("SigmaX");
    scalnames.Push("SigmaY");
    scalnames.Push("TauXY");
    // scalnames.Push("SigmaR");
    // scalnames.Push("SigmaT");
    // scalnames.Push("TauRT");
    // vecnames.Push("Flux");
    vecnames.Push("displacement");
    vecnames.Push("ExactDisplacement");
    vecnames.Push("ExactStress");
    vecnames.Push("StressAiry");
    vecnames.Push("StressTotal");
    // vecnames.Push("Stress");
    int count = 0;
    an.SetStep(count);
    an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
    an.PostProcess(3);

}

