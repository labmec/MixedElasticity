/*
  This target aims to be an example of usage of the current PZ version
*/
#include <TPZGeoMeshTools.h>
#include "TPZAnalyticSolution.h"
#include <TPZGmshReader.h>
#include "TPZCompMeshTools.h"
#include "pzlog.h"
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZHDivApproxCreator.h"

#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
#include "TPZTimer.h"
#include "TPZSimpleTimer.h"
#include "TPZVTKGenerator.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "Elasticity/TPZMixedElasticityND.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "TPZVTKGeoMesh.h"

enum EMatid  {ENone, EDomain, EBoundary};

/**
   @brief Creates a geometric mesh with elements of a given type on a unit square or cube (depending on the mesh dimension).
   @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
   @param[in] volId Material identifier for the volumetric region.
   @param[in] bcId Material identifier for the boundary.
*/
template<class tshape>
TPZGeoMesh*
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, EMatid bcId);

// The test function
template<class tshape>
void SolveFEMProblem(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily);

void InsertMaterials(int &dim, TPZHDivApproxCreator& hdivc, TPZAnalyticSolution *fAn);

int main() {

    const int xdiv = 10; //Number of elements in each direction
    const int pOrder = 1; // Polynomial degree
    // Family of HDiv approximation spaces. 
    // The possible choices are HDivFamily::EHDivStandard, HDivFamily::EHDivConstant and HDivFamily::EHDivKernel
    HDivFamily hdivfam = HDivFamily::EHDivConstant; 
    
    //Creates the geometric mesh for the given topology and solve the FEM problem.
    SolveFEMProblem<pzshape::TPZShapeQuad>(xdiv,pOrder,hdivfam);
    //   SolveFEMProblem<pzshape::TPZShapeTriangle>(xdiv,pOrder,hdivfam);
    //   SolveFEMProblem<pzshape::TPZShapeTetra>(xdiv,pOrder,hdivfam);
    //   SolveFEMProblem<pzshape::TPZShapeCube>(xdiv,pOrder,hdivfam);
   
    return 0;
}


template<class tshape>
void SolveFEMProblem(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily)
{

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    int DIM = tshape::Dimension;
    TPZVec<int> nDivs;

    if (DIM == 2) nDivs = {xdiv,xdiv};
    if (DIM == 3) nDivs = {xdiv,xdiv,xdiv};
    
    // Creates/import a geometric mesh  
    auto gmesh = CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary);

    // Creates an hdivApproxCreator object. It is an environment developped to
    // help creating H(div)-family possible approximation spaces.
    TPZHDivApproxCreator hdivCreator(gmesh);
    //Set the family of H(div) functions: Standard, Constant or Kernel
    hdivCreator.HdivFamily() = hdivfamily; 
    //Set the problem type to be solved: Only EDarcy and EElastic are currently available
    hdivCreator.ProbType() = ProblemType::EDarcy;
    //Includes the rigid body spaces (constant flux and pressure) if set as true
    hdivCreator.IsRigidBodySpaces() = false;
    //Set the default polynomial order
    hdivCreator.SetDefaultOrder(pOrder);
    //Set the extra polynomial order for the bubble functions. If zero, the polynomial degree
    //of the internal functions are the same as the default order
    hdivCreator.SetExtraInternalOrder(0);
    //Sets if the resulting problem should or not be condensed 
    hdivCreator.SetShouldCondense(false);
    // hdivCreator.SetShouldCondense(false);
    
    //Sets the type of hybridizantion desired.
    //The current options are HybridizationType::ENone, HybridizationType::EStandard 
    //and HybridizationType::ESemi (the last only works with H(div)-constant spaces) 
    hdivCreator.HybridType() = HybridizationType::ENone;
    // hdivCreator.HybridType() = HybridizationType::EStandard;

    // Prints gmesh mesh properties
    std::string vtk_name = "geoMesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

    //Creates an analytical solution to test
    TPZAnalyticSolution *gAnalytic = 0;
    if (hdivCreator.ProbType() == ProblemType::EDarcy){
        TLaplaceExample1 *lap = new TLaplaceExample1;
        lap->fExact = TLaplaceExample1::EHarmonic;
        gAnalytic = lap;
    } else if (hdivCreator.ProbType() == ProblemType::EElastic){
        if (DIM == 2){
            TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
            elas->gE = 1.e3;
            elas->gPoisson = 0.3;
            elas->fProblemType = TElasticity2DAnalytic::ELoadedBeam;
            elas->fPlaneStress = 0;
            gAnalytic = elas;
        } else if (DIM == 3){
            TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
            elas->fE = 1.;
            elas->fPoisson = 0.0;
            elas->fProblemType = TElasticity3DAnalytic::ELoadedBeam;
            gAnalytic = elas;
        }
    } else {
        DebugStop();
    }

    //Insert Materials
    InsertMaterials(DIM,hdivCreator,gAnalytic);

    //Gets the Multiphysics mesh from the HdivApproxCreator
    TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
    //Here the cmesh is printed
    std::string txt = "cmesh.txt";
    std::ofstream myfile(txt);
    cmesh->Print(myfile);

    //Create the analysis environment
    TPZLinearAnalysis an(cmesh,RenumType::ESloan);
    an.SetExact(gAnalytic->ExactSolution(),4);
    // if (hdivCreator.ProbType() == ProblemType::EDarcy){
        
    // } else if (hdivCreator.ProbType() == ProblemType::EElastic){
    //     an.SetExact(gAnalytic,solOrder);
    // } else {
    //     DebugStop();
    // }

    // Solve problem
    constexpr int nThreads{10};
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> stiffness(cmesh);   
    stiffness.SetNumThreads(nThreads);
    an.SetStructuralMatrix(stiffness);

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //Assemble and solve the problem.
    an.Run();
  
    //Printing results in a vtk file
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), an.Mesh());
    const std::string plotfile = "myfile";//sem o .vtk no final
    constexpr int vtkRes{0};//Resolution
    
    //Fields to be printed
    TPZVec<std::string> fields;
    if (hdivCreator.ProbType() == ProblemType::EDarcy){
        fields = {"Pressure","ExactPressure","Flux","ExactFlux"};
    } else if (hdivCreator.ProbType() == ProblemType::EElastic){
        fields = {"Displacement","SigmaX","SigmaY","TauXY"};
    } else {
        DebugStop();
    }
    
    auto vtk = TPZVTKGenerator( an.Mesh(), fields, plotfile, vtkRes);
    vtk.Do();

    //Compute error
    std::ofstream anPostProcessFile("postprocess.txt");
    TPZManVector<REAL,5> error;
    int64_t nelem = cmesh->NElements();
    cmesh->LoadSolution(cmesh->Solution());
    cmesh->ExpandSolution();
    cmesh->ElementSolution().Redim(nelem, 5);
    an.PostProcessError(error,false,anPostProcessFile);
    
    //Print Errors
    std::cout << "ERROR[0] = " << std::scientific << std::setprecision(15) << error[0] << std::endl;
    std::cout << "ERROR[1] = " << error[1] << std::endl;
    std::cout << "ERROR[2] = " << error[2] << std::endl;
    std::cout << "ERROR[3] = " << error[3] << std::endl;
    std::cout << "ERROR[4] = " << error[4] << std::endl;
    
}

//Create 
template <class tshape>
TPZGeoMesh*
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, EMatid bcId)
{
    
    MMeshType meshType;
    int dim = tshape::Dimension;

    switch (tshape::Type())
    {
    case ETriangle:
        meshType = MMeshType::ETriangular;
        break;
    case EQuadrilateral:
        meshType = MMeshType::EQuadrilateral;
        break;
    case ETetraedro:
        meshType = MMeshType::ETetrahedral;
        break;
    case ECube:
        meshType = MMeshType::EHexahedral;
        break;
        case EPrisma:
        meshType = MMeshType::EPrismatic;
        break;
    default:
        DebugStop();
    }

    TPZManVector<REAL,3> minX = {0,0,0};
    TPZManVector<REAL,3> maxX = {1,1,1};
    int nMats = 2*dim+1;

    //all bcs share the same id
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,bcId);
    matIds[0] = volId;
    // matIds[1] = bcId;
    // matIds[2] = EBoundary1;
    // matIds[3] = EBoundary1;
    // matIds[4] = EBoundary1;
    
    TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                        matIds, nDivs, meshType,createBoundEls);
    // TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(meshType,
    //                     volId,createBoundEls, bcId);
    
    return gmesh;
    
}

void InsertMaterials(int &dim, TPZHDivApproxCreator& hdivCreator,TPZAnalyticSolution *fAn){

    if (hdivCreator.ProbType() == ProblemType::EDarcy){
        TPZMixedDarcyFlow* matdarcy = new TPZMixedDarcyFlow(EDomain,dim);
        // matdarcy->SetConstantPermeability(1.);
        TLaplaceExample1* lapl = dynamic_cast<TLaplaceExample1*> (fAn) ;
        matdarcy->SetExactSol(lapl->ExactSolution(),4);

        hdivCreator.InsertMaterialObject(matdarcy);

        TPZFMatrix<STATE> val1(1,1,0.);
        TPZManVector<STATE> val2(1,0.);
        TPZBndCondT<STATE> *BCond1 = matdarcy->CreateBC(matdarcy, EBoundary, 0, val1, val2);
        BCond1->SetForcingFunctionBC(lapl->ExactSolution(),4);
        hdivCreator.InsertMaterialObject(BCond1);
    } else if (hdivCreator.ProbType() == ProblemType::EElastic){
        TElasticity2DAnalytic *elas2D;
        TElasticity3DAnalytic *elas3D;
        TPZMixedElasticityND* matelas;
        if (dim == 2) {
            elas2D = dynamic_cast<TElasticity2DAnalytic*> (fAn) ;
            matelas = new TPZMixedElasticityND(EDomain, elas2D->gE, elas2D->gPoisson, 0, 0, elas2D->fPlaneStress, dim);
            matelas->SetExactSol(elas2D->ExactSolution(),4);
            hdivCreator.InsertMaterialObject(matelas);

            TPZFMatrix<STATE> val1(dim,dim,0.);
            TPZManVector<STATE> val2(dim,0.);
            TPZBndCondT<STATE> *BCond1 = matelas->CreateBC(matelas, EBoundary, 0, val1, val2);
            BCond1->SetForcingFunctionBC(elas2D->ExactSolution(),4);
            hdivCreator.InsertMaterialObject(BCond1);
        }
        if (dim == 3) {
            elas3D = dynamic_cast<TElasticity3DAnalytic*> (fAn) ;
            matelas = new TPZMixedElasticityND(EDomain, elas3D->fE, elas3D->fPoisson, 0, 0, 0, dim);
            matelas->SetExactSol(elas3D->ExactSolution(),4);

            hdivCreator.InsertMaterialObject(matelas);

            TPZFMatrix<STATE> val1(dim,dim,0.);
            TPZManVector<STATE> val2(dim,0.);
            TPZBndCondT<STATE> *BCond1 = matelas->CreateBC(matelas, EBoundary, 0, val1, val2);
            BCond1->SetForcingFunctionBC(elas3D->ExactSolution(),4);
            hdivCreator.InsertMaterialObject(BCond1);
        }
        
    } else {
        DebugStop();//Material Not Implemented
    }
}



