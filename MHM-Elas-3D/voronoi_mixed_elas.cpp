#include <TPZGeoMeshTools.h>
#include "TPZAnalyticSolution.h"
#include <TPZGmshReader.h>
#include "TPZCompMeshTools.h"
#include "pzlog.h"

#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
#include "TPZTimer.h"
#include "fstream"
#include "TPZSimpleTimer.h"
#include "TPZVTKGenerator.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZAnalyticSolution.h"
#include "TPZHDivApproxCreator.h"
#include "Elasticity/TPZMixedElasticityND.h"
#include "TPZVTKGeoMesh.h"
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"

enum EMatid  {ENone, EDomain, EBoundary, EPont, EWrap, EIntface, EPressureHyb, EBCLeft, EBCRight, EZeroNeu};

/**
 @brief Creates a geometric mesh with elements of a given type on a unit square or cube (depending on the mesh dimension).
 @param[in] meshType element type to be created.
 @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
 @param[in] volId Material identifier for the volumetric region.
 @param[in] bcId Material identifier for the boundary.
 */
template<class tshape>
TPZGeoMesh*
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, EMatid bcId);

template<class tshape>
TPZGeoMesh*
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, EMatid bcLeft, EMatid bcRight, EMatid bcZeroNeumann);

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

/**
 @brief Reads the test mesh from gmsh
 @param[in] file_name the .msh mesh file.
 */
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

void FixBCsForHomogeneousTest(TPZGeoMesh* gmesh);

int main()
{
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif
  
  const int pord = 1;
  int DIM = 3;
  TPZVec<int> nDivs;
  
  const int ndiv = 3;
  if (DIM == 2) nDivs = {ndiv,ndiv};
  if (DIM == 3) nDivs = {ndiv,ndiv,ndiv};
  
  // Creates/import a geometric mesh
  TPZGeoMesh* gmesh = nullptr;
  const bool readFromGMesh = true;
  if(readFromGMesh){
    gmesh = ReadMeshFromGmsh("MHMesh1.msh");
    FixBCsForHomogeneousTest(gmesh);
  }
  else{
    if(DIM == 2)
      gmesh = CreateGeoMesh<pzshape::TPZShapeQuad>(nDivs, EDomain, EBCLeft, EBCRight, EZeroNeu);
    else
      gmesh = CreateGeoMesh<pzshape::TPZShapeTetra>(nDivs, EDomain, EBCLeft, EBCRight, EZeroNeu);
  }
  
  std::ofstream outgmesh("geomesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outgmesh);
  
  TPZHDivApproxCreator hdivCreator(gmesh);
  hdivCreator.HdivFamily() = HDivFamily::EHDivStandard;
  hdivCreator.ProbType() = ProblemType::EElastic;
  hdivCreator.IsRigidBodySpaces() = false;
  hdivCreator.SetDefaultOrder(pord);
  hdivCreator.SetExtraInternalOrder(0);
  hdivCreator.SetShouldCondense(false);
//  hdivCreator.HybridType() = HybridizationType::EStandard;
  hdivCreator.HybridType() = HybridizationType::ENone;
  
  TPZAnalyticSolution *gAnalytic = 0;
  TPZMixedElasticityND* matelastic = 0;
  if(DIM == 2)
  {
    TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
    elas->gE = 1.e3;
    elas->gPoisson = 0.0;
//    elas->fProblemType = TElasticity2DAnalytic::EThiago;
    elas->fPlaneStress = 0;
//    gAnalytic = elas;
    matelastic = new TPZMixedElasticityND(EDomain, elas->gE, elas->gPoisson, 0, 0, elas->fPlaneStress, DIM);
    
  }
  else if(DIM == 3)
  {
    TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
    elas->fE = 1.;//206.8150271873455;
    elas->fPoisson = 0.0;//0.3040039545229857;
//    elas->fProblemType = TElasticity3DAnalytic::EDispx;
//    gAnalytic = elas;
    matelastic = new TPZMixedElasticityND(EDomain, elas->fE, elas->fPoisson, 0, 0, 0 /*planestress*/, DIM);
  }
  
  
  //Insert Materials
//  matelastic->SetExactSol(gAnalytic->ExactSolution(),4);
  
  hdivCreator.InsertMaterialObject(matelastic);
  
  const int diri = 0, neu = 1;
  // Fixed on the left
  TPZFMatrix<STATE> val1(DIM,DIM,0.);
  TPZManVector<STATE> val2(DIM,0.);
  val2[0] = 0.;
  TPZBndCondT<STATE> *BCond1 = matelastic->CreateBC(matelastic, EBCLeft, diri, val1, val2);
  hdivCreator.InsertMaterialObject(BCond1);
  
  // pull on right
  val2[0] = 2.;
  TPZBndCondT<STATE> *BCond2 = matelastic->CreateBC(matelastic, EBCRight, neu, val1, val2);
  hdivCreator.InsertMaterialObject(BCond2);

  // zero neumann on rest
  val2[0] = 0.;
  TPZBndCondT<STATE> *BCond3 = matelastic->CreateBC(matelastic, EZeroNeu, neu, val1, val2);
  hdivCreator.InsertMaterialObject(BCond3);
  
  //Multiphysics mesh
  TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
//  std::string txt = "cmesh.txt";
//  std::ofstream myfile(txt);
//  cmesh->Print(myfile);
  
  //Create analysis environment
  TPZLinearAnalysis an(cmesh);
  
  //Solve problem
  SolveProblemDirect(an,cmesh);
   
  {
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), cmesh);
    TPZSimpleTimer postProc("Post processing2");
    const std::string plotfile = "solution";
    constexpr int vtkRes{0};
    
    
    TPZVec<std::string> fields = {
      // "ExactDisplacement",
      // "ExactStress",
      "Displacement",
      "SigmaX",
      "SigmaY",
      "TauXY"
    };
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    
    vtk.Do();
  }
  
  return 0;
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
  
  TPZManVector<REAL,3> minX = {-1,-1,-1};
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

template <class tshape>
TPZGeoMesh*
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, EMatid bcLeft, EMatid bcRight, EMatid bcZeroNeumann)
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
  
  TPZManVector<REAL,3> minX = {-1,-1,-1};
  TPZManVector<REAL,3> maxX = {1,1,1};
  int nMats = 2*dim+1;
  
  //all bcs share the same id
  constexpr bool createBoundEls{true};
  TPZVec<int> matIds(nMats,bcZeroNeumann);
  matIds[0] = volId;
  if(dim == 2){
    matIds[4] = bcLeft;
    matIds[2] = bcRight;
  }
  else {
    matIds[2] = bcLeft;
    matIds[4] = bcRight;
  }
  
  TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                                                           matIds, nDivs, meshType,createBoundEls);
  return gmesh;
  
}


template <class tshape>
TPZGeoMesh*
ReadMeshFromGmsh(std::string file_name)
{
  //read mesh from gmsh
  TPZGeoMesh *gmesh;
  gmesh = new TPZGeoMesh();
  {
    TPZGmshReader reader;
    // essa interface permite voce mapear os nomes dos physical groups para
    // o matid que voce mesmo escolher
    TPZManVector<std::map<std::string,int>,4> stringtoint(4);
    stringtoint[3]["Domain"] = 1;
    stringtoint[2]["Surfaces"] = 2;
    
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(file_name,gmesh);
  }
  
  return gmesh;
}


void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    //sets number of threads to be used by the solver
    constexpr int nThreads{16};
    // TPZSkylineStructMatrix<REAL> matskl(cmesh);
    // TPZSSpStructMatrix<STATE> matskl(cmesh);
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matskl(cmesh);
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    
    an.SetSolver(step);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    an.Assemble();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time Assemble = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    ///solves the system
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    an.Solve();
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time Solve = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count() << "[ms]" << std::endl;

}

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name)
{
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(4);
        stringtoint[3]["volume"] = EDomain;
        stringtoint[2]["bound"] = EZeroNeu;

        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh,false);
    }

    return gmesh;
}

void FixBCsForHomogeneousTest(TPZGeoMesh* gmesh){
  const int64_t nel = gmesh->NElements();
  for (int iel = 0; iel < nel; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if(gel->MaterialId() != EZeroNeu) continue;
    
    TPZManVector<REAL,3> centqsi(2,0.),cent(3,0.);
    gel->CenterPoint(gel->NSides()-1, centqsi);
    gel->X(centqsi, cent);
    if(fabs(cent[0]) < 1.e-10){
      gel->SetMaterialId(EBCLeft);
    }
    if(fabs(cent[0]-1.) < 1.e-10){
      gel->SetMaterialId(EBCRight);
    }
  }
}
