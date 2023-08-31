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
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZAnalyticSolution.h"
#include "TPZHDivApproxCreator.h"
#include "Elasticity/TPZMixedElasticityND.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZVTKGeoMesh.h"
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzsbstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "pzsubcmesh.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "meshpath_config.h"
#include "pzintel.h"
#include <algorithm>

#include "VoronoiAnalise.h"

enum EMatid  {ENone, EDomain, EBoundary, EMHM, EPont, EWrap, EIntface, EPressureHyb, EBCLeft, EBCRight, EZeroNeu};

const int global_nthread = 16;

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

void CheckBCs(TPZGeoMesh* gmesh);

void DivideGMesh(TPZGeoMesh *gmesh, int internaldiv, int skeletondiv);

void IdentifyMHMDomain(TPZGeoMesh *gmesh, TPZVec<int> &domain);

void Substructure(TPZCompMesh *cmesh, TPZVec<int> &domains);

void GroupElements(TPZCompMesh *submesh);

void CondenseElements(TPZCompMesh *submesh);

// find the element groups which link to identical subdomain
void IdentifyInteractionPlanes(TPZMultiphysicsCompMesh *cmesh, TPZVec<int> &domain, const int pord, REAL& smallestVoronoiInterfaceArea);

class InputParser{
    public:
        InputParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }
        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }
        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
};

// argv parameters:

// np: number of points used to create the voronoi mesh (1 to 5)
// analysol: type of analytical sol used for error computation (0 = EStretchX, 1 = EStretchY, 2 = EYotov)
// pord: polynomial order of the approximation in general
// pordskel: polynomial order of the skeleton
// nrefint: number of internal refinements
// useredspaceface: using reduced space on voronoi interface (true of false)
// pordface: polynomial order of function defined on voronoi interface
int main(int argc, char *argv[])
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG(std::string(MESHES_DIR) + "/" + "log4cxx.cfg");
#endif
    std::ofstream out("voronoi-results.txt",std::ios::app);
    
    TElasticity3DAnalytic::EDefState asol;
    bool useReducedSpaceOnFace = true;
    int pordface = 0, pord = 1, pordskel = 1;
    int nrefinternal = 0;
    std::string voronoi_np = "1";
//    if (argc > 1 && argc != 15) DebugStop();
    
    if (argc > 1) {
        InputParser input(argc, argv);
//        if(!input.cmdOptionExists("-np") || !input.cmdOptionExists("-analysol") || !input.cmdOptionExists("-pord") ||
//           !input.cmdOptionExists("-pordskel") || !input.cmdOptionExists("-nrefint")){
//            DebugStop();
//        }
        voronoi_np = input.getCmdOption("-np");
        const int input_asol = stoi(input.getCmdOption("-analysol"));
        
        if (input_asol == 0)
            asol = TElasticity3DAnalytic::EDispx;
        else if(input_asol == 1)
            asol = TElasticity3DAnalytic::EDispy;
        else if(input_asol == 2)
            asol = TElasticity3DAnalytic::EYotov;
        else
            DebugStop();
        
        if(input.cmdOptionExists("-useredspaceface")) useReducedSpaceOnFace = stoi(input.getCmdOption("-useredspaceface"));
        if(input.cmdOptionExists("-pordface")) pordface = stoi(input.getCmdOption("-pordface"));
        
        if(input.cmdOptionExists("-pord")) pord = stoi(input.getCmdOption("-pord"));
        if(input.cmdOptionExists("-pordskel")) pordskel = stoi(input.getCmdOption("-pordskel"));
        if(input.cmdOptionExists("-nrefint")) nrefinternal = stoi(input.getCmdOption("-nrefint"));
        
        out << "\n----------------- Starting new simulation -----------------" << std::endl;
        out << "MHMeshEquiTet_np" << voronoi_np << " | AnalySol = " << input_asol << " | useReducedSpaceOnFace = " << useReducedSpaceOnFace << " | pordface = " << pordface << " | pord = " << pord  << " | pordskel = " << pordskel << " | nrefint = " << nrefinternal << std::endl;
        std::cout << "\n----------------- Starting new simulation -----------------" << std::endl;
        std::cout << "MHMeshEquiTet_np" << voronoi_np << " | AnalySol = " << input_asol << " | useReducedSpaceOnFace = " << useReducedSpaceOnFace << " | pordface = " << pordface << " | pord = " << pord << " | pordskel = " << pordskel  << " | nrefint = " << nrefinternal << std::endl;
                            

    }

        
    
    std::string meshname;
    if(argc > 1){
//        meshname = "MHMesh_np" + voronoi_np + ".msh";
        meshname = "MHMeshEquiTet_np" + voronoi_np + ".msh";
    }
    else{
//        meshname = "MHMesh_np2.msh";
        meshname = "MHMeshEquiTet_np1.msh";
    }
        
    int nrefskel = 0;
            
    // Creates/import a geometric mesh
    TPZGeoMesh* gmesh = nullptr;
    const bool readFromGMesh = true;
    if(readFromGMesh){
        meshname = std::string(MESHES_DIR) + "/" + meshname;
        gmesh = ReadMeshFromGmsh(meshname);
    }
    else{
        TPZVec<int> nDivs;
        const int ndiv = 1;
        nDivs = {ndiv,ndiv,ndiv};
//        gmesh = CreateGeoMesh<pzshape::TPZShapeTetra>(nDivs, EDomain, EBCLeft, EBCRight, EZeroNeu);
        gmesh = CreateGeoMesh<pzshape::TPZShapeTetra>(nDivs, EDomain, EZeroNeu);
    }
    int DIM = gmesh->Dimension();
    if (nrefinternal || nrefskel) {
        if(nrefskel > nrefinternal) DebugStop();
        DivideGMesh(gmesh, nrefinternal, nrefskel);
    }
    CheckBCs(gmesh);

    {
        std::ofstream outgmesh("geomesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outgmesh);
        std::ofstream outgmesh2("geomesh.txt");
        gmesh->Print(outgmesh2);
    }

    TPZVec<int> domain(gmesh->NElements(),-1);
    IdentifyMHMDomain(gmesh, domain);
    {
        std::ofstream out("domains.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, domain, true);
    }
    
    TPZHDivApproxCreator hdivCreator(gmesh);
    hdivCreator.HdivFamily() = HDivFamily::EHDivStandard;
    hdivCreator.ProbType() = ProblemType::EElastic;
    hdivCreator.IsRigidBodySpaces() = true;
    hdivCreator.SetDefaultOrder(pord);
    hdivCreator.SetExtraInternalOrder(0);
    hdivCreator.SetShouldCondense(false);
    //  hdivCreator.HybridType() = HybridizationType::EStandard;
    hdivCreator.HybridType() = HybridizationType::ENone;
    
    TPZAnalyticSolution *gAnalytic = 0;
    TPZMixedElasticityND* matelastic = 0;
    if(DIM == 2) {
        DebugStop(); // Not interested in 2d for now
    }
    else if(DIM == 3) {
        TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
        elas->fE = 250.;//206.8150271873455;
//        elas->fPoisson = 0.25;//0.3040039545229857;
        elas->fPoisson = 0.;
        elas->fProblemType = TElasticity3DAnalytic::EYotov;
//        elas->fProblemType = TElasticity3DAnalytic::EDispx;
        if(argc > 1) elas->fProblemType = asol;
        gAnalytic = elas;
        matelastic = new TPZMixedElasticityND(EDomain, elas->fE, elas->fPoisson, 0, 0, 0 /*planestress*/, DIM);
    }
    
    
    //Insert Materials
    matelastic->SetExactSol(gAnalytic->ExactSolution(),4);
    matelastic->SetForcingFunction(gAnalytic->ForceFunc(), 3);
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
    val2[0] = 1.;
    TPZBndCondT<STATE> *BCond3 = matelastic->CreateBC(matelastic, EZeroNeu, diri, val1, val2);
    BCond3->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
    hdivCreator.InsertMaterialObject(BCond3);
    
    // Interface materials for MHM
    val2[0] = 0.;
    TPZBndCondT<STATE> *MHM = matelastic->CreateBC(matelastic, EMHM, diri, val1, val2);
    hdivCreator.InsertMaterialObject(MHM);
    
    //Multiphysics mesh
//    TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
    // Create by myself
    TPZMultiphysicsCompMesh *cmeshmulti = nullptr;
    {
        std::cout << "\n---------------- Creating Space -----------------" << std::endl;
        hdivCreator.CheckSetupConsistency();
        hdivCreator.SetMeshElementType();

        int lagLevelCounter = 1;
        TPZManVector<TPZCompMesh*,7> meshvec(5);
        hdivCreator.CreateAtomicMeshes(meshvec,lagLevelCounter);
        
        // possibly decrease polynomial order of skeleton
        TPZCompMesh* cmeshflux = meshvec[0];
        if(pord != pordskel){
            if(pordskel > pord) DebugStop(); // cannot happen!
            for (int i = 0; i < cmeshflux->NElements(); i++) {
                TPZCompEl* cel = cmeshflux->Element(i);
                if(!cel) continue;
                TPZGeoEl* gel = cel->Reference();
                if (gel->MaterialId() != EMHM) continue;
                
                TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(cel);
                if(!intel) DebugStop();
                
                intel->PRefine(pordskel);
            }
        }
        cmeshflux->ExpandSolution();
        hdivCreator.CreateMultiPhysicsMesh(meshvec,lagLevelCounter,cmeshmulti);
    }
    TPZMultiphysicsCompMesh *cmesh = cmeshmulti;
    
    
    
    REAL smallestVoronoiInterfaceArea = -1.;
    if(useReducedSpaceOnFace){
        if(pordface > pord) DebugStop();
        IdentifyInteractionPlanes(cmesh, domain, pordface, smallestVoronoiInterfaceArea);
    }
    out << "smallestVoronoiInterfaceArea = " << smallestVoronoiInterfaceArea << std::endl;
    

    
    // Loop over elements. Get the skeleton 2d ones, and change their porder
    
    
    
    Substructure(cmesh, domain);
//    cmesh->LoadReferences();
//    GroupElements(cmesh);
//    CondenseElements(cmesh);
    std::cout << "NEquations " << cmesh->NEquations() << std::endl;
    //Create analysis environment
#ifdef PZDEBUG
    {
        std::string txt = "multiphysicsmesh.txt";
        std::ofstream myfile(txt);
        cmesh->Print(myfile);
    }
#endif
    TPZLinearAnalysis an(cmesh);
    
    //Solve problem
    SolveProblemDirect(an,cmesh);
    if(0)
    {
        auto SetDisp = [](TPZCompMesh *cmesh)
        {
            cmesh->Solution().Zero();
            TPZBlock &block = cmesh->Block();
            int64_t ncon = cmesh->NConnects();
            TPZFMatrix<STATE> &sol = cmesh->Solution();
            for(int64_t ic = 0; ic<ncon; ic++) {
                TPZConnect &c = cmesh->ConnectVec()[ic];
                if(c.LagrangeMultiplier() == 4) {
                    int64_t seq = c.SequenceNumber();
                    int64_t pos = block.Position(seq);
                    sol(pos,0) = 1.;
                }
            }
            cmesh->LoadSolution(sol);
            if(0)
            {
                std::ofstream out("mesh.txt");
                cmesh->Print(out);
            }
            std::cout << "printed\n";
        };
        SetDisp(cmesh);
        if(0) {
            int64_t nel = cmesh->NElements();
            for(int64_t el = 0; el<nel; el++) {
                TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cmesh->Element(el));
                if(sub) {
                    SetDisp(sub);
                }
            }
        }
    }
    cmesh->TransferMultiphysicsSolution();
#ifdef PZDEBUG
    {
        std::string txt = "multiphysicsmesh.txt";
        std::ofstream myfile(txt);
        cmesh->Print(myfile);
    }
#endif
    

    // -------> Calculating error
    std::cout << "------- Starting PostProc Error -------" << std::endl;
    an.SetExact(gAnalytic->ExactSolution());
    an.SetThreadsForError(global_nthread);
    std::ofstream ErroOut("myerrors.txt", std::ios::app);
    TPZMaterial *mat = cmesh->FindMaterial(EDomain);
    TPZMatErrorCombinedSpaces<STATE> *materr = dynamic_cast<TPZMatErrorCombinedSpaces<STATE>*>(mat);
    TPZManVector<REAL, 10> Errors(materr->NEvalErrors());
    bool store_errors = true;
    Errors.Fill(0.);
    cmesh->ElementSolution().Redim(cmesh->NElements(), Errors.size());
    {
        int64_t nel = cmesh->NElements();
        for(int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(cel);
            if(submesh) {
                submesh->ElementSolution().Redim(submesh->NElements(), Errors.size());
            }
        }
    }
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    an.PostProcessError(Errors, store_errors, ErroOut);
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time PostProc Error = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count()/1000. << " s" << std::endl;
    
    std::cout << "Computed errors." << std::endl;
    // error_sigma - error_energy - error_div_sigma - error_u - error_r - error_as - energy_norm_exact_sol
    std::cout.precision(15);
    std::cout.setf(std::ios::fixed);
    std::cout << "Errors = ";
    std::cout << Errors << std::endl;
    out << "Errors = ";
    out << Errors << std::endl;
    out << "ErrorsFixedPrecision = ";
    out.precision(15);
    out.setf(std::ios::fixed);
    out << Errors << std::endl;

    //  for (int i = 0; i < Errors.size(); i++) {
    //    std::cout << Errors[i] << std::endl;
    //  }
//    cmesh->ElementSolution().Print("element error");
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
            "TauXY",
            "ElementSigmaError"
        };
        auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
        vtk.SetNThreads(global_nthread);
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
    constexpr int nThreads{global_nthread};

#if defined(__x86_64__) || defined(__x86_64)
    TPZSSpStructMatrix<STATE> matskl(cmesh);
//    TPZFStructMatrix<STATE> matskl(cmesh);
#elif defined(__arm__) || defined(__aarch64__)
    TPZSkylineStructMatrix<REAL> matskl(cmesh);
//    TPZFStructMatrix<> matskl(cmesh);
//    TPZBandStructMatrix<> matskl(cmesh);
#endif
    
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    
    an.SetSolver(step);
    
    std::cout << "------- Starting Assemble -------" << std::endl;
    std::cout << "Nequations = " << an.Mesh()->NEquations() << std::endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    an.Assemble();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time Assemble = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()/1000. << " s" << std::endl;
    
    ///solves the system
    std::cout << "------- Starting Solve -------" << std::endl;
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    an.Solve();
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time Solve = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count()/1000. << " s" << std::endl;
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
        stringtoint[2]["internal"] = EMHM;
        
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh,false);
    }
    
    return gmesh;
}

void CheckBCs(TPZGeoMesh* gmesh){
    const int64_t nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl* gel = gmesh->Element(iel);
        if(gel->Dimension() != 3) continue;
        int firsts = gel->FirstSide(2);
        for( int s = firsts; s<gel->NSides()-1; s++)
        {
            TPZGeoElSide gelside(gel,s);
            if(gelside.Neighbour() == gelside)
            {
                DebugStop();
            }
        }
    }
}

void DivideGMesh(TPZGeoMesh *gmesh, int internaldiv, int skeletondiv){
    for (int div=0; div<internaldiv; div++) {
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->MaterialId() == EMHM) continue; // == EMHM skips all interfaces between macro domains
            if(gel->HasSubElement()) continue;
            TPZManVector<TPZGeoEl *> subels;
            gel->Divide(subels);
        }
    }
    for (int div=0; div<skeletondiv; div++) {
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->MaterialId() != EMHM) continue; // != EMHM skips everyone that is not interface between macro domains
            if(gel->HasSubElement()) continue;
            TPZManVector<TPZGeoEl *> subels;
            gel->Divide(subels);
        }
    }
}

static std::set<int64_t> process;

void SetSubdomain(TPZGeoEl *gel, int domain, TPZVec<int> &alldomains)
{
    int64_t index = gel->Index();
//    std::cout << "Inserting " << index << " d " << domain << std::endl;
    alldomains[index] = domain;
    if(gel->HasSubElement()) {
        int nsubel = gel->NSubElements();
        if(nsubel){
//            std::cout << "nsubelements " << nsubel << std::endl;
        }
        else {
            std::cout << "I should stop\n";
        }
        for(int is=0; is<nsubel; is++) {
            SetSubdomain(gel->SubElement(is), domain, alldomains);
        }
    }
}


void ProcessElements(TPZGeoMesh *gmesh, TPZVec<int> &domain) {
    while(process.size()) {
        int64_t index = *process.begin();
        process.erase(index);
        TPZGeoEl *gel = gmesh->Element(index);
        int64_t gelindex = gel->Index();
//        std::cout << "Processing " << gelindex << " d " << domain[gelindex] << std::endl;
        int firstside = gel->FirstSide(2);
        int lastside = gel->NSides()-1;
        if(gel->Dimension() == 2) lastside++;
        int skeletonmatid = EMHM;
        int geldomain = domain[gel->Index()];
        if(geldomain == -1) DebugStop();
        if(gel->MaterialId() == EMHM) DebugStop();
        for(int side = firstside; side < lastside; side++) {
            TPZGeoElSide gelside(gel,side);
            if(gelside.HasNeighbour(skeletonmatid)) {
//                std::cout << "Element " << gelindex << " has skeleton neighbour along side " << side << "\n";
                continue;
            }
            TPZGeoElSide neighbour = gelside.Neighbour();
            int64_t indexneigh = neighbour.Element()->Index();
            if(domain[indexneigh] == -1) {
                process.insert(indexneigh);
                SetSubdomain(neighbour.Element(), geldomain, domain);
            }
        }
    }
}
void IdentifyMHMDomain(TPZGeoMesh *gmesh, TPZVec<int> &domain)
{
    int64_t nel = gmesh->NElements();
    domain.resize(nel);
    domain.Fill(-1);
    int currentdomain = 0;
    bool found = true;
    while(found) {
        found = false;
        for (int64_t el = 0; el < nel; el++) {
            if(domain[el] != -1) continue;
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->MaterialId() == EMHM) continue;
            if(gel->Father()) continue;
            SetSubdomain(gel, currentdomain, domain);
            found = true;
            process.insert(el);
            std::cout << "Inserting " << el << " d " << currentdomain << std::endl;
            ProcessElements(gmesh, domain);
            break;
        }
        currentdomain++;
    }
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() == EMHM) continue;
        if(domain[el] == -1) {
            DebugStop();
        }
    }

}

void MakeConnectsInternal(TPZSubCompMesh *submesh);

void Substructure(TPZCompMesh *cmesh, TPZVec<int> &domains)
{
//    std::cout << domains << std::endl;
    std::map<int,TPZSubCompMesh *> submeshes;
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) {
            DebugStop();
        }
        int64_t gelindex = gel->Index();
        int domain = domains[gelindex];
        if(domain == -1) continue;
        if(submeshes.find(domain) == submeshes.end()) {
            TPZSubCompMesh *sub = new TPZSubCompMesh(*cmesh);
            submeshes[domain] = sub;
        }
        TPZSubCompMesh *sub = submeshes[domain];
        sub->TransferElement(cmesh, el);
        
    }
    for (auto it : submeshes) {
        it.second->LoadReferences();
        it.second->ExpandSolution();
        GroupElements(it.second);
    }
    cmesh->ComputeNodElCon();
    for (auto it : submeshes) {
        MakeConnectsInternal(it.second);
    }
    for (auto it : submeshes) {
        CondenseElements(it.second);
    }
    for (auto it : submeshes) {
        TPZCompMesh *cmesh = it.second;
        int64_t neq = cmesh->NEquations();
        std::cout << "Configuring submesh neq = " << neq << std::endl;
#if defined(__x86_64__) || defined(__x86_64)
        it.second->SetAnalysisSparse(global_nthread);
#elif defined(__arm__) || defined(__aarch64__)
//        it.second->SetAnalysisFStruct(0);
        it.second->SetAnalysisSkyline();
#endif
    }
    cmesh->ComputeNodElCon();
    int64_t ncon = cmesh->NConnects();
    for(int64_t ic = 0; ic<ncon; ic++)
    {
        TPZConnect &c = cmesh->ConnectVec()[ic];
        if(c.NElConnected() == 0 && c.HasDependency()) c.RemoveDepend();
    }
    cmesh->CleanUpUnconnectedNodes();
}

void GroupElements(TPZCompMesh *cmesh) {
    // look for the boundary elements and group them with their neighbour
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) {
            DebugStop();
        }
        int64_t gelindex = gel->Index();
        if(gel->Dimension() == 2) {
            if(gel->MaterialId() == EMHM) continue;
            TPZGeoElSide gelside(gel);
            TPZGeoElSide neighbour = gelside.Neighbour();
            TPZGeoEl *neigh = neighbour.Element();
            TPZCompEl *cneigh = neigh->Reference();
            if(neigh->Dimension() != 3) DebugStop();
            int firstside = neigh->FirstSide(2);
            int lastside = neigh->NSides()-1;
            std::set<TPZCompEl *> grouped;
            grouped.insert(cneigh);
            for(int side=firstside; side<lastside; side++) {
                TPZGeoElSide gelside(neigh,side);
                TPZGeoElSide neighbour = gelside.Neighbour();
                if(neighbour.Element()->Dimension() == 2 && neighbour.Element()->MaterialId() != EMHM) {
                    TPZCompEl *cel = neighbour.Element()->Reference();
                    if(!cel) DebugStop();
                    grouped.insert(cel);
                }
            }
            if(grouped.size() == 1) DebugStop();
            TPZElementGroup *celgr = new TPZElementGroup(*cmesh);
            for(auto it : grouped) {
                celgr->AddElement(it);
            }
        }
    }
}

void CondenseElements(TPZCompMesh *submesh)
{
    submesh->ComputeNodElCon();
//    {
//        std::ofstream out("submesh.txt");
//        submesh->Print(out);
//    }
    int64_t nel = submesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = submesh->Element(el);
        if(!cel) continue;
        int ncon = cel->NConnects();
        bool found = false;
        for(int ic=0; ic<ncon; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if(c.LagrangeMultiplier() == 4) {
                c.IncrementElConnected();
                found = true;
                break;
            }
        }
        if(!found) continue;
//        for(int ic=0; ic<ncon; ic++) {
//            cel->Connect(ic).Print(*submesh);
//        }
        
        TPZCondensedCompEl *condense = new TPZCondensedCompElT<STATE>(cel);
    }
    submesh->CleanUpUnconnectedNodes();
}

void MakeConnectsInternal(TPZSubCompMesh *submesh) {
    int ncon = submesh->NConnects();
    int count = 0;
    bool found = false;
    for(int ic = 0; ic<ncon; ic++) {
        TPZConnect &c = submesh->Connect(ic);
        if(c.LagrangeMultiplier() == 4) {
            c.IncrementElConnected();
            count++;
        }
        if(count == 1) {
            submesh->MakeAllInternal();
            submesh->ComputeNodElCon();
            submesh->CleanUpUnconnectedNodes();
            found = true;
            break;
        }
    }
    if(!found) {
        std::cout << "No Lagrange multiplier of level 4\n";
        DebugStop();
    }
}

// find the element groups which link to identical subdomain
void IdentifyInteractionPlanes(TPZMultiphysicsCompMesh *cmesh, TPZVec<int> &domain, const int pord, REAL& smallestVoronoiInterfaceArea) {
    std::map<std::pair<int,int>,TPZFaceDefinition> faces;
    cmesh->LoadReferences();
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(sub) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        int matid = gel->MaterialId();
        if(matid != EMHM) continue;
        TPZGeoElSide gelside(gel);
        std::set<int> leftrightdomain;
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            int64_t elindex = neighbour.Element()->Index();
            int eldomain = domain[elindex];
            leftrightdomain.insert(eldomain);
            neighbour = neighbour.Neighbour();
        }
        if(leftrightdomain.size() != 2) {
            DebugStop();
        }
        std::pair<int,int> lr = {*leftrightdomain.begin(),*leftrightdomain.rbegin()};
        faces[lr].fLeftRightDomain = lr;
        faces[lr].fElementList.push_back(cel);
    }
    
    smallestVoronoiInterfaceArea = std::numeric_limits<REAL>::max();
    for(auto &it : faces) {
        TPZFaceDefinition &face = it.second;
        face.InitializeDataStructure(pord);
        if(face.fTotalArea < smallestVoronoiInterfaceArea) smallestVoronoiInterfaceArea = face.fTotalArea;
    }
    cmesh->ExpandSolution();
    for(auto &it : faces) {
        TPZFaceDefinition &face = it.second;
        for(auto cel : face.fElementList) {
            face.Project(cel);
        }
    }
    cmesh->CleanUpUnconnectedNodes();
}

