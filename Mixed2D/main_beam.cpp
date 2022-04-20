#include <fstream>

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include <TPZGeoElement.h>

#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZCompelDiscScaled.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"

#include <pz_config.h>

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "pzgeoel.h"
#include "pzgeopoint.h"
#include "tpzgeoblend.h"


#include "TPZInterfaceEl.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZLinearAnalysis.h"
#include <tpzarc3d.h>

#include "pzstepsolver.h"
#include "TPZStructMatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"

#include "TPZNullMaterial.h"
#include "meshgen.h"
#include "TPZGmshReader.h"

#include "TPZMixedElasticityND.h"
#include "TPZMixedElasticityNDAiry.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZCompElLagrange.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZSpStructMatrix.h"
#include "pzlog.h"
#include <iostream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzfunction.h"
#include "TPZReadGIDGrid.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

#include "TPZGenGrid2D.h"
#include "TPZExtendGridDimension.h"

#include "TPZAnalyticSolution.h"
#include "TPZRefPatternTools.h"
#include "TPZRefPatternDataBase.h"
#include "TPZCompElH1.h"
#include "AiryFunctionHoledPlate.h"
#include "AiryFunctionLoadedBeam.h"
#include "TPZCompElDiscStress.h"
#include "TPZCompElDiscStressBound.h"
#include "TPZCompMeshTools.h"

#include "TPZMixedElasticityCMeshCreator.h"
#include "TPZMixedElasticityUtils.h"

#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
#endif

TPZAnalyticSolution *gAnalytic = 0;
//------------------Elasticity Problem------------------------
// local enum for mesh types @ToDo these names might lead to confusion. We should consider changing.
enum EElementType {
    ETriangular = 0, ESquare = 1, ETrapezoidal = 2
};

/**
 * @brief Funcao para criar a malha computacional da tensao
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_Airy(TPZGeoMesh *gmesh, int pOrder);


void Error(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);


TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

//Variáveis globais do problema:

const int dim = 2; // Dimension of the problem
const int matID = 1; // Material of the volumetric element
const int matLagrange = -10; // Material of the Lagrange multipliers
const int dirichlet = 0, neumann = 1, mixed = 2, dirichletvar = 4, pointtype = 5; // Boundary conditions of the problem ->default: Dirichlet on left and right

using namespace std;


int main(int argc, char *argv[]) {
//    TPZMaterial::gBigNumber = 1.e16;

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    gRefDBase.InitializeRefPatterns(2);

    constexpr int pOrder{1};

    EElementType elementType = ESquare;
    int numthreads = 0;

    // std::string rootname;

    //Problem data:
    TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
    // elas->gE = 206.8150271873455;
    elas->gE = 1.;
    // elas->gPoisson = 0.3040039545229857; 
    elas->gPoisson = 0.0;
    elas->fProblemType = TElasticity2DAnalytic::EUniformLoadedBeam;
    elas->fPlaneStress = 0;
    gAnalytic = elas;

    unsigned int stressPOrder = pOrder; //Polynomial order of the approximation
    int stressInternalPOrder = stressPOrder; //k
    stressInternalPOrder += 2; //k+2

    int displacementPOrder = elementType == ETriangular ? stressInternalPOrder - 1 : stressInternalPOrder;
    int rotationPOrder = displacementPOrder;
    TPZGeoMesh *gmesh = ReadMeshFromGmsh("../../mesh/plate.msh");

    // std::set<int> matIdCorner= {-5};
    // TPZRefPatternTools::RefineDirectional(gmesh,matIdCorner);

    TPZMixedElasticityUtils utils;

    TPZMixedElasticityCMeshCreator mElastCreator;
    mElastCreator.SetDomainMaterialId(matID);
    std::set<int> bcDirichlet = {-2,-3,-4};
    std::set<int> bcNeumann = {-1};
    mElastCreator.SetDirichletMaterialId(bcDirichlet); 
    mElastCreator.SetNeumannMaterialId(bcNeumann);

    mElastCreator.CreateGeometricWrapBCElements(gmesh);


#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Prints the geometric mesh in txt format
    std::ofstream filegvtk("MalhaGeo.vtk"); //Prints the geometric mesh in vtk format
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk, true);
#endif
    //Creating computational mesh:
    TPZCompMesh *cmesh_S_HDiv = mElastCreator.CMesh_S(gmesh,stressPOrder,HDivFamily::EHDivStandard,gAnalytic); //Creates the computational mesh for the stress field

    mElastCreator.ChangeInternalOrder(cmesh_S_HDiv, stressInternalPOrder);
    auto charac_size = cmesh_S_HDiv->MaximumRadiusOfMesh() / sqrt(2.);
    if (cmesh_S_HDiv->ApproxSpace().HDivFam() == HDivFamily::EHDivConstant){
        displacementPOrder = 0;
    }
    TPZCompMesh *cmesh_Airy = CMesh_Airy(gmesh, stressPOrder); //Creates the computational mesh for the stress field
    TPZCompMesh *cmesh_U_HDiv = mElastCreator.CMesh_U(gmesh, displacementPOrder); //Creates the computational mesh for the displacement field
    TPZCompMesh *cmesh_P_HDiv = mElastCreator.CMesh_P(gmesh, rotationPOrder, charac_size); //Creates the computational mesh for the rotation field
    
    std::cout << "Mesh Characteristic size = " << charac_size << std::endl;
    //Multiphysics Cmesh
    TPZCompMesh *cmesh_m_HDiv = mElastCreator.CMesh_m(gmesh, stressInternalPOrder,gAnalytic);

#ifdef PZDEBUG
    {
        std::ofstream filecS("MalhaC_S.txt"); //Prints the stress computational mesh in txt format
        std::ofstream filecAiry("MalhaC_Airy.txt"); //Prints the displacement computational mesh in txt format
        std::ofstream filecU("MalhaC_U.txt"); //Prints the displacement computational mesh in txt format
        std::ofstream filecP("MalhaC_P.txt"); //Prints the rotation computational mesh in txt format
        cmesh_S_HDiv->Print(filecS);
        cmesh_Airy->Print(filecAiry);
        cmesh_U_HDiv->Print(filecU);
        cmesh_P_HDiv->Print(filecP);
    }
#endif

    TPZManVector<TPZCompMesh*, 4> meshvector_HDiv(4);
    meshvector_HDiv[0] = cmesh_S_HDiv;
    meshvector_HDiv[1] = cmesh_U_HDiv;
    meshvector_HDiv[2] = cmesh_P_HDiv;
    meshvector_HDiv[3] = cmesh_Airy;

    TPZBuildMultiphysicsMesh::AddElements(meshvector_HDiv, cmesh_m_HDiv);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector_HDiv, cmesh_m_HDiv);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector_HDiv, cmesh_m_HDiv);
    cmesh_m_HDiv->LoadReferences();
    
#ifdef PZDEBUG
    std::ofstream fileg1("MalhaGeo2.txt");
    gmesh->Print(fileg1); //Prints the geometric mesh in txt format

    std::ofstream filecm("MalhaC_m.txt");
    cmesh_m_HDiv->Print(filecm); //Prints the multi-physics computational mesh in txt format
#endif

    //Solving the system:
    bool optimizeBandwidth = true;
    cmesh_m_HDiv->InitializeBlock();
    
    TPZCompMesh *cmesh = cmesh_m_HDiv;
    TPZManVector<TPZCompMesh *> meshvector = meshvector_HDiv;
    if(cmesh_S_HDiv->ApproxSpace().HDivFam() != HDivFamily::EHDivKernel)
    {
        TPZCompMesh * cmesh_m_Hybrid;
        TPZManVector<TPZCompMesh*, 4> meshvector_Hybrid(4);
        TPZHybridizeHDiv hybridizer;
        bool group_element = false;
        tie(cmesh_m_Hybrid, meshvector_Hybrid) = hybridizer.Hybridize(cmesh_m_HDiv, meshvector_HDiv, group_element, -1.);
        cmesh_m_Hybrid->InitializeBlock();
        cmesh = cmesh_m_Hybrid;
        meshvector = meshvector_Hybrid;
#ifdef PZDEBUG
        std::ofstream out("MalhaC_hybrid.txt");
        cmesh->Print(out);
        std::ofstream filecS("MalhaC_S.txt"); //Prints the stress computational mesh in txt format
        std::ofstream filecAiry("MalhaC_Airy.txt"); //Prints the displacement computational mesh in txt format
        std::ofstream filecU("MalhaC_U.txt"); //Prints the displacement computational mesh in txt format
        std::ofstream filecP("MalhaC_P.txt"); //Prints the rotation computational mesh in txt format
        meshvector[0]->Print(filecS);
        meshvector[1]->Print(filecU);
        meshvector[2]->Print(filecP);
        meshvector[3]->Print(filecAiry);
#endif
    }


    // TPZCompMeshTools::CreatedCondensedElements(cmesh,true,false);
    TPZLinearAnalysis an(cmesh, optimizeBandwidth); //Creates the object that will manage the analysis of the problem
    
    utils.SolveProblem(an,matID);
    
    
    TPZMatErrorCombinedSpaces<STATE> *materr = dynamic_cast<TPZMatErrorCombinedSpaces<STATE> *>(cmesh_m_HDiv->FindMaterial(matID));
    TPZManVector<REAL, 6> Errors(materr->NEvalErrors());

    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);


    utils.PrintResults(an);
    

#ifdef PZDEBUG
    //Imprimindo vetor solução:
    {
        TPZFMatrix<STATE> solucao = cmesh->Solution(); //Pegando o vetor de solução, alphaj
        std::ofstream solout("sol.nb");
        solucao.Print("Sol", solout, EMathematicaInput); //Imprime na formatação do Mathematica

        std::ofstream fileAlpha("alpha.nb");
        an.Solution().Print("Alpha = ", fileAlpha, EMathematicaInput);
    }
#endif


    //    //Calculo do erro
    //    std::cout << "Computing Error " << std::endl;

    std::stringstream sout;
    // sout << rootname;
    switch (elementType) {
        case ETriangular:
            sout << "_tria";
            break;
        case ESquare:
            sout << "_quad";
            break;
        case ETrapezoidal:
            sout << "_trap";
            break;
    }
    sout << "_" << stressPOrder << "_Error.nb";
    ofstream ErroOut(sout.str(), std::ios::app);
    // ErroOut << "(* Type of simulation " << rootname << " *)\n";
    ErroOut << "(* Number of Condensed equations " << cmesh->NEquations() << " *)" << std::endl;
    ErroOut << "(* Number of equations before condensation " << cmesh->Solution().Rows() << " *)" << std::endl;
    ErroOut << "(*\n";
    an.SetExact(gAnalytic->ExactSolution());
    an.SetThreadsForError(numthreads);
    bool store_errors = true;
    cmesh->ElementSolution().Redim(cmesh->NElements(), Errors.size());
    std::cout << "Computing errors." << std::endl;
    an.PostProcessError(Errors, store_errors, ErroOut);
    std::cout << "Computed errors." << std::endl;
    ErroOut << "nelx ribporder internalporder n_condensed - n_total - error_sigma - error_energy - error_div_sigma - error_u - error_r - error_as\n";
    ErroOut << "*)\n";
    TPZManVector<STATE, 10> output(Errors.size() + 5, 0);
    output[0] = 0;
    output[1] = stressPOrder;
    output[2] = stressInternalPOrder;
    output[3] = cmesh->NEquations();
    output[4] = cmesh->Solution().Rows();
    for (int i = 0; i < Errors.size(); i++) {
        output[5 + i] = Errors[i];
    }
    ErroOut << "Error[[" << charac_size << "," << pOrder << "]] = {" << std::scientific << std::setprecision(15) << output << "};\n";

    std::cout << "Errors = " << std::scientific << std::setprecision(15) << Errors << std::endl;


    an.CleanUp();

    if(cmesh != cmesh_m_HDiv)
    {
        delete cmesh;
        for (int i = meshvector.size() - 1; i >= 0; i--) {
            meshvector[i]->CleanUp();
            delete meshvector[i];
        }
    }
    delete cmesh_m_HDiv;
    for (int i = meshvector_HDiv.size() - 1; i >= 0; i--) {
        meshvector_HDiv[i]->CleanUp();
        delete meshvector_HDiv[i];
    }
    delete gmesh;


    std::cout << "FINISHED!" << std::endl;

    return 0;
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
        
        stringtoint[2]["Domain"] = 1;
        stringtoint[1]["Lower"] = -1;
        stringtoint[1]["Right"] = -2;
        stringtoint[1]["Upper"] = -3;
        stringtoint[1]["Left"] = -4;
        stringtoint[1]["Hole"] = -5;
        stringtoint[0]["Corner"] = -6;

        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh);

        //Prints gmesh mesh properties
        std::string vtk_name = "geoMesh.vtk";
        std::ofstream vtkfile(vtk_name.c_str());

        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
    }

    return gmesh;
}

TPZCompMesh *CMesh_Airy(TPZGeoMesh *gmesh, int pOrder) {
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    cmesh->SetDimModel(dim); //Dimesion of the model

    std::set<int> matids = {matID,-1,-2,-3,-4,-5};

    for (std::set<int>::iterator it=matids.begin(); it!=matids.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        mat->SetDimension(dim);
        mat->SetNStateVariables(1);
        cmesh->InsertMaterialObject(mat);
    }
     
    //Creates the analytical functions
    TPZVec<REAL> axis(3,0.);
    REAL length = 2.;
    REAL height = 0.5;
    auto AiryFunction = new AiryFunctionLoadedBeam(axis,height,length);
    AiryFunction->SetElasticConstants(1., 0., 1., 1.);
    
    for (int64_t el = 0; el < gmesh->NElements(); el++) {
        TPZGeoEl *gel = gmesh->Element(el);

        if(!gel) DebugStop();
        auto matid = gel->MaterialId();

        if (matid == matID){
            // Create volumetric elements
            auto cel = new TPZCompElDiscStress<AiryFunctionLoadedBeam>(*cmesh,gel);
            cel->SetStressDef(AiryFunction);
        } else {
            auto cel = new TPZCompElDiscStressBound<AiryFunctionLoadedBeam>(*cmesh,gel);
            cel->SetStressDef(AiryFunction);
        }
    }  

    cmesh->AutoBuild();


    return cmesh;


}



