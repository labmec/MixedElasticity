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
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoblend.h"

#include "tpzintpoints.h"

#include "TPZInterfaceEl.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzautopointer.h"
#include "TPZBndCondT.h"
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
#include "pzshapepoint.h"
#include "TPZCompElH1.h"
using namespace pzshape;
#include "TPZCompElKernelHDiv.h"
#include "AiryFunctionHoledPlate.h"
#include "AiryFunctionLoadedBeam.h"
#include "TPZCompElDiscStress.h"

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
TPZCompMesh *CMesh_S(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional da tensao
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_Airy(TPZGeoMesh *gmesh, int pOrder);

/// change the order of the internal connect to the given order
void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional do deslocamento
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_U(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional da rotacao
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 * @param elementdimension tamanho de elementos
 */
TPZCompMesh *CMesh_P(TPZGeoMesh *gmesh, int pOrder, REAL elementdimension);

/**
 * @brief Funcao para criar a malha computacional multi-fisica
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder);

void Error(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);


TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

//Variáveis globais do problema:

const int dim = 2; // Dimension of the problem
const int matID = 1; // Material of the volumetric element
const int matLagrange = -10; // Material of the Lagrange multipliers
const int matBCbott = -1, matBCtop = -2, matBCleft = -3, matBCright = -4; // Materials of the boundary conditions
const int dirichlet = 0, neumann = 1, mixed = 2, dirichletvar = 4, pointtype = 5; // Boundary conditions of the problem ->default: Dirichlet on left and right

using namespace std;

enum EConfig {
    EThiago = 0, EAxiSymmetric = 1, EThiagoPlus = 2, EAxiSymmetricPlus = 3, EThiagoPlusPlus = 4
};

std::string ConfigRootname[5] = {
    "Mixed",
    "Mixed_AxiSymmetric",
    "MixedPlus",
    "Mixed_AxiSymmetricPlus",
    "MixedPlusPlus"
};


TPZCompMesh *CMesh_S(TPZGeoMesh *gmesh, int pOrder) {
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    cmesh->SetDimModel(dim); //Dimesion of the model

    //Definition of the approximation space:
    // cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivKernel);
    // cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
    cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivStandard);
    cmesh->SetAllCreateFunctionsHDiv(); //Creating H(div) functions:

    //Criando material cujo nSTATE = 2:

    TPZNullMaterial<> * material = new TPZNullMaterial<>(matID);
    material->SetNStateVariables(dim);
    material->SetDimension(dim);
    cmesh->InsertMaterialObject(material); //Insere material na malha



    //Boundary conditions:
    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE> val2s(2, 0.), val2(2, 0.);
    // val2s[0] = 10.0; // vx -> 0
    // val2s[1] = 0.0; // vy -> 0

    auto * BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    BCond0->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha

    auto * BCond1 = material->CreateBC(material, -2, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    BCond1->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    auto * BCond2 = material->CreateBC(material, -3, dirichlet, val1, val2s); //Cria material que implementa a condicao de contorno esquerda
    BCond2->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha

    auto * BCond3 = material->CreateBC(material, -4, dirichlet, val1, val2s); //Cria material que implementa a condicao de contorno direita
    BCond3->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha

    auto * BCond5 = material->CreateBC(material, -5, dirichlet, val1, val2s); //Cria material que implementa a condicao de contorno direita
    BCond5->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond5); //Insere material na malha

    // auto * BCond6 = material->CreateBC(material, -6, pointtype, val1, val2s); //Cria material que implementa a condicao de contorno direita
    // cmesh->InsertMaterialObject(BCond6); //Insere material na malha

    val2s[0] = 0;
    val2s[1] = 0;
    cmesh->InsertMaterialObject(material->CreateBC(material, matLagrange, dirichlet, val1, val2s)); //Insere material na malha

    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    int ncel = cmesh->NElements();
    for (int i = 0; i < ncel; i++) {
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if (!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *> (compEl);
        if (facel)DebugStop();
    }

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;

}

TPZCompMesh *CMesh_U(TPZGeoMesh *gmesh, int pOrder) {
    if (pOrder == 0) {
        // Constant-per-element
        TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
        cmesh->SetDimModel(dim); //Insere dimensão do modelo

        cmesh->SetAllCreateFunctionsDiscontinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);

        
        TPZNullMaterial<> * material = new TPZNullMaterial<>(matID);
        material->SetDimension(dim);
        material->SetNStateVariables(dim);

        cmesh->InsertMaterialObject(material); //Insere material na malha

        std::set<int> materialids;
        materialids.insert(matID);
        {
            gmesh->ResetReference();
            int64_t nel = gmesh->NElements();
            for (int64_t el = 0; el < nel; el++) {
                TPZGeoEl *gel = gmesh->Element(el);
                if (!gel)continue;
                int matid = gel->MaterialId();
                if (materialids.find(matid) == materialids.end()) {
                    continue;
                }
                new TPZCompElDiscScaled(*cmesh, gel);
                gel->ResetReference();
            }
        }

        int ncon = cmesh->NConnects();
        for (int i = 0; i < ncon; i++) {
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }

        int64_t nelem = cmesh->NElements();
        for (int64_t el = 0; el < nelem; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZCompElDiscScaled *disc = dynamic_cast<TPZCompElDiscScaled *> (cel);
            if (!disc) {
                continue;
            }
            disc->SetTotalOrderShape();
            disc->SetTrueUseQsiEta();
        }

        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
        return cmesh;
    } else {
        //Criando malha computacional:
        TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
        cmesh->SetDimModel(dim); //Insere dimensão do modelo

        cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
        cmesh->ApproxSpace().CreateDisconnectedElements(true);


        TPZNullMaterial<> * material = new TPZNullMaterial<>(matID);
        material->SetDimension(dim);
        material->SetNStateVariables(dim);

        cmesh->InsertMaterialObject(material); //Insere material na malha

        //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha

        int ncel = cmesh->NElements();
        for (int i = 0; i < ncel; i++) {
            TPZCompEl * compEl = cmesh->ElementVec()[i];
            if (!compEl) continue;
            TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *> (compEl);
            if (facel)DebugStop();

        }
        std::set<int> materialids;
        materialids.insert(matID);
        cmesh->AutoBuild(materialids);
        cmesh->LoadReferences();
        cmesh->ApproxSpace().CreateDisconnectedElements(false);
        cmesh->AutoBuild();

        int ncon = cmesh->NConnects();
        for (int i = 0; i < ncon; i++) {
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }

        //    cmesh->AdjustBoundaryElements();
        //    cmesh->CleanUpUnconnectedNodes();

        return cmesh;
    }
}

TPZCompMesh *CMesh_P(TPZGeoMesh *gmesh, int pOrder, REAL elementdim) {
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo

    cmesh->SetAllCreateFunctionsDiscontinuous();

    //    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    //Criando material cujo nSTATE = 1:
    TPZNullMaterial<> *material = new TPZNullMaterial<>(matID); //criando material que implementa a formulacao fraca do problema modelo
    material->SetDimension(dim);
    if(dim == 3)
    {
        material->SetNStateVariables(3);
    }

    cmesh->InsertMaterialObject(material); //Insere material na malha



    std::set<int> materialids;
    materialids.insert(matID);
    //materialids.insert(3);
    {
        gmesh->ResetReference();
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel)continue;
            int matid = gel->MaterialId();
            if (materialids.find(matid) == materialids.end()) {
                continue;
            }
            new TPZCompElDiscScaled(*cmesh, gel);
            gel->ResetReference();
        }
    }

    //cmesh->LoadReferences();
    //    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    //    cmesh->AutoBuild();


    int ncon = cmesh->NConnects();
    for (int i = 0; i < ncon; i++) {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }

    int64_t nelem = cmesh->NElements();
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZCompElDiscScaled *disc = dynamic_cast<TPZCompElDiscScaled *> (cel);
        if (!disc) {
            continue;
        }
        disc->SetTotalOrderShape();
        disc->SetFalseUseQsiEta();
        disc->SetConstC(elementdim);
        disc->SetScale(1./elementdim);
    }

    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    return cmesh;

}

TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder) {

    //Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();

    // Criando material:

    // example is initialized in the calling method
    //    example.fProblemType = TElasticityExample1::EThiago;
    //    example.fStressState   = TElasticityExample1::EPlaneStrain;

    REAL E = 0.; //* @param E elasticity modulus
    REAL nu = 0.; //* @param nu poisson coefficient

    TElasticity2DAnalytic * analytic2D = 0;
    TElasticity3DAnalytic * analytic3D = 0;
    if(dim == 2)
    {
        analytic2D = dynamic_cast<TElasticity2DAnalytic *>(gAnalytic);
        if(!analytic2D) DebugStop();
        TPZVec<REAL> x;
        x.Resize(3);
        x.Fill(0.);
        analytic2D->Elastic<REAL>(x, E, nu);
    }
    else if(dim == 3)
    {
        analytic3D = dynamic_cast<TElasticity3DAnalytic *>(gAnalytic);
        if(!analytic3D) DebugStop();
        TPZManVector<REAL,3> x(3,0.);
        analytic3D->Elastic(x, E, nu);
    }

    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = 0.; //* @param fx forcing function \f$ -x = fx \f$
    int plainStress = 0; //* @param plainstress = 1 \f$ indicates use of plainstress
    if(dim == 2)
    {
        if (analytic2D->fPlaneStress == 0) {
            plainStress = 0;
        } else {
            plainStress = 1;
        }
    }
    // TPZMixedElasticityND * material = new TPZMixedElasticityND(matID, E, nu, fx, fy, plainStress, dim);
    TPZMixedElasticityNDAiry * material = new TPZMixedElasticityNDAiry(matID, E, nu, fx, fy, plainStress, dim);

    material->SetForcingFunction(gAnalytic->ForceFunc(),4);
    material->SetExactSol(gAnalytic->ExactSolution(),4);     

    cmesh->InsertMaterialObject(material);


    //Condições de contorno:

    TPZFMatrix<REAL> val1(dim, dim, 0.);
    TPZManVector<REAL> val2(dim, 0.);
    val2[1] = 0.;
    auto * BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    BCond0->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    val2[1] = 0.;
    auto * BCond1 = material->CreateBC(material, -3, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    BCond1->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    //Right
    // val2[0] = 0.; // vy -> 0
    // val2[1] = 0.; // vy -> 0
    // val1(0, 0) = 0;
    // val1(1, 1) = material->BigNumber();
    auto * BCond2 = material->CreateBC(material, -4, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    BCond2->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    //Left
    // val2[0] = 1.; // vy -> 0
    // val2[1] = 0.; // vy -> 0
    // val1(0, 0) = 0;
    // val1(1, 1) = material->BigNumber();
    auto * BCond3 = material->CreateBC(material, -2, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond3->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    //Hole
    val2[1] = 0.;
    auto * BCond5 = material->CreateBC(material, -5, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond5->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond5); //Insere material na malha

    // //Point
    // TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
    // auto * BCPoint = material->CreateBC(material, -6, pointtype, val3, val2); //Cria material que implementa um ponto para a pressão
    // cmesh->InsertMaterialObject(BCPoint); //Insere material na malha

    auto * BCond7 = material->CreateBC(material, matLagrange, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond7->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond7); //Insere material na malha


    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;

}

/// change the order of the internal connect to the given order

void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder) {
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) continue;

        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != cmesh->Dimension()) {
            continue;
        }
        int nc = cel->NConnects();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel) DebugStop();

        intel->ForceSideOrder(gel->NSides() - 1, pOrder);

        if (cmesh->ApproxSpace().HDivFam() == HDivFamily::EHDivKernel){
            TPZCompElKernelHDiv<TPZShapeQuad> *celKernel = dynamic_cast<TPZCompElKernelHDiv<TPZShapeQuad> *> (cel);
            for (int icon = 0; icon < celKernel->NConnects(); icon++)
            {
                TPZConnect &c = celKernel->Connect(icon);
                int nShapeF = celKernel->NConnectShapeF(icon,c.Order());
                c.SetNShape(nShapeF);
                int64_t seqnum = c.SequenceNumber();
                int nvar = 1;
                TPZMaterial * mat = celKernel->Material();
                if (mat) nvar = mat->NStateVariables();
                celKernel->Mesh()->Block().Set(seqnum, nvar * nShapeF);
                celKernel->AdjustIntegrationRule();
            }    
        }
        
    }
    cmesh->ExpandSolution();
}



int main(int argc, char *argv[]) {
//    TPZMaterial::gBigNumber = 1.e16;

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    gRefDBase.InitializeRefPatterns(2);

    constexpr int pOrder{1};

    EConfig conf = EThiago;
    EElementType elementType = ESquare;
    int numthreads = 0;

    std::string rootname;

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
    stressInternalPOrder += 0; //k+2

    int displacementPOrder = elementType == ETriangular ? stressInternalPOrder - 1 : stressInternalPOrder;
    int rotationPOrder = displacementPOrder;
    TPZGeoMesh *gmesh = ReadMeshFromGmsh("../../mesh/plate.msh");

    // std::set<int> matIdCorner= {-5};
    // TPZRefPatternTools::RefineDirectional(gmesh,matIdCorner);
    

#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Prints the geometric mesh in txt format
    std::ofstream filegvtk("MalhaGeo.vtk"); //Prints the geometric mesh in vtk format
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk, true);
#endif
    //Creating computational mesh:
    TPZCompMesh *cmesh_S_HDiv = CMesh_S(gmesh, stressPOrder); //Creates the computational mesh for the stress field
    ChangeInternalOrder(cmesh_S_HDiv, stressInternalPOrder);
    auto charac_size = cmesh_S_HDiv->MaximumRadiusOfMesh() / sqrt(2.);
    if (cmesh_S_HDiv->ApproxSpace().HDivFam() == HDivFamily::EHDivConstant){
        displacementPOrder = 0;
    }
    TPZCompMesh *cmesh_Airy = CMesh_Airy(gmesh, stressPOrder); //Creates the computational mesh for the stress field
    TPZCompMesh *cmesh_U_HDiv = CMesh_U(gmesh, displacementPOrder); //Creates the computational mesh for the displacement field
    TPZCompMesh *cmesh_P_HDiv = CMesh_P(gmesh, rotationPOrder, charac_size); //Creates the computational mesh for the rotation field
    
    std::cout << "Mesh Characteristic size = " << charac_size << std::endl;
    //Multiphysics Cmesh
    TPZCompMesh *cmesh_m_HDiv = CMesh_m(gmesh, stressInternalPOrder);

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
    meshvector_HDiv[1] = cmesh_Airy;
    meshvector_HDiv[2] = cmesh_U_HDiv;
    meshvector_HDiv[3] = cmesh_P_HDiv;
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
        TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
        TPZHybridizeHDiv hybridizer;
        bool group_element = true;
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
        meshvector[1]->Print(filecAiry);
        meshvector[2]->Print(filecU);
        meshvector[3]->Print(filecP);
#endif
    }
    
    TPZLinearAnalysis an(cmesh, optimizeBandwidth); //Creates the object that will manage the analysis of the problem
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>>  matskl(cmesh);

    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);

    std::cout << "Assemble matrix with NDoF = " << cmesh->NEquations() << "." << std::endl;
    an.Assemble(); //Assembles the global stiffness matrix (and load vector)
    std::cout << "Assemble finished." << std::endl;

    TPZMatErrorCombinedSpaces<STATE> *materr = dynamic_cast<TPZMatErrorCombinedSpaces<STATE> *>(cmesh_m_HDiv->FindMaterial(matID));
    TPZManVector<REAL, 6> Errors(materr->NEvalErrors());
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
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);

    if (1) {
        std::string plotfile;
        {
            std::stringstream sout;
            sout << rootname << "results.vtk";
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
        // vecnames.Push("Stress");
        int count = 0;
        an.SetStep(count);
        an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
        an.PostProcess(0);
    }

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

    //   matids.clear();
    //   matids.insert(-1);
    //   TPZManVector<STATE,3> result;
    //  result = cmesh_m_HDiv->Integrate("state",matids);
    //  std::cout << "Sigma Y"  << result << std::endl;


    //    //Calculo do erro
    //    std::cout << "Computing Error " << std::endl;

    std::stringstream sout;
    sout << rootname;
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
    ErroOut << "(* Type of simulation " << rootname << " *)\n";
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

    TPZNullMaterial<> * material = new TPZNullMaterial<>(matID);
    material->SetNStateVariables(dim);
    material->SetDimension(dim);
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
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

        if (matid != matID) continue;

        // Create volumetric elements
        auto cel = new TPZCompElDiscStress<AiryFunctionLoadedBeam>(*cmesh,gel);
        cel->SetStressDef(AiryFunction);
    }
    // cmesh->SetAllCreateFunctionsDiscontinuous();

    

    cmesh->AutoBuild();

    // //Creates the analytical functions
    // TPZVec<REAL> axis(3,0.);
    // REAL length = 2.;
    // REAL height = 0.5;    
    
    // std::function<AiryFunctionLoadedBeam()> airyFunction = [=](){
    //     auto aux = new AiryFunctionLoadedBeam(axis,height,length);
    //     aux->SetElasticConstants(1., 0., 1., 1.);
    //     return *aux;
    // };

    // //Sets the external shape functions (Airy) to the CompEls
    // for (auto cel:cmesh->ElementVec())
    // {
    //     auto disc = dynamic_cast<TPZCompElDiscStress *> (cel);
    //     disc->SetExternalShapeFunction(airyFunction);
    // }
    

    return cmesh;


}