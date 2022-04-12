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

// @proposition - Pedro
// enum ELocalMeshType {
//     ETriang = 0, ESquare = 1, ETrapezoid = 3
// };


/**
 * @brief Funcao para criar a malha computacional da tensao
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_S(TPZGeoMesh *gmesh, int pOrder);

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
 * @brief Funcao para criar a malha computacional da tensao analitica - função de Airy
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_Airy(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional multi-fisica
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder);

void CreateCondensedElements(TPZCompMesh *cmesh);

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

//Analytical solution
constexpr int solOrder{2};
auto exactSol = [](const TPZVec<REAL> &loc, TPZVec<STATE>&disp, TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];
    REAL q = 1.;
    REAL E = 1.;
    REAL I = 1.;
    REAL l = 1.;//Comprimento
    REAL c = 1.;//metade da altura
    REAL Poisson = 0.;

    REAL deflection =  5.*q*l*l*l*l/(24.*E*I)*(1+12.*c*c/(5.*l*l)*(4./5.+Poisson/2));
    disp[0] = q/(2.*E*I)*((l*l*x-x*x*x/3)*y + x*(2.*y*y*y/3.-2.*c*c*y) + Poisson*x*(y*y*y/3-c*c*y+2.*c*c*c/3));
    disp[1] = -q/(2.*E*I)*(y*y*y*y/12. - c*c*y*y/2. + 2.*c*c*c*y/3. + Poisson*((l*l-x*x)*y*y/2.+y*y*y*y/6.-c*c*y*y/5.)) 
              -q/(2.*E*I)*(l*l*x*x/2. - x*x*x*x/12. - c*c*x*x/5. + (1.+0.5*Poisson)*c*c*x*x) + deflection;


    // gradDisp(0,0) = -2*x;
    // gradDisp(1,0) = 2.*y;
    
};



TPZCompMesh *CMesh_S(TPZGeoMesh *gmesh, int pOrder) {
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    cmesh->SetDimModel(dim); //Dimesion of the model

    //Definition of the approximation space:
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
    BCond0->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha

    auto * BCond1 = material->CreateBC(material, -2, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    BCond1->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    auto * BCond2 = material->CreateBC(material, -3, dirichlet, val1, val2s); //Cria material que implementa a condicao de contorno esquerda
    BCond2->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha

    auto * BCond3 = material->CreateBC(material, -4, dirichlet, val1, val2s); //Cria material que implementa a condicao de contorno direita
    BCond3->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha

    auto * BCond5 = material->CreateBC(material, -5, neumann, val1, val2s); //Cria material que implementa a condicao de contorno direita
    BCond5->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
    cmesh->InsertMaterialObject(BCond5); //Insere material na malha

    // auto * BCond6 = material->CreateBC(material, -6, neumann, val1, val2s); //Cria material que implementa a condicao de contorno direita
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
        TPZManVector<REAL,3> x(3,0.);
        analytic2D->Elastic(x, E, nu);
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
    TPZMixedElasticityND * material = new TPZMixedElasticityND(matID, E, nu, fx, fy, plainStress, dim);

    if (TElasticityExample1::fStressState == TElasticityExample1::EAxiSymmetric) {
        material->SetAxisSymmetric();
    }

    cmesh->InsertMaterialObject(material);


    //Condições de contorno:

    TPZFMatrix<REAL> val1(dim, dim, 0.);
    TPZManVector<REAL> val2(dim, 0.);
    val2[1] = 0.;
    auto * BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    BCond0->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
    // BCond0->SetForcingFunctionBC(exactSol,4);
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    val2[1] = 0.;
    auto * BCond1 = material->CreateBC(material, -3, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    BCond1->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
    // BCond1->SetForcingFunctionBC(exactSol,4);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    //Right
    // val2[0] = 0.; // vy -> 0
    // val2[1] = 0.; // vy -> 0
    // val1(0, 0) = 0;
    // val1(1, 1) = material->BigNumber();
    auto * BCond2 = material->CreateBC(material, -4, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    BCond2->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
    // BCond2->SetForcingFunctionBC(exactSol,4);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    //Left
    // val2[0] = 1.; // vy -> 0
    // val2[1] = 0.; // vy -> 0
    // val1(0, 0) = 0;
    // val1(1, 1) = material->BigNumber();
    auto * BCond3 = material->CreateBC(material, -2, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond3->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
    // BCond3->SetForcingFunctionBC(exactSol,4);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    //Hole
    val2[1] = 0.;
    auto * BCond5 = material->CreateBC(material, -5, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond5->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
    // BCond5->SetForcingFunctionBC(exactSol,4);
    cmesh->InsertMaterialObject(BCond5); //Insere material na malha

    // //Point
    // TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
    // auto * BCPoint = material->CreateBC(material, -6, pointtype, val3, val2); //Cria material que implementa um ponto para a pressão
    // cmesh->InsertMaterialObject(BCPoint); //Insere material na malha




    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;

}

void CreateCondensedElements(TPZCompMesh *cmesh) {
    int64_t nel = cmesh->NElements();
    cmesh->ComputeNodElCon();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *> (cel);
        if (!mphys) continue;
        TPZGeoEl *gel = mphys->Reference();
        if (gel->Dimension() != cmesh->Dimension()) continue;

        int nshape = 3;
        int nstate = 1;
        int order = 1;
        int64_t newconnectindex = cmesh->AllocateNewConnect(nshape, nstate, order);
        cmesh->ConnectVec()[newconnectindex].SetLagrangeMultiplier(2);
        cmesh->ConnectVec()[newconnectindex].IncrementElConnected();
        cmesh->ConnectVec()[newconnectindex].IncrementElConnected();
        REAL delx = abs(gel->Node(1).Coord(0) - gel->Node(0).Coord(0));
        REAL dely = abs(gel->Node(1).Coord(1) - gel->Node(0).Coord(1));
        int constridf = 0;
        if (delx > dely) {
            constridf = 1;
        }
        int numel = mphys->NumberOfCompElementsInsideThisCompEl();
        if (numel != 3) DebugStop();
        TPZManVector<int> numconnects(numel, 0);
        for (int i = 0; i < numel; i++) numconnects[i] = mphys->Element(i)->NConnects();
        // removed the equation delay associated with the anti symmetric tensor
        TPZManVector<TPZCompElLagrange::TLagrange, 4> EquationDelay(3);
        // do not condense the first displacement in x, first displacement in y and displacement in x of the second connect
        EquationDelay[0].fConnect[0] = mphys->ConnectIndex(numconnects[0]);
        EquationDelay[0].fIdf[0] = 0;
        EquationDelay[0].fConnect[1] = newconnectindex;
        EquationDelay[0].fIdf[1] = 0;
        EquationDelay[1].fConnect[0] = mphys->ConnectIndex(numconnects[0]);
        EquationDelay[1].fIdf[0] = 1;
        EquationDelay[1].fConnect[1] = newconnectindex;
        EquationDelay[1].fIdf[1] = 1;
        EquationDelay[2].fConnect[0] = mphys->ConnectIndex(numconnects[0] + 1);
        EquationDelay[2].fIdf[0] = constridf;
        EquationDelay[2].fConnect[1] = newconnectindex;
        EquationDelay[2].fIdf[1] = 2;
        //        EquationDelay[3].fConnect[0] = mphys->ConnectIndex(numconnects[0] + numconnects[1] + numconnects[2] - 1);
        //        EquationDelay[3].fIdf[0] = 0;
        //        EquationDelay[3].fConnect[1] = newconnectindex;
        //        EquationDelay[3].fIdf[1] = 3;
        auto *ellag = new TPZCompElLagrange(*cmesh, EquationDelay);
        TPZElementGroup *elgr = new TPZElementGroup(*cmesh);
        elgr->AddElement(cel);
        elgr->AddElement(ellag);
        TPZCondensedCompEl *condensed = new TPZCondensedCompEl(elgr, false);
    }
    cmesh->InitializeBlock();
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
    }
    cmesh->ExpandSolution();
}



int main(int argc, char *argv[]) {
//    TPZMaterial::gBigNumber = 1.e16;

#ifdef LOG4CXX
    InitializePZLOG();
#endif


    gRefDBase.InitializeRefPatterns(2);

    EConfig conf = EThiago;
    int initial_p = 1;
    int final_p = 1;
    int initial_h = 0;
    int final_h = 0;
    bool plotting = true;
    EElementType elementType = ESquare;
    int numthreads = 0;

    switch (argc) {
        case 9:
            numthreads = atoi(argv[8]);
        case 8:
            elementType = EElementType(atoi(argv[7]));
        case 7:
            plotting = atoi(argv[6]);
        case 6:
            final_h = atoi(argv[5]);
        case 5:
            initial_h = atoi(argv[4]);
        case 4:
            final_p = atoi(argv[3]);
        case 3:
            initial_p = atoi(argv[2]);
        case 2:
            conf = EConfig(atoi(argv[1]));
    };
    int n_ref_p = final_p - initial_p + 1;
    int n_ref_h = final_h - initial_h + 1;

#ifdef USING_MKL
    mkl_set_dynamic(0); // disable automatic adjustment of the number of threads
    mkl_set_num_threads(numthreads);
#endif

    std::string rootname;
    double hx = 2, hy = 2; //Dimensões em x e y do domínio
    double x0 = -1;
    double y0 = -1;

    //Problem data:
    switch (conf) {
        case EThiago:
        case EThiagoPlus:
        case EThiagoPlusPlus:
            if(dim == 2)
            {
                TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
                elas->gE = 206.8150271873455;
                elas->gPoisson = 0.3040039545229857;
                elas->fProblemType = TElasticity2DAnalytic::EHoledPlate;
                elas->fPlaneStress = 1;
                gAnalytic = elas;
            }
            else
                DebugStop();
            hx = 1;
            hy = 1;
            x0 = 0;
            y0 = 0;

            rootname = ConfigRootname[conf] + "_Stretchx";
            break;
        case EAxiSymmetric:
        case EAxiSymmetricPlus:
        {
            if(dim != 2) DebugStop();
            TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
            elas->gE = 100;
            elas->gPoisson = 0.;
            elas->fProblemType = TElasticity2DAnalytic::ERot;
            elas->fPlaneStress = 0;
            // should be axisymetric
            DebugStop();
            gAnalytic = elas;

            hx = 2;
            hy = 2;
            x0 = 1;
            y0 = -1;
            rootname = ConfigRootname[conf] + "_Test1";
        }
            break;
        default:
            DebugStop();
            break;
    }


    for (unsigned int pref = initial_p - 1; pref < final_p; ++pref) {
        for (unsigned int href = initial_h; href <= final_h; ++href) {
            unsigned int h_level = 1 << href;
            unsigned int nelx = h_level, nely = h_level; //Number of elements in x and y directions
            std::cout << "********* " << "Number of h refinements: " << href << " (" << nelx << "x" << nely << " elements). p order: " << pref + 1 << ". *********" << std::endl;
            unsigned int nx = nelx + 1, ny = nely + 1; //Number of nodes in x and y directions
            unsigned int stressPOrder = pref + 1; //Polynomial order of the approximation
            int stressInternalPOrder = stressPOrder; //k
            if (conf == EThiagoPlus || conf == EAxiSymmetricPlus) {
                stressInternalPOrder += 1; //k+1
            }
            if (conf == EThiagoPlusPlus) {
                stressInternalPOrder += 2; //k+2
            }
            int displacementPOrder = elementType == ETriangular ? stressInternalPOrder - 1 : stressInternalPOrder;
            int rotationPOrder = displacementPOrder;
            TPZGeoMesh *gmesh = ReadMeshFromGmsh("../../mesh/holed_plate.msh");

            std::set<int> matIdCorner= {-5};
            TPZRefPatternTools::RefineDirectional(gmesh,matIdCorner);
            // TPZRefPatternTools::RefineDirectional(gmesh,matIdCorner);
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
            TPZCompMesh *cmesh_U_HDiv = CMesh_U(gmesh, displacementPOrder); //Creates the computational mesh for the displacement field
            TPZCompMesh *cmesh_P_HDiv = CMesh_P(gmesh, rotationPOrder, hx / nelx); //Creates the computational mesh for the rotation field


            TPZCompMesh *cmesh_m_HDiv = CMesh_m(gmesh, stressInternalPOrder);
#ifdef PZDEBUG
            {
                std::ofstream filecS("MalhaC_S.txt"); //Prints the stress computational mesh in txt format
                std::ofstream filecU("MalhaC_U.txt"); //Prints the displacement computational mesh in txt format
                std::ofstream filecP("MalhaC_P.txt"); //Prints the rotation computational mesh in txt format
                cmesh_S_HDiv->Print(filecS);
                cmesh_U_HDiv->Print(filecU);
                cmesh_P_HDiv->Print(filecP);
            }
#endif

            TPZManVector<TPZCompMesh*, 3> meshvector_HDiv(3);
            meshvector_HDiv[0] = cmesh_S_HDiv;
            meshvector_HDiv[1] = cmesh_U_HDiv;
            meshvector_HDiv[2] = cmesh_P_HDiv;
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
            if(1)
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
                std::ofstream filecU("MalhaC_U.txt"); //Prints the displacement computational mesh in txt format
                std::ofstream filecP("MalhaC_P.txt"); //Prints the rotation computational mesh in txt format
                meshvector[0]->Print(filecS);
                meshvector[1]->Print(filecU);
                meshvector[2]->Print(filecP);

#endif
            }
            TPZLinearAnalysis an(cmesh, optimizeBandwidth); //Creates the object that will manage the analysis of the problem
#ifdef USING_MKL
            TPZSymetricSpStructMatrix matskl(cmesh);
#else
            TPZSkylineStructMatrix<STATE> matskl(cmesh); // asymmetric case ***
#endif
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
            TElasticityExample1 example;
            an.SetExact(example.Exact());
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

            if (plotting) {
                std::string plotfile;
                {
                    std::stringstream sout;
                    sout << rootname << ".vtk";
                    plotfile = sout.str();
                }
                TPZStack<std::string> scalnames, vecnames;
                scalnames.Push("SigmaX");
                scalnames.Push("SigmaY");
                scalnames.Push("TauXY");
                // vecnames.Push("Flux");
                vecnames.Push("displacement");
                // vecnames.Push("Stress");
                int count = href * n_ref_p + pref - (initial_p - 1);
                an.SetStep(count);
                an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
                an.PostProcess(2);
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
            ErroOut << "(* Number of elements " << h_level << " *)" << std::endl;
            ErroOut << "(* Type of Element ";
            switch (elementType) {
                case ETriangular:
                    ErroOut << "triangular ";
                    break;
                case ESquare:
                    ErroOut << "square ";
                    break;
                case ETrapezoidal:
                    ErroOut << "trapezoidal ";
                    break;
            }
            ErroOut << " *)\n";
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
            output[0] = h_level;
            output[1] = stressPOrder;
            output[2] = stressInternalPOrder;
            output[3] = cmesh->NEquations();
            output[4] = cmesh->Solution().Rows();
            for (int i = 0; i < Errors.size(); i++) {
                output[5 + i] = Errors[i];
            }
            ErroOut << "Error[[" << href + 1 << "," << pref + 1 << "]] = {" << output << "};\n";

            std::cout << "Errors = " << Errors << std::endl;


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

        }
    }

  
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
        // stringtoint[1]["Hole"] = -5;
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