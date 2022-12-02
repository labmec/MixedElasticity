#include <fstream>

#include "pzgeoel.h"
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

#include "TPZMixedElasticityND.h"
#include "TPZInterfaceEl.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzautopointer.h"
#include "TPZBndCondT.h"
#include "TPZLinearAnalysis.h"
#include <tpzarc3d.h>

#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"

#include "TPZNullMaterial.h"
#include "meshgen.h"
#include "TPZGmshReader.h"

#include "pzmixedelasmat.h"
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
#include "TPZGeoMeshTools.h"

#include "TPZGenGrid2D.h"
#include "TPZExtendGridDimension.h"

#include "TPZAnalyticSolution.h"

#include "TPZMHMixedHybridMeshControl.h"

#include "TPZMultiphysicsCompMesh.h"


//#define USENEWVTK
#ifdef USENEWVTK
#include "TPZVTKGenerator.h"
#endif

#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
static LoggerPtr loggerconfig(Logger::getLogger("pz.hdiv.vecconfig"));
#endif

// -------------------- Global Variables --------------------
//constexpr int Globnx{65}, Globny{65}, Globnz{33};
//constexpr REAL Globpartsize{156.25};

constexpr int Globnx{17}, Globny{17}, Globnz{9};
constexpr REAL Globpartsize{625};

//------------------Elasticity Problem------------------------

// local enum for mesh types @ToDo these names might lead to confusion. We should consider changing.
enum EElementType {
    ETriangular = 0, ESquare = 1, ETrapezoidal = 2, ETetraheda = 3
};

// @proposition - Pedro
// enum ELocalMeshType {
//     ETriang = 0, ESquare = 1, ETrapezoid = 3
// };

/**
 * @brief Funcao para criar a malha geometrica do problema a ser simulado
 * @note A malha sera unidim5ensional formada por nel elementos de tamanho elsize
 * @param uNDiv number of divisions ortogonal to the plates performed on the domain
 * @param vNDiv number of divisions parallel to the plates performed on the domain
 * @param nel numero de elementos
 * @param elsize tamanho dos elementos
 */
TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType);

/**
* @brief Funcao para criar a malha geometrica do problema a ser simulado
* @note A malha sera tridimensional formada por nel elementos de tamanho elsize
* @param nelx, nely number of elements in the x and y direction (the number of elements in the z direction = nelx
* @param hx, hy size of the domain in x and y (size in z = hx)
* @param x0, y0 bottom left point coordinate (bottom left in z = 0)
* @param meshType = triangle, quadrilateral or trapeze
*/
TPZGeoMesh *CreateGMesh3D(int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType);

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
 * @brief Funcao para criar a malha computacional multi-fisica
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMesh_AxiS(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional multi-fisica ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_Girk(TPZGeoMesh *gmesh, int pOrder);

void Error(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);

void InsertMaterialObjects(TPZCompMesh &cmeshref);

//Variáveis globais do problema:

const int dim = 3; // Dimension of the problem
const int matID = 1; // Material of the volumetric element
const int matLagrange = -10; // Material of the Lagrange multipliers
const int matBCbott = -1, matBCtop = -2, matBCleft = -3, matBCright = -4, matBCfront = -5, matBCback = -6; // Materials of the boundary conditions
const int dirichlet = 0, neumann = 1, mixed = 2, dirichletvar = 4, pointtype = 5; // Boundary conditions of the problem ->default: Dirichlet on left and right

using namespace std;

//enum EConfig {
//    EThiago = 0, EAxiSymmetric = 1, EThiagoPlus = 2, EAxiSymmetricPlus = 3, EThiagoPlusPlus = 4
//};

std::string ConfigRootname[5] = {
    "Mixed",
    "Mixed_AxiSymmetric",
    "MixedPlus",
    "Mixed_AxiSymmetricPlus",
    "MixedPlusPlus"
};

TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType) {
    //Creating geometric mesh, nodes and elements.
    //Including nodes and elements in the mesh object:
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);

    //Auxiliary vector to store coordinates:
    //TPZVec <REAL> coord(3, 0.);
    TPZVec<REAL> gcoord1(3, 0.);
    TPZVec<REAL> gcoord2(3, 0.);
    gcoord1[0] = x0;
    gcoord1[1] = y0;
    gcoord1[2] = 0;
    gcoord2[0] = x0 + hx;
    gcoord2[1] = y0 + hy;
    gcoord2[2] = 0;
    //Inicialização dos nós:

    TPZManVector<int> nelem(2, 1);
    nelem[0] = nelx;
    nelem[1] = nely;

    TPZGenGrid2D gengrid(nelem, gcoord1, gcoord2);

    switch (meshType) {
        case ETriangular:
            gengrid.SetElementType(MMeshType::ETriangular);
            break;
        case ESquare:
            gengrid.SetElementType(MMeshType::EQuadrilateral);
            break;
        case ETrapezoidal:
            gengrid.SetDistortion(0.25);
            break;
    }

    gengrid.Read(gmesh, matID);
    gengrid.SetBC(gmesh, 4, matBCbott);
    gengrid.SetBC(gmesh, 5, matBCright);
    gengrid.SetBC(gmesh, 6, matBCtop);
    gengrid.SetBC(gmesh, 7, matBCleft);

    gmesh->BuildConnectivity();
    {
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
    }

    //Printing geometric mesh:

    //ofstream bf("before.vtk");
    //TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;



}

TPZGeoMesh *CreateGMesh3D(int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType) {
//Creating geometric mesh, nodes and elements.
//Including nodes and elements in the mesh object:
    TPZGeoMesh *gmesh2D = CreateGMesh(nelx, nely, hx, hy, x0, y0, meshType);
    TPZExtendGridDimension extend(gmesh2D, hx/nelx);
    TPZGeoMesh *gmesh3D = extend.ExtendedMesh(nelx,-5,-6);
//    delete gmesh2D;
    return gmesh3D;
}

TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel, TPZCompMesh &mesh) {
    if (!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh, gel);

    return NULL;
}

TPZCompMesh *CMesh_S(TPZGeoMesh *gmesh, int pOrder) {
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    cmesh->SetDimModel(dim); //Dimesion of the model

    //Definition of the approximation space:
    cmesh->SetAllCreateFunctionsHDiv(); //Creating H(div) functions:

    //Criando material cujo nSTATE = 2:

    TPZNullMaterial<STATE> * material = new TPZNullMaterial<STATE>(matID);
    material->SetNStateVariables(dim);
    material->SetDimension(dim);
    cmesh->InsertMaterialObject(material); //Insere material na malha



    //Boundary conditions:
    TPZFMatrix<STATE> val1(2, 2, 0.);
    
    TPZManVector<STATE> val2s(2, 0.), val2(2, 0.);
    val2s[0] = 10.0; // vx -> 0
    val2s[1] = 0.0; // vy -> 0

    TPZBndCondT<STATE> * BCond0 = material->CreateBC(material, matBCbott, neumann, val1, val2); //Cria material que implementa a condição de contorno inferior
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha

    auto * BCond1 = material->CreateBC(material, matBCtop, neumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    auto * BCond2 = material->CreateBC(material, matBCleft, neumann, val1, val2s); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha

    auto * BCond3 = material->CreateBC(material, matBCright, neumann, val1, val2s); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha

    auto * BCond4 = material->CreateBC(material, 2, neumann, val1, val2s); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond4); //Insere material na malha

    auto * BCond5 = material->CreateBC(material, -5, neumann, val1, val2s); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond5); //Insere material na malha

    auto * BCond6 = material->CreateBC(material, -6, neumann, val1, val2s); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond6); //Insere material na malha

    val2s[0] = 0;
    val2s[1] = 0;
    cmesh->InsertMaterialObject(material->CreateBC(material, matLagrange, neumann, val1, val2s)); //Insere material na malha

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

TPZCompMesh *CMesh_Girk(TPZGeoMesh *gmesh, int pOrder) {

    //Criando malha computacional:
    int bc_inte_order = 10;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    //    TElasticityExample1 example;

    // Criando material:


    REAL E = 20.59; //* @param E elasticity modulus
    REAL nu = 0.; //* @param nu poisson coefficient


    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = -32.69; //* @param fx forcing function \f$ -x = fx \f$
    int plain = 1.; //* @param plainstress = 1 \f$ indicates use of plainstress

    TPZMixedElasticityMaterialLocal * material1 = new TPZMixedElasticityMaterialLocal(1, E, nu, fx, fy, plain, dim);
    //TPZMixedElasticityMaterial * material2 = new TPZMixedElasticityMaterial(3,E,nu,fx,fy,plain,dim);
    //material1->SetAxisSymmetric();
    //material2->SetAxisSymmetric();
    //material->SetForcingFunction(example.ForcingFunction());
    // Inserindo material na malha
    //    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (f_source);
    //    TPZAutoPointer<TPZFunction<STATE> > pp = new TPZDummyFunction<STATE> (p_exact);
    //    TPZAutoPointer<TPZFunction<STATE> > vp = new TPZDummyFunction<STATE> (v_exact);

    //    material->SetForcingFunction(fp);
    //    material->SetForcingFunctionExactPressure(pp);
    //    material->SetForcingFunctionExact(vp);

    cmesh->InsertMaterialObject(material1);
    //cmesh->InsertMaterialObject(material2);

    //Condições de contorno:

    TPZFMatrix<REAL> val1(2, 2, 0.);
    TPZManVector<REAL,2> val2(2, 0.);
    REAL x;
    val1(1, 1) = material1->BigNumber();

    auto * BCond0 = material1->CreateBC(material1, -3, mixed, val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    //BCond0->SetForcingFunction(solucao_exact,bc_inte_order);
    //BCond0->SetForcingFunction(example.ValueFunction());
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha

    auto * BCond1 = material1->CreateBC(material1, -2, neumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    //BCond1->SetForcingFunction(example.ValueFunction());
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond1->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    val2[1] = 43.533;
    val2.Fill(0.);
    //val1(1,0) = material1->gBigNumber;
    //val1(0,1) = material1->gBigNumber;
    auto * BCond2 = material1->CreateBC(material1, -1, mixed, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(example.ValueFunction());
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    //Cond2->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha

    val2.Fill(0.);
    auto * BCond3 = material1->CreateBC(material1, 2, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //BCond3->SetForcingFunction(example.ValueFunction());
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond3->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha

    //Ponto

    //    TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
    //    val4(0,0)=0.0;
    //
    //    TPZMaterial * BCPoint = material->CreateBC(material, matPoint, pointtype, val3, val4); //Cria material que implementa um ponto para a pressão
    //    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha





    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;


}

TPZCompMesh *CMesh_AxiS(TPZGeoMesh *gmesh, int pOrder, TElasticityExample1 &example) {

    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    //    TElasticityExample1 example;

    // Criando material:

    example.fProblemType = TElasticityExample1::EThiago;
    example.fStressState = TElasticityExample1::EPlaneStrain;
    REAL E = 206.815026; //* @param E elasticity modulus
    REAL nu = 0.; //* @param nu poisson coefStressStatfStressStateficient


    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = -20.; //* @param fx forcing function \f$ -x = fx \f$
    int plain = 0.; //* @param plainstress = 1 \f$ indicates use of plainstress

    TPZMixedElasticityMaterialLocal * material = new TPZMixedElasticityMaterialLocal(matID, E, nu, fx, fy, plain, dim);
    material->SetAxisSymmetric();
    material->SetPlaneStrain();
    //material->SetForcingFunction(example.ForcingFunction());
    // Inserindo material na malha
    //    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (f_source);
    //    TPZAutoPointer<TPZFunction<STATE> > pp = new TPZDummyFunction<STATE> (p_exact);
    //    TPZAutoPointer<TPZFunction<STATE> > vp = new TPZDummyFunction<STATE> (v_exact);

    //    material->SetForcingFunction(fp);
    //    material->SetForcingFunctionExactPressure(pp);
    //    material->SetForcingFunctionExact(vp);

    cmesh->InsertMaterialObject(material);


    //Condições de contorno:

    TPZFMatrix<REAL> val1(2, 2, 0.);
    TPZManVector<REAL,2> val2(2, 0.);

    auto * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    //BCond0->SetForcingFunction(solucao_exact,bc_inte_order);
    //BCond0->SetForcingFunction(example.ValueFunction());
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha

    auto * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    //BCond1->SetForcingFunction(example.ValueFunction());
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond1->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    auto * BCond2 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(example.ValueFunction());
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    //Cond2->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha

    auto * BCond3 = material->CreateBC(material, matBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //BCond3->SetForcingFunction(example.ValueFunction());
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond3->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha

    //Ponto

    //    TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
    //    val4(0,0)=0.0;
    //
    //    TPZMaterial * BCPoint = material->CreateBC(material, matPoint, pointtype, val3, val4); //Cria material que implementa um ponto para a pressão
    //    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha





    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;

}


void InsertMaterialObjects(TPZCompMesh &cmeshref)
{
    // ----------------------- Getting E,nu data --------------------------
    // --------------------------------------------------------------------
    TPZManVector<TPZFMatrix<STATE>,Globny> edata(Globny), nudata(Globny);
    for (int iy = 0; iy < Globny; iy++) {
        edata[iy].Resize(Globnz, Globnx);
        nudata[iy].Resize(Globnz, Globnx);
    }
    REAL tempE = 0., tempNu = 0.;;
    const std::string foldername = "../data/" + to_string(Globnx-1);
    std::string basee = foldername + "/e_", basenu = foldername + "/nu_";
    for (int iy = 0; iy < Globny; iy++) {
        std::string ename = basee + to_string(iy+1) + ".txt";
        std::string nuname = basenu + to_string(iy+1) + ".txt";
        std::ifstream inE(ename), inNu(nuname);
        if(!inE || !inNu) DebugStop(); // Files must exist
        for(int iz = 0; iz < Globnz ; iz++) {
            for(int ix = 0 ; ix < Globnx ; ix++){
                inE >> tempE;
                inNu >> tempNu;
                edata[iy](iz,ix) = tempE;
                nudata[iy](iz,ix) = tempNu;
            }
        }
    }
    
    auto ConstLawFunctionLambda = [edata,nudata] (const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv) {
        int rounded_x = static_cast<int>(floor(x[0]/Globpartsize));
        int rounded_y = static_cast<int>(floor(x[1]/Globpartsize));
        int rounded_z = static_cast<int>(floor(x[2]/Globpartsize));
        STATE eval = edata[rounded_y](rounded_z,rounded_x)/3.e9;
        STATE nuval = nudata[rounded_y](rounded_z,rounded_x);
        result[0] = eval;
        result[1] = nuval;
    };
    
    std::function<void(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)> myfunc = ConstLawFunctionLambda;
    
    
    TPZCompMesh *cmesh = &cmeshref;
    int dim = cmesh->Dimension();
    REAL E = 0.; //* @param E elasticity modulus
    REAL nu = 0.; //* @param nu poisson coefficient


    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = 0.; //* @param fx forcing function \f$ -x = fx \f$
    int plainStress = 0; //* @param plainstress = 1 \f$ indicates use of plainstress
    TPZMixedElasticityND * material = new TPZMixedElasticityND(matID, E, nu, fx, fy, plainStress, dim);
    material->SetElasticityFunction(myfunc);
//    material->SetBigNumber(1.e8);

  
#ifdef LOG4CXX
    if(loggerconfig->isDebugEnabled())
    {
        REAL E,nu;
        material->GetElasticity(E, nu);
        TPZMixedElasticityND::TElasticityAtPoint elast(E,nu);
        TPZFMatrix<STATE> A;
        material->ElasticityModulusTensor(A, elast);
        std::stringstream sout;
        A.Print("A = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerconfig,sout.str())
    }
#endif


    material->SetBodyForce(0., 0.,-9.81/1.e6);
    cmesh->InsertMaterialObject(material);

    //Condições de contorno:

    TPZFMatrix<REAL> val1(dim, dim, 0.);
    TPZManVector<REAL,3> val2(dim, 0.);
    {
        TPZManVector<REAL,3> val2here(dim, 0.);
        val2here[2] = 0.;
        auto * BCond1 = material->CreateBC(material, matBCtop, neumann, val1, val2here); //Cria material que implementa a condicao de contorno superior
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond1);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    {
        auto * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond0);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    
    {
        auto * BCond2 = material->CreateBC(material, matBCleft, neumann, val1, val2); //Cria material que implementa a condicao de contorno esquerda
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond2);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    {
        auto * BCond3 = material->CreateBC(material, matBCright, neumann, val1, val2); //Cria material que implementa a condicao de contorno direita
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond3);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    {
        auto * BCond5 = material->CreateBC(material, matBCfront, neumann, val1, val2); //Cria material que implementa a condicao de contorno direita
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond5);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    {
        auto * BCond6 = material->CreateBC(material, matBCback, neumann, val1, val2); //Cria material que implementa a condicao de contorno direita
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond6);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    {
        auto * BCond4 = material->CreateBC(material, matLagrange, neumann, val1, val2);
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond4);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
}

TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder) {

    //Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();

    // Criando material:
    InsertMaterialObjects(*cmesh);

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
    }
    cmesh->ExpandSolution();
}

static TPZAutoPointer<TPZFunction<STATE> > ConstitutiveLawFunction()
{
    
    TPZAutoPointer<TPZFunction<STATE> > result;
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(TElasticity2DAnalytic::ElasticDummy,4);
    result = TPZAutoPointer<TPZFunction<STATE> >(dummy);
    return result;
}

int main(int argc, char *argv[]) {
//    TPZMaterial::fBigNumber = 1.e16;
        
    // --------------------------

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif 
    bool plotting = false;
	EElementType elementType = ETetraheda;
    int numthreads = 16;

    std::stringstream basename;
    basename << "_EYotov";
    
    std::stringstream rootname;
    rootname << basename.str();
    
    TPZGeoMesh *gmesh = 0;
    std::cout << "\n----------- Creating gmesh -----------" << std::endl;
    const int nrefint = 1;
    if(dim == 2)
    {
        DebugStop();
    }
    else if(dim == 3)
    {
        TPZManVector<REAL,3> minX = {0.,0.,0.};
        TPZManVector<REAL,3> maxX = {10000.,10000.,5000.};
        TPZManVector<int,5> matids = {1,matBCbott,matBCleft,matBCfront,matBCright,matBCback,matBCtop};
        
        TPZManVector<int> nDivs = {(Globnx-1)/(nrefint*2),(Globny-1)/(nrefint*2),(Globnz-1)/(nrefint*2)};
        MMeshType meshType = MMeshType::ETetrahedral;
        rootname << "_tetra_";
        bool createBoundEls = true;
        gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                                                     matids, nDivs, meshType, createBoundEls);
        std::ofstream filegvtk("GMeshInicial.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk, true);
    }
    TPZAutoPointer<TPZMultiphysicsCompMesh> cmesh_m_HDiv;
    const bool isRefine = true;
    if(isRefine){
        TPZCheckGeom check(gmesh);
        check.UniformRefine(nrefint);
        rootname << "Ref1_";
        std::ofstream filegvtk("GMeshInicialRef.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk, true);
    }
    std::cout << "Number of geo elements in mesh = " << gmesh->NElements() << std::endl;
    
    TPZVec<int64_t> coarseindices(gmesh->NElements());
    int64_t count = 0;
    int coarselevel = 0;
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel && gel->Dimension() == gmesh->Dimension() && gel->Level()==coarselevel)
        {
            coarseindices[count++] = el;
        }
    }
    coarseindices.resize(count);
    TPZMHMixedMeshControl control(gmesh);
    // std::cout << "coarse Indices = " << coarseindices <<std::endl;
    control.DefinePartitionbyCoarseIndices(coarseindices);
    
    // Set porder of approximations
    control.SetInternalPOrder(1);
    control.SetSkeletonPOrder(1);
    
    control.fMaterialIds = {1};
    control.fMaterialBCIds = {-1,-2,-3,-4,-5,-6};
    if(dim == 2)
    {
        DebugStop();
    } else if(dim == 3)
    {
        control.SetProblemType(TPZMHMeshControl::EElasticity3D);
    }
    InsertMaterialObjects(control.CMesh());
    //            for(auto bcmatid : control.fMaterialBCIds)
    //            {
    //                control.CMesh()->MaterialVec().erase(bcmatid);
    //            }
    bool substruct = true;
    //            control.DivideSkeletonElements(1);
    
    std::cout << "\n----------- Creating computational mesh -----------" << std::endl;
    control.BuildComputationalMesh(substruct);
    
    rootname << "_Sub";
    cmesh_m_HDiv = control.CMesh();
    
#ifdef PZDEBUG
    {
        std::ofstream fileg("MalhaGeo.txt"); //Prints the geometric mesh in txt format
        std::ofstream filegvtk("MalhaGeo7.vtk"); //Prints the geometric mesh in vtk format
        gmesh->Print(fileg);
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk, true);
    }
#endif
    
    cmesh_m_HDiv->CleanUpUnconnectedNodes();
    cmesh_m_HDiv->SaddlePermute();
    
#ifdef PZDEBUG
    {
        std::ofstream fileg1("MalhaGeo2.txt");
        gmesh->Print(fileg1); //Prints the geometric mesh in txt format
        
        std::ofstream filecm("MalhaC_m.txt");
        cmesh_m_HDiv->Print(filecm); //Prints the multi-physics computational mesh in txt format
        
        //                std::ofstream out("MHMCompMesh.txt");
        //                TPZMHMeshControl::PrintDiagnostics(out);
    }
#endif
    
    //Solving the system:
    bool optimizeBandwidth = true;
    cmesh_m_HDiv->InitializeBlock();
    
    TPZCompMesh *cmesh = cmesh_m_HDiv.operator->();
    TPZLinearAnalysis an(cmesh, optimizeBandwidth); //Creates the object that will manage the analysis of the problem
                                                    // #ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matskl(cmesh);
    // #else
    //             TPZSkylineStructMatrix<STATE> matskl(cmesh); // asymmetric case ***
    // #endif
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    
    std::cout << "\n----------- Assemble -----------" << std::endl;
    std::cout << "Number of equations condensed = " << cmesh->NEquations() << std::endl;
    std::cout << "Number of equations = " << cmesh->Solution().Rows() << std::endl;
    an.Assemble(); //Assembles the global stiffness matrix (and load vector)
    std::cout << "Assemble finished." << std::endl;
    
    TPZMaterial *mat = cmesh_m_HDiv->FindMaterial(matID);
    TPZMatErrorCombinedSpaces<STATE> *materr = dynamic_cast<TPZMatErrorCombinedSpaces<STATE>*>(mat);
    TPZManVector<REAL, 10> Errors(materr->NEvalErrors());
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
    
    std::cout << "\n----------- Solve -----------" << std::endl;
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
    cmesh_m_HDiv->LoadSolutionFromMultiPhysics();
    plotting = true;
    if (plotting) {
#ifdef USENEWVTK
        const std::string plotfile = rootname.str();
        TPZStack<std::string,10> scalnames;
        scalnames.Push("Young_Modulus");
        scalnames.Push("Poisson");
        scalnames.Push("Displacement");
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        scalnames.Push("TauXY");
        
        constexpr int vtkRes{0}; //resolucao do vtk
        auto vtk = TPZVTKGenerator(cmesh_m_HDiv, scalnames, plotfile, vtkRes, 3);
        vtk.SetNThreads(64);
        vtk.Do();
#else
        cout << "\n------------------ Plotting VTK ------------------" << endl;
        std::string plotfile;
        {
            std::stringstream sout;
            sout << rootname.str() << ".vtk";
            plotfile = sout.str();
        }
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        scalnames.Push("TauXY");
        scalnames.Push("Young_Modulus");
        scalnames.Push("Poisson");
        //                vecnames.Push("Flux");
        vecnames.Push("displacement");
        //                vecnames.Push("Stress");
        an.SetStep(1);
        an.DefineGraphMesh(cmesh_m_HDiv->Dimension(), scalnames, vecnames, plotfile);
        an.PostProcess(0);
        
        cout << "\n------------------ Finished Plotting VTK ------------------" << endl;
#endif
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
    
    an.CleanUp();

    std::cout << "FINISHED!" << std::endl;

    return 0;
}
