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

#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
static LoggerPtr loggerconfig(Logger::getLogger("pz.hdiv.vecconfig"));
#endif

TPZAnalyticSolution *gAnalytic = 0;
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

void CreateCondensedElements(TPZCompMesh *cmesh);

void Error(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);

void InsertMaterialObjects(TPZCompMesh &cmeshref);

TPZGeoMesh* Create3DTetMesh(const int href);

//Variáveis globais do problema:

const int dim = 3; // Dimension of the problem
const int matID = 1; // Material of the volumetric element
const int matLagrange = -10; // Material of the Lagrange multipliers
const int matBCbott = -1, matBCtop = -2, matBCleft = -3, matBCright = -4; // Materials of the boundary conditions
const int dirichlet = 0, neumann = 1, mixed = 2, dirichletvar = 4, pointtype = 5; // Boundary conditions of the problem ->default: Dirichlet on left and right

using namespace std;

// Compute the area
REAL AxiArea(TPZGeoMesh * gmesh, std::set<int> matids);
//Função principal do programa:

// integrate r sig.n on the bottom
STATE IntegrateBottom(TPZCompMesh *cmesh, int matid);

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





void AddMultiphysicsInterfaces(TPZCompMesh &cmesh) {
    TPZGeoMesh *gmesh = cmesh.Reference();
    std::set<int> velmatid;
    velmatid.insert(matLagrange);
    velmatid.insert(matBCtop);
    velmatid.insert(matBCbott);
    velmatid.insert(matBCleft);
    velmatid.insert(matBCright);

    // volumetric element variables to be exported to the interface: stress tensor
    TPZManVector<int64_t, 4> LeftElIndices(4);
    LeftElIndices[0] = 0;
    LeftElIndices[1] = 1;
    LeftElIndices[2] = 2;
    LeftElIndices[3] = 3;
    // 1-D element variables to be exported to the interface: Lagrange multipliers
    TPZManVector<int64_t, 2> RightElIndices(2);
    RightElIndices[0] = 0;
    RightElIndices[1] = 1;

    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel) {
            continue;
        }
        int matid = gel->MaterialId();
        if (velmatid.find(matid) != velmatid.end()) {
            int nsides = gel->NSides();
            TPZGeoElSide gelside(gel, nsides - 1);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {

                if (neighbour.Element()->Dimension() == 2 && gelside.Element()->MaterialId() == matLagrange) {
                    // create an interface element
                    TPZCompElSide celside = gelside.Reference();
                    TPZCompElSide celneigh = neighbour.Reference();
                    if (!celside || !celneigh) {
                        DebugStop();
                    }
                    std::cout << "Created an element between volumetric element " << neighbour.Element()->Index() <<
                            " side " << neighbour.Side() <<
                            " and interface element " << gelside.Element()->Index() << std::endl;
                    TPZGeoElBC gelbc(gelside, matID);
                    TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(cmesh, gelbc.CreatedElement(), celneigh, celside);
                    intf->SetLeftRightElementIndices(LeftElIndices, RightElIndices);
                }
                neighbour = neighbour.Neighbour();
            }
        }
    }
}

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
    TPZCompMesh *cmesh = &cmeshref;
    int dim = cmesh->Dimension();
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
    material->SetForcingFunction(gAnalytic->ForceFunc(),5);
    // Inserindo material na malha
    //    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (f_source);
    //    TPZAutoPointer<TPZFunction<STATE> > pp = new TPZDummyFunction<STATE> (p_exact);
    //    TPZAutoPointer<TPZFunction<STATE> > vp = new TPZDummyFunction<STATE> (v_exact);

    //    material->SetForcingFunction(fp);
    //    material->SetForcingFunctionExactPressure(pp);
    //    material->SetForcingFunctionExact(vp);

    cmesh->InsertMaterialObject(material);


    //Condições de contorno:

    TPZFMatrix<REAL> val1(dim, dim, 0.);
    TPZManVector<REAL,3> val2(dim, 0.);


    {
        auto * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
        BCond1->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond1);
        //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
        //BCond1->SetForcingFunction(solucao_exact,bc_inte_order);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    {
        auto * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond0);
        BCond0->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    
    {
        auto * BCond2 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond2);
        BCond2->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
        //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
        //Cond2->SetForcingFunction(solucao_exact,bc_inte_order);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    {
        auto * BCond3 = material->CreateBC(material, matBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond3);
        BCond3->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    {
        auto * BCond5 = material->CreateBC(material, -5, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond5);
        BCond5->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    {
        auto * BCond6 = material->CreateBC(material, -6, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond6);
        BCond6->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    {
        auto * BCond4 = material->CreateBC(material, matLagrange, neumann, val1, val2); //Cria material que implementa a condicao de contorno direita
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond4);
        BCond4->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }
    //Ponto

    //    TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
    //    val4(0,0)=0.0;
    //
    //    TPZMaterial * BCPoint = material->CreateBC(material, matPoint, pointtype, val3, val4); //Cria material que implementa um ponto para a pressão
    //    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha

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


    InsertMaterialObjects(*cmesh);



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
        auto ellagrange = new TPZCompElLagrange(*cmesh, EquationDelay);
        TPZElementGroup *elgr = new TPZElementGroup(*cmesh);
        elgr->AddElement(cel);
        elgr->AddElement(ellagrange);
        TPZCondensedCompEl *condensed = new TPZCondensedCompEl(elgr, false);
    }
    cmesh->InitializeBlock();
}

void CreateCondensedElements2(TPZCompMesh *cmesh) {
    int64_t nel = cmesh->NElements();
    cmesh->ComputeNodElCon();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *> (cel);
        if (!mphys) continue;
        TPZGeoEl *gel = mphys->Reference();
        if (gel->Dimension() != cmesh->Dimension()) continue;

        int numel = mphys->NumberOfCompElementsInsideThisCompEl();
        if (numel != 3) DebugStop();
        TPZManVector<int> numconnects(numel, 0);
        for (int i = 0; i < numel; i++) numconnects[i] = mphys->Element(i)->NConnects();
        cel->Connect(numconnects[0]).IncrementElConnected();
        cel->Connect(numconnects[0] + 1).IncrementElConnected();
        cel->Connect(numconnects[0] + numconnects[1] + gel->NCornerNodes() - 1).IncrementElConnected();
        TPZCondensedCompEl *condensed = new TPZCondensedCompEl(cel);
    }
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

// Compute the area

REAL AxiArea(TPZGeoMesh * gmesh, std::set<int> matids) {
    int64_t nel = gmesh->NElements();
    REAL result = 0.;
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int matid = gel->MaterialId();
        if (matids.find(matid) == matids.end()) continue;
        int nsides = gel->NSides();
        int order = 1;
        TPZIntPoints *intpoints = gel->CreateSideIntegrationRule(nsides - 1, order);
        int np = intpoints->NPoints();
        TPZManVector<REAL, 3> point(2, 0.), x(3, 0.);
        REAL weight;

        /** @brief Compute a decomposition of the gradient of the mapping function, as a rotation matrix (Jacobian) and orthonormal basis (axes)  **/
        //        void Jacobian(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const;


        TPZFNMatrix<6, REAL> jac(2, 2), jacinv(2, 2), axes(2, 3);
        REAL detjac;
        //         TPZManVector<REAL,3> co(3);
        //         std::cout << " gel index "  << el << std::endl;
        //         for(int n=0; n<3; n++)
        //         {
        //             gel->NodePtr(n)->GetCoordinates(co);
        //             std::cout << co << " \n ";
        //         }
        REAL elarea = 0.;
        for (int ip = 0; ip < np; ip++) {
            //        virtual void Point(int i, TPZVec<REAL> &pos, REAL &w) const = 0;
            intpoints->Point(ip, point, weight);
            gel->Jacobian(point, jac, axes, detjac, jacinv);
            gel->X(point, x);
            elarea += abs(detjac) * weight * x[0];
        }
        //         std::cout << "elarea " << elarea << std::endl;
        result += elarea;
        delete intpoints;
    }
    return result;
}

// integrate r sig.n on the bottom

STATE IntegrateBottom(TPZCompMesh *cmesh, int targetmatid) {
    int64_t nel = cmesh->NElements();
    STATE integSigy = 0.;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) continue;
        int matid = gel->MaterialId();
        if (targetmatid != matid) {
            continue;
        }
        if (gel->Dimension() != 1) DebugStop();
        TPZGeoElSide gelside(gel, 2);
        TPZGeoElSide neighbour = gelside.Neighbour();
        int varindex = neighbour.Element()->Reference()->Material()->VariableIndex("SigmaY");
        if (neighbour.Element()->Dimension() != 2 || neighbour.Element()->Reference() == 0) DebugStop();
        TPZTransform<REAL> tr(1);
        gelside.SideTransform3(neighbour, tr);
        tr = neighbour.Element()->SideToSideTransform(neighbour.Side(), neighbour.Element()->NSides() - 1).Multiply(tr);
        TPZIntPoints *rule = gel->CreateSideIntegrationRule(2, 7);
        TPZCompEl *neighcel = neighbour.Element()->Reference();
        int np = rule->NPoints();
        for (int ip = 0; ip < np; ip++) {
            TPZManVector<REAL, 3> locpoint(1), volpoint(2);
            REAL weight;
            rule->Point(ip, locpoint, weight);
            tr.Apply(locpoint, volpoint);
            TPZManVector<STATE, 3> sol(1), x(3);
            neighcel->Solution(volpoint, varindex, sol);
            gel->X(locpoint, x);
            integSigy += sol[0] * weight * x[0];
        }

        delete rule;
    }
    return integSigy;
}

int main(int argc, char *argv[]) {
//    TPZMaterial::gBigNumber = 1.e16;

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif 
    EConfig conf = EThiago;
    int initial_p = 1;
    int final_p = 1;
    int initial_h = 0;
    int final_h = 0;
    bool plotting = false;
//    EElementType elementType = ESquare;
	EElementType elementType = ETetraheda;
    int numthreads = 8;

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

//#ifdef PZ_USING_MKL
//    mkl_set_dynamic(0); // disable automatic adjustment of the number of threads
//    mkl_set_num_threads(numthreads);
//#endif

    std::stringstream basename;
    double hx = 2, hy = 2; //Dimensões em x e y do domínio
    double x0 = -1;
    double y0 = -1;

    //Problem data:
    basename << ConfigRootname[conf];
    switch (conf) {
        case EThiago:
        case EThiagoPlus:
        case EThiagoPlusPlus:
            if(dim == 2)
            {
                TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
                elas->gE = 206.8150271873455;
                elas->gPoisson = 0.3040039545229857;
                elas->fProblemType = TElasticity2DAnalytic::ERot;
                elas->fPlaneStress = 0;
                gAnalytic = elas;
            }
            else if(dim == 3)
            {
                TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
                elas->fE = 250.;//206.8150271873455;
                elas->fPoisson = 0.25;//0.3040039545229857;
                elas->fProblemType = TElasticity3DAnalytic::EYotov;
                basename << ConfigRootname[conf] << "_EYotov";
                gAnalytic = elas;
            }
            else
                DebugStop();
            hx = 1;
            hy = 1;
            x0 = 0;
            y0 = 0;

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
            basename << "_Test1";
        }
            break;
        default:
            DebugStop();
            break;
    }
    //    TPZManVector<STATE, 2> force(2);
    //    TPZFNMatrix<4, STATE> sigma(2, 2);
    //    TPZManVector<REAL, 3> x(3, 0.);
    //    x[0] = x0 + hx / 2.;
    //    x[1] = y0 + hy / 2.;
    //    TElasticityExample1::Sigma(x, sigma);
    //    TElasticityExample1::Force(x, force);

    for (unsigned int pref = initial_p - 1; pref < final_p; ++pref) {
        for (unsigned int href = initial_h; href <= final_h; ++href) {
            std::stringstream rootname;
            rootname << basename.str();
            unsigned int h_level = 1 << href;
            int nelx = h_level, nely = h_level; //Number of elements in x and y directions
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
            TPZGeoMesh *gmesh = 0;
            if(dim == 2)
            {
                gmesh = CreateGMesh(nelx, nely, hx, hy, x0, y0, elementType); //Creates the geometric mesh
            }
            else if(dim == 3)
            {
//				if (elementType == ESquare)
//					gmesh = CreateGMesh3D(nelx, nely, hx, hy, x0, y0, elementType); //Creates the geometric mesh
//				else if (elementType == ETetraheda)
//					gmesh = Create3DTetMesh(href);
				TPZManVector<REAL,3> minX = {0.,0.,0.};
				TPZManVector<REAL,3> maxX = {1.,1.,1.};
				TPZManVector<int,5> matids = {1,-1,-1,-1,-1,-1,-1};
				TPZManVector<int> nDivs = {nelx,nelx,nelx};
				MMeshType meshType = MMeshType::ETetrahedral;
				rootname << "_tetra_";
				bool createBoundEls = true;
				gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
						matids, nDivs, meshType, createBoundEls);
                std::ofstream filegvtk("GMeshInicial.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk, true);
            }
            TPZAutoPointer<TPZMultiphysicsCompMesh> cmesh_m_HDiv;
            if(1){
                TPZCheckGeom check(gmesh);
                check.UniformRefine(0);
                rootname << "Ref1_";
                std::ofstream filegvtk("GMeshInicialRef.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk, true);
            }
            
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
            std::cout << "coarse Indices = " << coarseindices <<std::endl;
            control.DefinePartitionbyCoarseIndices(coarseindices);
            
            control.SetInternalPOrder(pref+1);
            control.SetSkeletonPOrder(pref+1);

            control.fMaterialIds = {1};
            control.fMaterialBCIds = {-1,-2,-3,-4,-5,-6};
            if(dim == 2)
            {
                control.SetProblemType(TPZMHMeshControl::EElasticity2D);
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
//            TPZCompMesh *cmesh_P_HDiv = CMesh_P(gmesh, rotationPOrder, hx / nelx); //Creates the computational mesh for the rotation field

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

            std::cout << "Assemble matrix with NDoF = " << cmesh->NEquations() << std::endl;
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
            cmesh_m_HDiv->LoadSolutionFromMultiPhysics();
            plotting = true;
            if (plotting) {
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
//                vecnames.Push("Flux");
                vecnames.Push("displacement");
//                vecnames.Push("Stress");
                int count = href * n_ref_p + pref - (initial_p - 1);
                an.SetStep(count);
                an.DefineGraphMesh(cmesh_m_HDiv->Dimension(), scalnames, vecnames, plotfile);
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
            sout << rootname.str();
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
            ErroOut << "(* Type of simulation " << rootname.str() << " *)\n";
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
            cmesh->ExpandElementSolution(Errors.size());
//            cmesh->ElementSolution().Redim(cmesh->NElements(), Errors.size());
            std::cout << "Computing errors." << std::endl;
            Errors.Fill(0.);
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
        }
    }

    //
    //
    //
    //    //Pós-processamento (paraview):
    //    std::cout << "Post Processing " << std::endl;
    //    std::string plotfile("ElasticityTest.vtk");
    //    TPZStack<std::string> scalnames, vecnames;
    //    vecnames.Push("Displacement");
    //    vecnames.Push("Stress");
    //    vecnames.Push("Rotation");
    ////    vecnames.Push("V_exact");
    ////    vecnames.Push("P_exact");
    //    //        vecnames.Push("V_exactBC");
    //
    //
    //    int postProcessResolution = 3; //  keep low as possible
    //
    //    int dim = gmesh->Dimension();
    //    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    //    an.PostProcess(postProcessResolution,dim);

    std::cout << "FINISHED!" << std::endl;

    return 0;
}

TPZGeoMesh* Create3DTetMesh(const int href) {
	TPZManVector<REAL,3> minX = {0.,0.,0.};
	TPZManVector<REAL,3> maxX = {1.,1.,1.};
	TPZManVector<int,5> matids = {1,-1,-1,-1,-1,-1,-1};
	//	int signedhref = (int) href;
	TPZManVector<int> nDivs = {href,href,href};
	MMeshType tetmeshType = MMeshType::ETetrahedral;
	bool createBoundEls = true;
	return TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
												matids, nDivs,
												tetmeshType, createBoundEls);
}
