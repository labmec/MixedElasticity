#include <fstream>

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include <TPZGeoElement.h>

#include "pzmat2dlin.h"
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

#include "TPZMatElasticity2D.h"
#include "TPZInterfaceEl.h"
#include "pzdiscgal.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include <tpzarc3d.h>

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzl2projection.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"
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

#include "pzgengrid.h"
#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
#endif

//------------------Elasticity Problem------------------------

enum EElementType {
    ETriangular = 0, ESquare = 1, ETrapezoidal = 2
};

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
TPZCompMesh *CMesh_AxiS(TPZGeoMesh *gmesh, int pOrder, TElasticityExample1 &example);

/**
 * @brief Funcao para criar a malha computacional multi-fisica ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_Girk(TPZGeoMesh *gmesh, int pOrder);

void CreateCondensedElements(TPZCompMesh *cmesh);

void Error(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);

//Variáveis globais do problema:

const int dim = 2; // Dimension of the problem
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




//void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget) {
//    TPZGeoMesh *gmesh = cmesh.Reference();
//    int64_t nel = gmesh->NElements();
//    for (int64_t el = 0; el < nel; el++) {
//        TPZGeoEl *gel = gmesh->Element(el);
//        if (gel->MaterialId() != matfrom) {
//            continue;
//        }
//
//        int nsides = gel->NSides();
//
//        TPZGeoElSide gelside(gel, nsides - 1);
//        TPZStack<TPZCompElSide> celstack;
//        gelside.EqualLevelCompElementList(celstack, 0, 0);
//        if (celstack.size() != 2) {
//            DebugStop();
//        }
//        gel->SetMaterialId(mattarget);
//        int64_t index;
//        new TPZMultiphysicsInterfaceElement(cmesh, gel, index, celstack[1], celstack[0]);
//    }
//}

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
                    int64_t index;
                    TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(cmesh, gelbc.CreatedElement(), index, celneigh, celside);
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

    TPZGenGrid gengrid(nelem, gcoord1, gcoord2);

    switch (meshType) {
        case ETriangular:
            gengrid.SetElementType(ETriangle);
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
    /*
     int i,j;
     int64_t id, index;
    for(i = 0; i < ny; i++){
        for(j = 0; j < nx; j++){
            id = i*nx + j;
            coord[0] = (x0)+(j)*hx/(nx - 1);
            coord[1] = (y0)+(i)*hy/(ny - 1);
            //          coord[0]=gcoord1[2*i+j];
            //          coord[1]=gcoord2[2*i+j];
            //using the same coordinate x for z
            //            coord[2] = 0.;
            //          newcoord[0]=cos(theta)*coord[0]+sin(theta)*coord[1];
            //          newcoord[1]=sin(theta)*coord[0]-cos(theta)*coord[1];
            //          Get the index in the mesh nodes vector for the new node
            index = gmesh->NodeVec().AllocateNewElement();
            //Set the value of the node in the mesh nodes vector
            gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
        }
    }
    
    
    //  Ponto 1
    //  TPZVec<int64_t> pointtopology(1);
    //  pointtopology[0] = 1;
    //
    //  gmesh->CreateGeoElement(EPoint,pointtopology,matPoint,id);
    
    
    //Vetor auxiliar para armazenar as conecções entre elementos:
    
    TPZVec <int64_t> connect(4,0);
    
    
    //Conectividade dos elementos:
    
    for(i = 0; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            index = (i)*(nx - 1)+ (j);
            connect[0] = (i)*ny + (j);
            connect[1] = connect[0]+1;
            connect[2] = connect[1]+(nx);
            connect[3] = connect[0]+(nx);
            gmesh->CreateGeoElement(EQuadrilateral,connect,matID,id);
            if(id != index) DebugStop();
        }
    }
    
    
    //Gerando informação da vizinhança:
    
    gmesh->BuildConnectivity();
    
    for(i = 0; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            index = (i)*(nx - 1)+ (j);
            TPZGeoEl *gel = gmesh->Element(index);
            if (j==0) {
                TPZGeoElBC gbcl(gel,7,matBCleft);
            }
            if(j == nx-2)
            {
                TPZGeoElBC gbcr(gel,5,matBCright);
            }
            if (i==0) {
                TPZGeoElBC gbcb(gel,4,matBCbott);
            }
            if (i==ny-2) {
                TPZGeoElBC gbct(gel,6,matBCtop);
            }
        }
    }
     */

    //    int64_t idx;
    //
    //    int nx = nelx + 1;
    //    int ny = nely + 1;
    //
    //    TPZVec<int64_t> nodint(2);
    //    for (int i = 0; i < nelx; ++i) {
    //        for (int j = 0; j < nely; ++j) {
    //            if (i != nelx - 1) {
    //                //right edge
    //                nodint[0] = j * nx + i + 1;
    //                nodint[1] = nodint[0] + nx;
    //                gmesh->CreateGeoElement(EOned, nodint, matLagrange, idx);
    //            }
    //            if (j != nely - 1) {
    //                //top edge
    //                nodint[0] = (j + 1) * nx + i;
    //                nodint[1] = nodint[0] + 1;
    //                gmesh->CreateGeoElement(EOned, nodint, matLagrange, idx);
    //            }
    //        }
    //    }
    //
    //    gmesh->BuildConnectivity();

    {
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
    }

    //Printing geometric mesh:

    //ofstream bf("before.vtk");
    //TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;



}

TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index) {
    if (!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh, gel, index);

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
    REAL E = TElasticityExample1::fElast; //* @param E elasticity modulus
    REAL nu = TElasticityExample1::fNu; //* @param nu Poisson coefficient
    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = 0.; //* @param fx forcing function \f$ -x = fx \f$
    int planeStress = 0.; //* @param plainstress = 1 \f$ indicates use of plain stress

    TPZMaterial * material = new TPZMixedElasticityMaterial(matID, E, nu, fx, fy, planeStress, dim);
    cmesh->InsertMaterialObject(material); //Insere material na malha

    //TPZMaterial * material2 = new TPZMixedElasticityMaterial(3,E,nu,fx,fy,plain,dim);
    //cmesh->InsertMaterialObject(material2); //Insere material na malha


    //Dimensões do material (para H1 e descontinuo):
    //TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
    //material->SetMaterial(xkin, xcin, xfin);


    //Boundary conditions:
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
    TPZFMatrix<STATE> val2s(2, 1, 0.);
    val2s(0, 0) = 10.0; // vx -> 0
    val2s(1, 0) = 0.0; // vy -> 0

    TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, neumann, val1, val2); //Cria material que implementa a condição de contorno inferior
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha

    TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, neumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    TPZMaterial * BCond2 = material->CreateBC(material, matBCleft, neumann, val1, val2s); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha

    TPZMaterial * BCond3 = material->CreateBC(material, matBCright, neumann, val1, val2s); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha

    TPZMaterial * BCond4 = material->CreateBC(material, 2, neumann, val1, val2s); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond4); //Insere material na malha

    val2s(0, 0) = 0;
    val2s(1, 0) = 0;
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

        //Criando material cujo nSTATE = 2:
        REAL E = TElasticityExample1::fElast; //* @param E elasticity modulus
        REAL nu = TElasticityExample1::fNu; //* @param nu poisson coefficient
        REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
        REAL fy = -20.; //* @param fx forcing function \f$ -x = fx \f$
        int plain = 0.; //* @param plainstress = 1 \f$ indicates use of plainstress

        TPZMaterial * material = new TPZMixedElasticityMaterial(matID, E, nu, fx, fy, plain, dim);

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
                int64_t index;
                new TPZCompElDiscScaled(*cmesh, gel, index);
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

        //Criando material cujo nSTATE = 2:
        REAL E = TElasticityExample1::fElast; //* @param E elasticity modulus
        REAL nu = TElasticityExample1::fNu; //* @param nu poisson coefficient
        REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
        REAL fy = -20.; //* @param fx forcing function \f$ -x = fx \f$
        int plain = 0.; //* @param plainstress = 1 \f$ indicates use of plainstress

        TPZMaterial * material = new TPZMixedElasticityMaterial(matID, E, nu, fx, fy, plain, dim);

        cmesh->InsertMaterialObject(material); //Insere material na malha

        //TPZMaterial * material2 = new TPZMixedElasticityMaterial(3,E,nu,fx,fy,plain,dim);

        //cmesh->InsertMaterialObject(material2); //Insere material na malha

        //Dimensões do material (para H1 e descontinuo):
        //    TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
        //    material->SetMaterial(xkin, xcin, xfin);

        //Condições de contorno:



        //    //    Ponto de pressao:
        //    //
        //    TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
        //    ////
        //    TPZMaterial * BCPoint = material->CreateBC(material, matPoint, pointtype, val3, val4); //Cria material que implementa um ponto para a pressao
        //    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha



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
    TPZMaterial *material = new TPZMatPoisson3d(matID, dim); //criando material que implementa a formulacao fraca do problema modelo

    cmesh->InsertMaterialObject(material); //Insere material na malha

    //TPZMaterial *material2 = new TPZMatPoisson3d(3,dim);//criando material que implementa a formulacao fraca do problema modelo

    //cmesh->InsertMaterialObject(material2); //Insere material na malha
    //    //Dimensões do material (para H1 e descontínuo):
    //    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    //    material->SetMaterial(xkin, xcin, xfin);

    //Condições de contorno:


    //    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    //
    //        val2(0,0) = 0.0; // px -> 0
    //        val2(1,0) = 0.0; // py -> 0
    //
    //        TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //        cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    //
    //        TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    //        cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    //
    //        TPZMaterial * BCond2 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //        cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    //
    //        TPZMaterial * BCond3 = material->CreateBC(material, matBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //        cmesh->InsertMaterialObject(BCond3); //Insere material na malha


    //    Ponto de pressao:
    //
    //    TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
    //    ////
    //    TPZMaterial * BCPoint = material->CreateBC(material, matPoint, pointtype, val3, val4); //Cria material que implementa um ponto para a pressao
    //    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
    //


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
            int64_t index;
            new TPZCompElDiscScaled(*cmesh, gel, index);
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

    TPZMixedElasticityMaterial * material1 = new TPZMixedElasticityMaterial(1, E, nu, fx, fy, plain, dim);
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

    TPZFMatrix<REAL> val1(2, 2, 0.), val2(2, 1, 0.);
    REAL x;
    val2(0, 0) = 0; // vx -> 0 //val2 represent norm stress;
    val2(1, 0) = 0; // vy -> 0
    val1(0, 0) = 0;
    val1(1, 1) = material1->gBigNumber;

    TPZMaterial * BCond0 = material1->CreateBC(material1, -3, mixed, val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    //BCond0->SetForcingFunction(solucao_exact,bc_inte_order);
    //BCond0->SetForcingFunction(example.ValueFunction());
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha

    TPZMaterial * BCond1 = material1->CreateBC(material1, -2, neumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    //BCond1->SetForcingFunction(example.ValueFunction());
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond1->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    val2(1, 0) = 43.533;
    val2.Zero();
    //val1(1,0) = material1->gBigNumber;
    //val1(0,1) = material1->gBigNumber;
    TPZMaterial * BCond2 = material1->CreateBC(material1, -1, mixed, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(example.ValueFunction());
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    //Cond2->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha

    val2.Zero();
    TPZMaterial * BCond3 = material1->CreateBC(material1, 2, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
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

    TPZMixedElasticityMaterial * material = new TPZMixedElasticityMaterial(matID, E, nu, fx, fy, plain, dim);
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

    TPZFMatrix<REAL> val1(2, 2, 0.), val2(2, 1, 0.);
    val2(0, 0) = 0; // vx -> 0 //val2 represent norm stress;
    val2(1, 0) = 0; // vy -> 0
    //val1(0,0) = 100.0;

    TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    //BCond0->SetForcingFunction(solucao_exact,bc_inte_order);
    //BCond0->SetForcingFunction(example.ValueFunction());
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha

    TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    //BCond1->SetForcingFunction(example.ValueFunction());
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond1->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    TPZMaterial * BCond2 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(example.ValueFunction());
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    //Cond2->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha

    TPZMaterial * BCond3 = material->CreateBC(material, matBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
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

    REAL E = TElasticityExample1::fElast; //* @param E elasticity modulus
    REAL nu = TElasticityExample1::fNu; //* @param nu poisson coefficient


    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = 0.; //* @param fx forcing function \f$ -x = fx \f$
    int plainStress = 0; //* @param plainstress = 1 \f$ indicates use of plainstress
    if (TElasticityExample1::fStressState == TElasticityExample1::EPlaneStrain) {
        plainStress = 0;
    } else if (TElasticityExample1::fStressState == TElasticityExample1::EPlaneStress) {
        plainStress = 1;
    } else if (TElasticityExample1::fStressState == TElasticityExample1::EAxiSymmetric) {
        plainStress = 0;
    } else {
        DebugStop();
    }

    TPZMixedElasticityMaterial * material = new TPZMixedElasticityMaterial(matID, E, nu, fx, fy, plainStress, dim);

    if (TElasticityExample1::fStressState == TElasticityExample1::EAxiSymmetric) {
        material->SetAxisSymmetric();
    }
    TElasticityExample1 example;
    material->SetForcingFunction(example.ForcingFunction());
    // Inserindo material na malha
    //    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (f_source);
    //    TPZAutoPointer<TPZFunction<STATE> > pp = new TPZDummyFunction<STATE> (p_exact);
    //    TPZAutoPointer<TPZFunction<STATE> > vp = new TPZDummyFunction<STATE> (v_exact);

    //    material->SetForcingFunction(fp);
    //    material->SetForcingFunctionExactPressure(pp);
    //    material->SetForcingFunctionExact(vp);

    cmesh->InsertMaterialObject(material);


    //Condições de contorno:

    TPZFMatrix<REAL> val1(2, 2, 0.), val2(2, 1, 0.);
    val2(0, 0) = 0; // vx -> 0 //val2 represent norm stress;
    val2(1, 0) = 0; // vy -> 0
    //val1(0,0) = 100.0;

    TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    //BCond0->SetForcingFunction(solucao_exact,bc_inte_order);
    BCond0->SetForcingFunction(example.ValueFunction());
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha

    TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    BCond1->SetForcingFunction(example.ValueFunction());
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond1->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    TPZMaterial * BCond2 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    BCond2->SetForcingFunction(example.ValueFunction());
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    //Cond2->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha

    TPZMaterial * BCond3 = material->CreateBC(material, matBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond3->SetForcingFunction(example.ValueFunction());
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond3->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha

    TPZMaterial * BCond4 = material->CreateBC(material, matLagrange, neumann, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond4->SetForcingFunction(example.ValueFunction());
    cmesh->InsertMaterialObject(BCond4); //Insere material na malha

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
        int64_t elindex;
        new TPZCompElLagrange(*cmesh, EquationDelay, elindex);
        int64_t groupindex;
        TPZElementGroup *elgr = new TPZElementGroup(*cmesh, groupindex);
        elgr->AddElement(cel);
        elgr->AddElement(cmesh->Element(elindex));
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
    TPZMaterial::gBigNumber = 1.e16;

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    EConfig conf = EThiago;
    int initial_p = 1;
    int final_p = 1;
    int initial_h = 1;
    int final_h = 9;
    bool plotting = false;
    EElementType elementType = ESquare;
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
            TElasticityExample1::fProblemType = TElasticityExample1::EThiago;
            TElasticityExample1::fStressState = TElasticityExample1::EPlaneStrain;
            TElasticityExample1::fElast = 206.8150271873455;
            TElasticityExample1::fNu = 0.3040039545229857;
            hx = 1;
            hy = 1;
            x0 = 0;
            y0 = 0;

            rootname = ConfigRootname[conf] + "_Thiago";
            break;
        case EAxiSymmetric:
        case EAxiSymmetricPlus:
            TElasticityExample1::fProblemType = TElasticityExample1::Etest1;
            TElasticityExample1::fStressState = TElasticityExample1::EAxiSymmetric;
            TElasticityExample1::fElast = 100.;
            TElasticityExample1::fNu = 0.;
            hx = 2;
            hy = 2;
            x0 = 1;
            y0 = -1;
            rootname = ConfigRootname[conf] + "_Test1";
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
            TPZGeoMesh *gmesh = CreateGMesh(nelx, nely, hx, hy, x0, y0, elementType); //Creates the geometric mesh

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


            //TPZCompMesh *cmesh_m_HDiv = CMesh_Girk(gmesh, RibpOrder); //Creates the multi-physics computational mesh
            TPZCompMesh *cmesh_m_HDiv = CMesh_m(gmesh, stressInternalPOrder);
            //TPZCompMesh *cmesh_m_HDiv = CMesh_AxiS(gmesh, InternalpOrder,  Example);
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

            //            AddMultiphysicsInterfaces(*cmesh_m_HDiv);

            //CreateCondensedElements(cmesh_m_HDiv);

#ifdef PZDEBUG
            std::ofstream fileg1("MalhaGeo2.txt");
            gmesh->Print(fileg1); //Prints the geometric mesh in txt format

            std::ofstream filecm("MalhaC_m.txt");
            cmesh_m_HDiv->Print(filecm); //Prints the multi-physics computational mesh in txt format
#endif

            //Solving the system:
            bool optimizeBandwidth = true;
            cmesh_m_HDiv->InitializeBlock();

            TPZCompMesh * cmesh_m_Hybrid;
            TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
            TPZHybridizeHDiv hybridizer;
            tie(cmesh_m_Hybrid, meshvector_Hybrid) = hybridizer.Hybridize(cmesh_m_HDiv, meshvector_HDiv, true, -1.);
            cmesh_m_Hybrid->InitializeBlock();

            TPZAnalysis an(cmesh_m_Hybrid, optimizeBandwidth); //Creates the object that will manage the analysis of the problem
#ifdef USING_MKL
            TPZSymetricSpStructMatrix matskl(cmesh_m_Hybrid);
#else
            TPZSkylineStructMatrix matskl(cmesh_m_Hybrid); // asymmetric case ***
#endif
            matskl.SetNumThreads(numthreads);
            an.SetStructuralMatrix(matskl);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            an.SetSolver(step);

            std::cout << "Assemble matrix with NDoF = " << cmesh_m_Hybrid->NEquations() << "." << std::endl;
            an.Assemble(); //Assembles the global stiffness matrix (and load vector)
            std::cout << "Assemble finished." << std::endl;

            TPZManVector<REAL, 6> Errors(cmesh_m_HDiv->FindMaterial(matID)->NEvalErrors());
            TElasticityExample1 example;
            an.SetExact(example.Exact());
            //            an.PostProcessError(Errors,std::cout);

#ifdef PZDEBUG
            //Imprimir Matriz de rigidez Global:
            if (false) {
                std::ofstream filestiff("stiffness.nb");
                an.Solver().Matrix()->Print("K1 = ", filestiff, EMathematicaInput);

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
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector_Hybrid, cmesh_m_Hybrid);

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
                vecnames.Push("Flux");
                vecnames.Push("displacement");
                vecnames.Push("Stress");
                int count = href * n_ref_p + pref - (initial_p - 1);
                an.SetStep(count);
                an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
                an.PostProcess(2);
            }

#ifdef PZDEBUG
            //Imprimindo vetor solução:
            {
                TPZFMatrix<STATE> solucao = cmesh_m_Hybrid->Solution(); //Pegando o vetor de solução, alphaj
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
            ErroOut << "(* Number of Condensed equations " << cmesh_m_Hybrid->NEquations() << " *)" << std::endl;
            ErroOut << "(* Number of equations before condensation " << cmesh_m_Hybrid->Solution().Rows() << " *)" << std::endl;
            ErroOut << "(*\n";
            an.SetExact(example.Exact());
            an.SetThreadsForError(numthreads);
            bool store_errors = true;
            cmesh_m_Hybrid->ElementSolution().Redim(cmesh_m_Hybrid->NElements(), Errors.size());
            std::cout << "Computing errors." << std::endl;
            an.PostProcessError(Errors, store_errors, ErroOut);
            std::cout << "Computed errors." << std::endl;
            ErroOut << "nelx ribporder internalporder n_condensed - n_total - error_sigma - error_energy - error_div_sigma - error_u - error_r - error_as\n";
            ErroOut << "*)\n";
            TPZManVector<STATE, 10> output(Errors.size() + 5, 0);
            output[0] = h_level;
            output[1] = stressPOrder;
            output[2] = stressInternalPOrder;
            output[3] = cmesh_m_Hybrid->NEquations();
            output[4] = cmesh_m_Hybrid->Solution().Rows();
            for (int i = 0; i < Errors.size(); i++) {
                output[5 + i] = Errors[i];
            }
            ErroOut << "Error[[" << href + 1 << "," << pref + 1 << "]] = {" << output << "};\n";

            std::cout << "Errors = " << Errors << std::endl;


            an.CleanUp();

            delete cmesh_m_Hybrid;
            for (int i = meshvector_Hybrid.size() - 1; i >= 0; i--) {
                meshvector_Hybrid[i]->CleanUp();
                delete meshvector_Hybrid[i];
            }

            delete cmesh_m_HDiv;
            for (int i = meshvector_HDiv.size() - 1; i >= 0; i--) {
                meshvector_HDiv[i]->CleanUp();
                delete meshvector_HDiv[i];
            }
            delete gmesh;

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
