#include <fstream>

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include <TPZGeoElement.h>

#include "pzmat2dlin.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"

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

#include "pzgengrid.h"
#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
#endif

//------------------Problema Elasticidade------------------------


/**
 * @brief Funcao para criar a malha geometrica do problema a ser simulado
 * @note A malha sera unidim5ensional formada por nel elementos de tamanho elsize
 * @param uNDiv number of divisions ortogonal to the plates performed on the domain
 * @param vNDiv number of divisions parallel to the plates performed on the domain
 * @param nel numero de elementos
 * @param elsize tamanho dos elementos
 */
TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy);

/**
 * @brief Funcao para criar a malha computacional da velocidade a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_S(TPZGeoMesh *gmesh, int pOrder);

/// change the order of the internal connect to the given order
void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional da pressão a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_U(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional da pressão a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_P(TPZGeoMesh *gmesh, int pOrder);


/**
 * @brief Funcao para criar a malha computacional multi-fisica ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_Girk(TPZGeoMesh *gmesh, int pOrder);
void CreateCondensedElements(TPZCompMesh *cmesh);


//Variáveis globais do problema:

const int dim = 2; //Dimensão do problema
const int matID = 1; //Materia do elemento volumétrico
const int matID2 = 3; //Materia do elemento volumétrico
const int matMoment = 2;
const int matBCbott = -1, matBCZeroTraction = -2, matBCsymmetry = -3; //Materiais das condições de contorno
const int matPoint = -5; //Materia de um ponto
const int dirichlet = 0, neumann = 1, mixed = 2, pointtype = 5, dirichletvar = 4;

using namespace std;

// Compute the area
REAL AxiArea(TPZGeoMesh * gmesh, std::set<int> matids);
//Função principal do programa:

// integrate r sig.n on the bottom
STATE IntegrateBottom(TPZCompMesh *cmesh, int matid);

// integrate forces along the line between both domains
void IntegrateForceOnInterface(TPZCompMesh *cmesh, TPZVec<STATE> &ForceCartesian, STATE &Shear, STATE &Moment);

// get the equations of a given matid
TPZEquationFilter GetEquationFilter(TPZCompMesh *cmesh, int matid);

int main(int argc, char *argv[]) {

    TPZMaterial::gBigNumber = 1.e16;

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    //Dados do problema:


    int RibpOrder = 2; //Ordem polinomial de aproximação
    int InternalpOrder = 2;
    //double elsizex=hx/nelx, elsizey=hy/nely; //Tamanho dos elementos
    //int nel = elsizex*elsizey; //Número de elementos a serem utilizados

    //Gerando malha geométrica:

    TPZGmshReader Girkmann;

    Girkmann.fPZMaterialId[2]["ELASTICITY1"] = matID;
    Girkmann.fPZMaterialId[2]["ELASTICITY2"] = matID2;
    Girkmann.fPZMaterialId[1]["SUPPORT"] = matBCbott;
    Girkmann.fPZMaterialId[1]["ZERO"] = matBCZeroTraction;
    Girkmann.fPZMaterialId[1]["SYMMETRY"] = matBCsymmetry;
    Girkmann.fPZMaterialId[1]["MOMENT"] = matMoment;
    Girkmann.fPZMaterialId[0]["POINTBC"] = matPoint;

    /*
    Girkmann.fPZMaterialId[2]["ELASTICITY1"] = 1;
    Girkmann.fPZMaterialId[1]["B"] = -1;
    Girkmann.fPZMaterialId[1]["ZERO"] = -2;
    Girkmann.fPZMaterialId[1]["N"] = -3;
     */

#ifdef MACOSX
    TPZGeoMesh *gmesh2 = Girkmann.GeometricGmshMesh("../Girkmann.msh");
#else
    TPZGeoMesh *gmesh2 = Girkmann.GeometricGmshMesh("Girkmann.msh");
#endif
    {
        std::ofstream out("Girkmann.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh2, out, true);
    }

    std::set<int> matids;
    matids.insert(matID);
    REAL Shellvol = AxiArea(gmesh2, matids);
    std::cout << "computed shell volume " << Shellvol << std::endl;
    std::cout << "Shell weight " << Shellvol * 32.69 << std::endl;
    matids.clear();
    matids.insert(matID2);
    REAL Stiffvol = AxiArea(gmesh2, matids);
    std::cout << "Stiffener volume " << Stiffvol << std::endl;
    std::cout << "Stiffener weight " << Stiffvol * 32.69 << std::endl;
    std::cout << "Total Volume " << Stiffvol + Shellvol << std::endl;
    std::cout << "Total weight " << 32.69 * (Stiffvol + Shellvol) << std::endl;
    REAL vol = Stiffvol + Shellvol;
    matids.clear();
    matids.insert(matBCbott);
    REAL bottomarea = AxiArea(gmesh2, matids);
    std::cout << "computed bottom area " << bottomarea << std::endl;
    std::cout << "Reaction force " << bottomarea * 43.5525 << std::endl;
    std::cout << "Theoretical reaction force " << vol * 32.69 / bottomarea << " error " << vol * 32.69 - bottomarea * 43.5525 << std::endl;

    TPZGeoMesh *gmesh = gmesh2;


    TElasticityExample1 Example;
    TPZCompMesh *cmesh_S = CMesh_S(gmesh, RibpOrder); //Função para criar a malha computacional da tensão
    ChangeInternalOrder(cmesh_S, InternalpOrder);
    TPZCompMesh *cmesh_U = CMesh_U(gmesh, InternalpOrder); //Função para criar a malha computacional da deslocamento
    TPZCompMesh *cmesh_P = CMesh_P(gmesh, RibpOrder); //Função para criar a malha computacional da rotação
    //TPZCompMesh *cmesh_m = CMesh_Girktest(gmesh, RibpOrder); //Função para criar a malha computacional multifísica
    TPZCompMesh *cmesh_m = CMesh_Girk(gmesh, InternalpOrder);
#ifdef PZDEBUG
    {
        std::ofstream filecS("MalhaC_S.txt"); //Impressão da malha computacional da tensão (formato txt)
        std::ofstream filecU("MalhaC_U.txt"); //Impressão da malha computacional da deslocamento (formato txt)
        std::ofstream filecP("MalhaC_P.txt"); //Impressão da malha computacional da rotação (formato txt)
        cmesh_S->Print(filecS);
        cmesh_U->Print(filecU);
        cmesh_P->Print(filecP);
    }
#endif

    TPZManVector<TPZCompMesh *, 3> meshvector(3);
    meshvector[0] = cmesh_S;
    meshvector[1] = cmesh_U;
    meshvector[2] = cmesh_P;
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    cmesh_m->LoadReferences();
    //    CreateCondensedElements(cmesh_m);


#ifdef PZDEBUG
    std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
    gmesh->Print(fileg1);

    std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
    cmesh_m->Print(filecm);
#endif

    //Resolvendo o Sistema:
    int numthreads = 0;

    std::cout << "Number of equations " << cmesh_m->NEquations() << std::endl;

    bool optimizeBandwidth = true; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema

#ifdef USING_MKL2
    TPZSymetricSpStructMatrix matskl(cmesh_m);
#else
    TPZSkylineStructMatrix matskl(cmesh_m); //caso nao simetrico ***
#endif

    //    matids.clear();
    //    matids.insert(matID);
    //    matids.insert(matBCsymmetry);
    //    matids.insert(matBCZeroTraction);
    //    matids.insert(matMoment);
    //    matids.insert(matBCbott);
    //    matskl.SetMaterialIds(matids);
    //    matskl.EquationFilter() = GetEquationFilter(cmesh_m, matID);

    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);


    //  std::cout << "Assemble matrix with NDoF = " << cmesh_m->NEquations() << std::endl;

    an.Assemble(); //Assembla a matriz de rigidez (e o vetor de carga) global

    if (0) {
        cmesh_S->Solution().Zero();
        cmesh_P->Solution().Zero();
        cmesh_U->Solution().Zero();
        int64_t nr = cmesh_U->Solution().Rows();
        for (int64_t ir = 1; ir < nr; ir += 2) {
            cmesh_U->Solution()(ir, 0) = 43.5525;
        }
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
        TPZAutoPointer<TPZMatrix<STATE> > stiff = an.Solver().Matrix();
        /**
         * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
         * @param x Is x on the above operation
         * @param y Is y on the above operation
         * @param z Is z on the above operation
         * @param alpha Is alpha on the above operation
         * @param beta Is beta on the above operation
         * @param opt Indicates if is Transpose or not
         */

        TPZFMatrix<STATE> result;
        stiff->MultAdd(cmesh_m->Solution(), an.Rhs(), result, 1., -1.);
        std::cout << "Norm residual " << Norm(result) << std::endl;
        std::set<int> excludematids;
        excludematids.insert(matBCZeroTraction);
        DebugStop();
        //        an.PrintVectorByElement(std::cout, result, excludematids);
    }
#ifdef PZDEBUG0
    //Imprimir Matriz de rigidez Global:
    {
        std::ofstream filestiff("stiffness.txt");
        an.Solver().Matrix()->Print("K1 = ", filestiff, EMathematicaInput);

        std::ofstream filerhs("rhs.txt");
        an.Rhs().Print("R = ", filerhs, EMathematicaInput);
    }
#endif
    an.Solve();

    {
        std::ofstream file("file.txt");
        an.Solution().Print("sol=", file, EMathematicaInput);

    }
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh_m);

    STATE integral = IntegrateBottom(cmesh_m, -1);
    std::cout << "Integrated SigY " << integral << std::endl;

    TPZManVector<STATE, 2> forces(2);
    STATE moment, shear;
    IntegrateForceOnInterface(cmesh_m, forces, shear, moment);

    std::cout << "Integrated forces " << forces << " shear " << shear << " moment " << moment << std::endl;

    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("SigmaX");
    scalnames.Push("SigmaY");
    scalnames.Push("TauXY");
    vecnames.Push("Flux");
    vecnames.Push("displacement");
    vecnames.Push("Stress");
    an.DefineGraphMesh(2, scalnames, vecnames, "mixedelast.vtk");
    an.PostProcess(0);

#ifdef PZDEBUG2
    //Imprimindo vetor solução:
    {
        TPZFMatrix<STATE> solucao = cmesh_m->Solution(); //Pegando o vetor de solução, alphaj
        std::ofstream solout("sol.txt");
        solucao.Print("Sol", solout, EMathematicaInput); //Imprime na formatação do Mathematica

        std::ofstream fileAlpha("alpha.txt");
        an.Solution().Print("Alpha = ", fileAlpha, EMathematicaInput);
    }
#endif



    std::cout << "FINISHED!" << std::endl;

    return 0;
}

TPZCompMesh *CMesh_S(TPZGeoMesh *gmesh, int pOrder) {

    //Criando malha computacional:

    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo


    //Definição do espaço de aprximação:

    cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:


    //Criando material cujo nSTATE = 2:

    REAL E = 1.; //* @param E elasticity modulus
    REAL nu = 0.; //* @param nu poisson coefficient
    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = 0.; //* @param fx forcing function \f$ -x = fx \f$
    int plain = 0.; //* @param plainstress = 1 \f$ indicates use of plainstress

    TPZMaterial * material = new TPZMixedElasticityMaterial(matID, E, nu, fx, fy, plain, dim);
    cmesh->InsertMaterialObject(material); //Insere material na malha
    //material->SetAxisSymmetric();
    TPZMaterial * material2 = new TPZMixedElasticityMaterial(3, E, nu, fx, fy, plain, dim);
    cmesh->InsertMaterialObject(material2); //Insere material na malha


    //Dimensões do material (para H1 e descontinuo):
    //TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
    //material->SetMaterial(xkin, xcin, xfin);


    //Condições de contorno:

    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
    TPZFMatrix<STATE> val2s(2, 1, 0.);
    val2s(0, 0) = 10.0; // vx -> 0
    val2s(1, 0) = 0.0; // vy -> 0

    TPZMat2dLin *matbot = new TPZMat2dLin(matBCbott);
    {
        TPZFNMatrix<4, STATE> xk(2, 2, 0.), xc(2, 2, 0.), xf(2, 1, 0.);
        matbot->SetMaterial(xk, xc, xf);
    }
    cmesh->InsertMaterialObject(matbot); //Insere material na malha

    TPZMaterial * BCond1 = material->CreateBC(material, matBCZeroTraction, neumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    TPZMaterial * BCond2 = material->CreateBC(material, matBCsymmetry, neumann, val1, val2s); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha


    TPZMaterial * BCond4 = material->CreateBC(material, matMoment, neumann, val1, val2s); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond4); //Insere material na malha

    /*
    TPZMaterial * BCond0 = material1->CreateBC(material1, -2, neumann, val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    //BCond0->SetForcingFunction(solucao_exact,bc_inte_order);
    //BCond0->SetForcingFunction(example.ValueFunction());
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    val2(1,0)=-10;
    TPZMaterial * BCond1 = material1->CreateBC(material1, -3, neumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    //BCond1->SetForcingFunction(example.ValueFunction());
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond1->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    //val2(1,0) = 43.533;
    val2(1,0)=0;
    //val1(1,0) = material1->gBigNumber;
    //val1(0,1) = material1->gBigNumber;
    TPZMaterial * BCond2 = material1->CreateBC(material1, -1, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(example.ValueFunction());
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    //Cond2->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
     */
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

    // @omar::

    //Criando malha computacional:

    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo


    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    //Criando material cujo nSTATE = 2:

    REAL E = 1.; //* @param E elasticity modulus
    REAL nu = 0.; //* @param nu poisson coefficient
    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = 0.; //* @param fx forcing function \f$ -x = fx \f$
    int plain = 0.; //* @param plainstress = 1 \f$ indicates use of plainstress

    TPZMaterial * material = new TPZMixedElasticityMaterial(matID, E, nu, fx, fy, plain, dim);

    cmesh->InsertMaterialObject(material); //Insere material na malha

    TPZMaterial * material2 = new TPZMixedElasticityMaterial(matID2, E, nu, fx, fy, plain, dim);
    // material->SetAxisSymmetric();
    cmesh->InsertMaterialObject(material2); //Insere material na malha

    //Dimensões do material (para H1 e descontinuo):
    //    TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
    //    material->SetMaterial(xkin, xcin, xfin);

    //Condições de contorno:

    std::set<int> matids;
    matids.insert(matID);
    matids.insert(matID2);

    cmesh->AutoBuild(matids);


    int ncon = cmesh->NConnects();
    for (int i = 0; i < ncon; i++) {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }


    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha


    //    cmesh->AdjustBoundaryElements();
    //    cmesh->CleanUpUnconnectedNodes();

    return cmesh;

}

TPZCompMesh *CMesh_P(TPZGeoMesh *gmesh, int pOrder) {

    // @omar::

    //pOrder--; // Space restriction apapapa

    //Criando malha computacional:

    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo

    // @omar::
    cmesh->SetAllCreateFunctionsDiscontinuous();

    //    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);


    //Criando material cujo nSTATE = 1:

    TPZMaterial *material = new TPZMatPoisson3d(matID, dim); //criando material que implementa a formulacao fraca do problema modelo

    cmesh->InsertMaterialObject(material); //Insere material na malha

    TPZMaterial *material2 = new TPZMatPoisson3d(matID2, dim); //criando material que implementa a formulacao fraca do problema modelo
    // material->SetAxisSymmetric();
    cmesh->InsertMaterialObject(material2); //Insere material na malha

    std::set<int> materialids;
    materialids.insert(matID);
    materialids.insert(matID2);
    cmesh->AutoBuild(materialids);
    cmesh->LoadReferences();
    //    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    //    cmesh->AutoBuild();


    // @omar::
    int ncon = cmesh->NConnects();
    for (int i = 0; i < ncon; i++) {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }

    int64_t nelem = cmesh->NElements();
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *> (cel);
        if (!disc) {
            continue;
        }
        disc->SetTotalOrderShape();
    }

    //    cmesh->AdjustBoundaryElements();
    //    cmesh->CleanUpUnconnectedNodes();
    cmesh->CleanUpUnconnectedNodes();
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


    REAL E = 20.59e6; //* @param E elasticity modulus
    REAL nu = 0.; //* @param nu poisson coefficient


    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = -32.69; //* @param fx forcing function \f$ -x = fx \f$
    int plain = 1.; //* @param plainstress = 1 \f$ indicates use of plainstress

    TPZMixedElasticityMaterial * material1 = new TPZMixedElasticityMaterial(matID, E, nu, fx, fy, plain, dim);
    TPZMixedElasticityMaterial * material2 = new TPZMixedElasticityMaterial(matID2, E, nu, fx, fy, plain, dim);
    material1->SetAxisSymmetric();
    material2->SetAxisSymmetric();

    cmesh->InsertMaterialObject(material1);
    cmesh->InsertMaterialObject(material2);

    //Condições de contorno:

    TPZFMatrix<REAL> val1(2, 2, 0.), val2(2, 1, 0.);
    REAL x;
    val2(0, 0) = 0; // vx -> 0 //val2 represent norm stress;
    val2(1, 0) = 0; // vy -> 0
    val1(0, 0) = 0;
    val1(1, 1) = material1->gBigNumber;

    // boundary condition for the symmetry line: zero vertical traction
    TPZMaterial * BCond0 = material1->CreateBC(material1, matBCsymmetry, mixed, val1, val2);
    cmesh->InsertMaterialObject(BCond0);

    // zero traction all around
    TPZMaterial * BCond1 = material1->CreateBC(material1, matBCZeroTraction, neumann, val1, val2);
    cmesh->InsertMaterialObject(BCond1);

    val1.Zero();
    val2.Zero();

    // bottom support
    // sliding the in the radial direction
    val1(0, 0) = material1->gBigNumber;
    // very soft spring in the vertical direction
    REAL weight = 10000;
    val1(1, 1) = weight;
    val2.Zero();
    val2(1, 0) = weight * 43.552515;
    //    val2(1,0) = 665.512*weight;
    TPZMaterial * BCond2 = material1->CreateBC(material1, matBCbott, mixed, val1, val2);
    cmesh->InsertMaterialObject(BCond2);

    // fake material - allows to identify the elements along which we should integrate
    val2.Zero();
    TPZMaterial * BCond3 = material1->CreateBC(material1, matMoment, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond3);

    //Ponto

    // add a value on the diagonal of the vertical displament to avoid singular matrix
    TPZFMatrix<REAL> val3(2, 2, 0.), val4(2, 1, 0.);
    val3(1, 1) = 1.0;

    TPZMaterial * BCPoint = material1->CreateBC(material1, matPoint, pointtype, val3, val4);
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

        int64_t newconnectindex = cmesh->AllocateNewConnect(4, 1, 1);
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
        for (int i = 0; i < numel; i++) {
            TPZGeoEl *gel = mphys->Reference();
            TPZCompEl *cel = mphys->Element(i);
            if (!cel) {
                std::cout << "Geometric element " << gel->Index() << " with matid " << gel->MaterialId()
                        << " has no element in approximation space number " << i << std::endl;
                DebugStop();
            }
            numconnects[i] = mphys->Element(i)->NConnects();
        }
        TPZManVector<TPZCompElLagrange::TLagrange, 4> EquationDelay(4);
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
        EquationDelay[3].fConnect[0] = mphys->ConnectIndex(numconnects[0] + numconnects[1] + numconnects[2] - 1);
        EquationDelay[3].fIdf[0] = 0;
        EquationDelay[3].fConnect[1] = newconnectindex;
        EquationDelay[3].fIdf[1] = 3;
        int64_t elindex;
        new TPZCompElLagrange(*cmesh, EquationDelay, elindex);
        int64_t groupindex;
        TPZElementGroup *elgr = new TPZElementGroup(*cmesh, groupindex);
        elgr->AddElement(cel);
        elgr->AddElement(cmesh->Element(elindex));
        TPZCondensedCompEl *condensed = new TPZCondensedCompEl(elgr);
    }
    cmesh->ExpandSolution();
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
        cel->Connect(numconnects[0] + numconnects[1] + numconnects[2] - 1).IncrementElConnected();
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
    STATE integDispy = 0.;
    REAL area = 0.;
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
        int varbound = cel->Material()->VariableIndex("state");
        int varindex = neighbour.Element()->Reference()->Material()->VariableIndex("SigmaY");
        int vardisp = neighbour.Element()->Reference()->Material()->VariableIndex("displacement");
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
            TPZManVector<STATE, 3> sol(1), x(3), soldisp(2), solbound(2);
            neighcel->Solution(volpoint, varindex, sol);
            neighcel->Solution(volpoint, vardisp, soldisp);
            gel->X(locpoint, x);
            cel->Solution(locpoint, varbound, solbound);
            solbound[0] /= x[0];
            solbound[1] /= x[0];
            TPZFNMatrix<4, REAL> jac(1, 1), jacinv(1, 1), axes(1, 3);
            REAL detjac;
            gel->Jacobian(locpoint, jac, axes, detjac, jacinv);
            std::cout << "At x " << x << " SigmaY = " << sol[0] << std::endl;
            std::cout << "solbound = " << solbound << std::endl;
            integSigy -= solbound[1] * weight * x[0] * fabs(detjac);
            integDispy += soldisp[1] * weight * x[0] * fabs(detjac);
            area += weight * x[0] * fabs(detjac);
        }

        delete rule;
    }
    std::cout << "Average vertical displacement " << integDispy / area << std::endl;
    return integSigy;
}

// integrate forces along the line between both domains

void IntegrateForceOnInterface(TPZCompMesh *cmesh, TPZVec<STATE> &ForceCartesian, STATE &Shear, STATE &Moment) {
    int targetmatid = matMoment;
    REAL cosa = cos(2. * M_PI / 9.);
    REAL sina = sin(2. * M_PI / 9.);
    TPZManVector<REAL, 2> normal(2), tangent(2);
    normal[0] = -cosa;
    normal[1] = sina;
    tangent[0] = sina;
    tangent[1] = cosa;
    TPZFNMatrix<4, STATE> tensorL(2, 2), tensorR(2, 2);
    REAL averageR = 15. / sina;
    ForceCartesian.Fill(0.);
    Shear = 0.;
    Moment = 0.;
    REAL ringwidth = 0;
    int64_t nel = cmesh->NElements();
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
        TPZGeoElSide gelside(gel, 2), neighbourL, neighbourR;
        {
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->MaterialId() == matID) {
                    neighbourL = neighbour;
                }
                if (neighbour.Element()->MaterialId() == matID2) {
                    neighbourR = neighbour;
                }
                neighbour = neighbour.Neighbour();
            }
        }
        if (!neighbourL || !neighbourR) {
            DebugStop();
        }
        int varbound = cel->Material()->VariableIndex("state");
        int varindex = neighbourL.Element()->Reference()->Material()->VariableIndex("Stress");
        int vardisp = neighbourL.Element()->Reference()->Material()->VariableIndex("displacement");
        if (neighbourL.Element()->Dimension() != 2 || neighbourL.Element()->Reference() == 0) DebugStop();
        if (neighbourR.Element()->Dimension() != 2 || neighbourR.Element()->Reference() == 0) DebugStop();
        TPZTransform<REAL> trL(1), trR(1);
        gelside.SideTransform3(neighbourL, trL);
        gelside.SideTransform3(neighbourR, trR);
        trL = neighbourL.Element()->SideToSideTransform(neighbourL.Side(), neighbourL.Element()->NSides() - 1).Multiply(trL);
        trR = neighbourR.Element()->SideToSideTransform(neighbourR.Side(), neighbourR.Element()->NSides() - 1).Multiply(trR);
        TPZIntPoints *rule = gel->CreateSideIntegrationRule(2, 7);
        TPZCompEl *neighcelL = neighbourL.Element()->Reference();
        TPZCompEl *neighcelR = neighbourR.Element()->Reference();
        int np = rule->NPoints();
        for (int ip = 0; ip < np; ip++) {
            TPZManVector<REAL, 3> locpoint(1), volpointL(2), volpointR(2), x(3);
            REAL weight;
            rule->Point(ip, locpoint, weight);
            gel->X(locpoint, x);
            trL.Apply(locpoint, volpointL);
            trR.Apply(locpoint, volpointR);
            TPZManVector<STATE, 3> solL(4), soldispL(2), solR(4), soldispR(2), sigbound(2);
            cel->Solution(locpoint, varbound, sigbound);
            sigbound[0] /= x[0];
            sigbound[1] /= x[0];
            neighcelL->Solution(volpointL, varindex, solL);
            neighcelR->Solution(volpointR, varindex, solR);
            neighcelL->Solution(volpointL, vardisp, soldispL);
            neighcelR->Solution(volpointR, vardisp, soldispR);
            TPZFNMatrix<4, REAL> jac(1, 1), jacinv(1, 1), axes(1, 3);
            REAL detjac;
            gel->Jacobian(locpoint, jac, axes, detjac, jacinv);
            ringwidth += fabs(detjac) * weight;
            REAL normax = normal[0] * axes(0, 0) + normal[1] * axes(0, 1);
            std::cout << "normal*axes " << normax << std::endl;
            std::cout << "At x " << x << " Stress Left = " << solL << std::endl;
            std::cout << "At x " << x << " Stress Right = " << solR << std::endl;
            TPZFNMatrix<4, STATE> sigL(2, 2), sigR(2, 2);
            TPZMixedElasticityMaterial::FromVoigt(solL, sigL);
            TPZMixedElasticityMaterial::FromVoigt(solR, sigR);
            TPZManVector<STATE, 2> sigLNormal(2), sigRNormal(2);
            sigLNormal[0] = sigL(0, 0) * normal[0] + sigL(0, 1) * normal[1];
            sigLNormal[1] = sigL(1, 0) * normal[0] + sigL(1, 1) * normal[1];
            sigRNormal[0] = sigR(0, 0) * normal[0] + sigR(0, 1) * normal[1];
            sigRNormal[1] = sigR(1, 0) * normal[0] + sigR(1, 1) * normal[1];
            STATE diffnormal = sqrt((sigLNormal[0] - sigRNormal[0])*(sigLNormal[0] - sigRNormal[0])+(sigLNormal[1] - sigRNormal[1])*(sigLNormal[1] - sigRNormal[1]));
            std::cout << "sigLNormal " << sigLNormal << " sigRNormal " << sigRNormal << " diff normal " << diffnormal << std::endl;
            std::cout << "sig Normal by boundary " << sigbound << std::endl;
            ForceCartesian[0] += sigbound[0] * weight * x[0] * fabs(detjac);
            ForceCartesian[1] += sigbound[1] * weight * x[0] * fabs(detjac);
            STATE sigNN = sigbound[0] * normal[0] + sigbound[1] * normal[1];
            STATE sigNT = sigbound[0] * tangent[0] + sigbound[1] * tangent[1];
            REAL radius = sqrt(x[0] * x[0] + x[1] * x[1]);
            Shear += sigNT * x[0] * weight * fabs(detjac);
            Moment += sigNN * x[0] * weight * fabs(detjac)*(radius - averageR);
        }

        delete rule;
    }
    //    Shear *= 2.*M_PI;
    //    Moment *= 2.*M_PI;
    std::cout << "ringwidth = " << ringwidth << std::endl;
}

// get the equations of a given matid

TPZEquationFilter GetEquationFilter(TPZCompMesh *cmesh, int matid) {
    int64_t neq = cmesh->NEquations();
    int64_t ncon = cmesh->NConnects();
    TPZVec<int> numelcon(ncon, 0);
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if (gel->MaterialId() == matid) {
            int ncon = cel->NConnects();
            for (int ic = 0; ic < ncon; ic++) {
                int64_t cindex = cel->ConnectIndex(ic);
                numelcon[cindex]++;
            }
        }
    }
    TPZStack<int64_t> active;
    for (int64_t ic = 0; ic < ncon; ic++) {
        if (numelcon[ic] > 0) {
            int64_t seqnum = cmesh->ConnectVec()[ic].SequenceNumber();
            int blsize = cmesh->ConnectVec()[ic].NDof();
            int64_t pos = cmesh->Block().Position(seqnum);
            for (int idf = 0; idf < blsize; idf++) {
                active.Push(pos + idf);
            }
        }
    }
    TPZEquationFilter filter(neq);
    filter.SetActiveEquations(active);
    return filter;

}
