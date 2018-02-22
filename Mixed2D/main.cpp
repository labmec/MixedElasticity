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
TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy, double x0, double y0);

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

const int dim = 2; //Dimensão do problema
const int matID = 1; //Materia do elemento volumétrico
const int matBCbott = -1, matBCtop = -2, matBCleft = -3, matBCright = -4; //Materiais das condições de contorno
const int dirichlet = 0, neumann = 1, mixed = 2, pointtype=5, dirichletvar=4; //Condições de contorno do problema ->default Dirichlet na esquerda e na direita



using namespace std;

// Compute the area
REAL AxiArea(TPZGeoMesh * gmesh, std::set<int> matids);
//Função principal do programa:

// integrate r sig.n on the bottom
STATE IntegrateBottom(TPZCompMesh *cmesh, int matid);

enum EConfig {EThiago, EAxiSymmetric, EThiagoPlus};

std::string ConfigRootname[3] = {
    "Mixed",
    "Mixed_AxiSymmetric",
    "MixedPlus"
};

int main(int argc, char *argv[])
{
    
    TPZMaterial::gBigNumber = 1.e16;
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    EConfig conf = EAxiSymmetric;
    int maxrefp = 1;
    int maxrefh = 7;
    bool plotting = false;
    
    
    std::string rootname;
    double hx=2,hy=2; //Dimensões em x e y do domínio
    double x0 = -1;
    double y0 = -1;


    //Dados do problema:
    switch (conf) {
        case EThiago:
            TElasticityExample1::fProblemType = TElasticityExample1::EPoly;
            TElasticityExample1::fStressState   = TElasticityExample1::EPlaneStrain;
            TElasticityExample1::fElast = 206.815026;
            TElasticityExample1::fNu = 0.30400395;
            hx = 2;
            hy = 2;
            x0 = -1;
            y0 = -1;
            
            rootname = "../" + ConfigRootname[0] + "_Poly";
            break;
        case EAxiSymmetric:
            TElasticityExample1::fProblemType = TElasticityExample1::Etest1;
            TElasticityExample1::fStressState = TElasticityExample1::EAxiSymmetric;
            TElasticityExample1::fElast = 100.;
            TElasticityExample1::fNu = 0.;
            hx = 2;
            hy = 2;
            x0 = 1;
            y0 = -1;
            rootname = "../" + ConfigRootname[1] + "_Test1";
            break;
        case EThiagoPlus:
            TElasticityExample1::fProblemType = TElasticityExample1::EThiago;
            TElasticityExample1::fStressState   = TElasticityExample1::EPlaneStrain;
            TElasticityExample1::fElast = 206.815026;
            TElasticityExample1::fNu = 0.30400395;
            hx = 2;
            hy = 2;
            x0 = -1;
            y0 = -1;
            rootname = "../" + ConfigRootname[2]+"_Thiago";
            break;
        default:
            DebugStop();
            break;
    }
    TPZManVector<STATE,2> displ(2), force(2);
    TPZFNMatrix<4,STATE> sigma(2,2);
    TPZManVector<REAL,3> x(3,0.);
    x[0]=x0+hx/2.;
    x[1]=y0+hy/2.;
    TElasticityExample1::Sigma(x,sigma);
    TElasticityExample1::Force(x,force);
    
    for (int iref = 0; iref < maxrefh; iref++) {
        for (int pref = 0; pref < maxrefp; pref++)
        {
            int h_level = 1 << iref;
            int nelx=h_level, nely=h_level; //Número de elementos em x e y
            int nx=nelx+1 ,ny=nely+1; //Número de nos em x  y
            int RibpOrder = pref+1; //Ordem polinomial de aproximação
            int InternalpOrder = pref+1;
            if (conf == EThiagoPlus) {
                InternalpOrder = pref+2;
            }
            TPZGeoMesh *gmesh = CreateGMesh(nx, ny, hx, hy, x0, y0); //Função para criar a malha geometrica
            
#ifdef PZDEBUG
            std::ofstream fileg("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
            std::ofstream filegvtk("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
            gmesh->Print(fileg);
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
            
            
            //Gerando malha computacional:
            //   TPZGeoMesh *gmesh=CreateGMesh(nx,ny,hx,hy);
            
            TPZCompMesh *cmesh_S = CMesh_S(gmesh, RibpOrder); //Função para criar a malha computacional da tensão
            ChangeInternalOrder(cmesh_S,InternalpOrder);
            TPZCompMesh *cmesh_U = CMesh_U(gmesh, InternalpOrder); //Função para criar a malha computacional da deslocamento
            TPZCompMesh *cmesh_P = CMesh_P(gmesh, RibpOrder); //Função para criar a malha computacional da rotação
                                                                   //TPZCompMesh *cmesh_m = CMesh_Girk(gmesh, RibpOrder); //Função para criar a malha computacional multifísica
            TPZCompMesh *cmesh_m = CMesh_m(gmesh, InternalpOrder);
            //TPZCompMesh *cmesh_m = CMesh_AxiS(gmesh, InternalpOrder,  Example);
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
            
            TPZManVector<TPZCompMesh *, 2> meshvector(3);
            meshvector[0] = cmesh_S;
            meshvector[1] = cmesh_U;
            meshvector[2] = cmesh_P;
            TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
            TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
            TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
            cmesh_m->LoadReferences();
            CreateCondensedElements(cmesh_m);
            
            
#ifdef PZDEBUG
            std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
            gmesh->Print(fileg1);
            
            std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
            cmesh_m->Print(filecm);
#endif
            
            //Resolvendo o Sistema:
            int numthreads = 0;
            
            bool optimizeBandwidth = true;
            TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
            TPZSkylineStructMatrix matskl(cmesh_m); //caso nao simetrico ***
            matskl.SetNumThreads(numthreads);
            an.SetStructuralMatrix(matskl);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            an.SetSolver(step);
            
            
            //  std::cout << "Assemble matrix with NDoF = " << cmesh_m->NEquations() << std::endl;
            
            an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
            
            TPZManVector<REAL,3> Errors;
            TElasticityExample1 example;
            an.SetExact(example.Exact());
//            an.PostProcessError(Errors,std::cout);
            
#ifdef PZDEBUG
            //Imprimir Matriz de rigidez Global:
            {
                std::ofstream filestiff("stiffness.nb");
                an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
                
                std::ofstream filerhs("rhs.nb");
                an.Rhs().Print("R = ",filerhs,EMathematicaInput);
            }
#endif
            /*
             REAL sumrhs = 0.;
             TPZFMatrix<STATE> &rhs = an.Rhs();
             for(long i=0; i< rhs.Rows(); i++)
             {
             sumrhs += rhs(i,0);
             }
             std::cout << " sumrhs "  << sumrhs << std::endl;
             std::cout << "Solving Matrix " << std::endl;
             */
            an.Solve();
            
#ifdef PZDEBUG
            {
                std::ofstream file("file.txt");
                an.Solution().Print("sol=",file,EMathematicaInput);
                
            }
#endif
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector,cmesh_m);

            if(plotting)
            {
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
                int count = iref*maxrefp+pref;
                an.SetStep(count);
                an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
                an.PostProcess(2);
            }
#ifdef PZDEBUG
            //Imprimindo vetor solução:
            {
                TPZFMatrix<STATE> solucao=cmesh_m->Solution();//Pegando o vetor de solução, alphaj
                std::ofstream solout("sol.nb");
                solucao.Print("Sol",solout,EMathematicaInput);//Imprime na formatação do Mathematica
                
                std::ofstream fileAlpha("alpha.nb");
                an.Solution().Print("Alpha = ",fileAlpha,EMathematicaInput);
            }
#endif
            
            //   matids.clear();
            //   matids.insert(-1);
            //   TPZManVector<STATE,3> result;
            //  result = cmesh_m->Integrate("state",matids);
            //  std::cout << "Sigma Y"  << result << std::endl;
            
            
            //    //Calculo do erro
            //    std::cout << "Computing Error " << std::endl;
            
            std::stringstream sout;
            sout << rootname << "Error.nb";
            ofstream ErroOut(sout.str(),std::ios::app);
            ErroOut << "(* Number of elements " << h_level << "*)" << std::endl;
            ErroOut << "(* Number of Condensed equations " << cmesh_m->NEquations() << "*)" << std::endl;
            ErroOut << "(* Number of equations before condensation " << cmesh_m->Solution().Rows() << "*)" << std::endl;
            ErroOut << "(*\n";
            an.SetExact(example.Exact());
            an.PostProcessError(Errors,ErroOut);
            ErroOut << "nelx ribporder internalporder - error_u - error_energy - error_sigma\n";
            ErroOut << "*)\n";
            TPZManVector<STATE,10> output(Errors.size()+5,0);
            output[0] = h_level;
            output[1] = pref+1;
            output[2] = InternalpOrder;
            output[3] = cmesh_m->NEquations();
            output[4] = cmesh_m->Solution().Rows();
            for (int i=0; i<Errors.size(); i++) {
                output[5+i] = Errors[i];
            }
            ErroOut << "Error[[" << iref+1 << "," << pref+1 << "]] = {" << output << "};";
            
            std::cout << "Errors = " << Errors << std::endl;
            
            
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


TPZGeoMesh *CreateGMesh(int nx, int ny, double hx, double hy, double x0, double y0)
{
    
    int i,j;
    long id, index;
    
    
    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    
    //Vetor auxiliar para armazenar coordenadas:
    
    TPZVec <REAL> coord (3,0.);
    //    TPZVec <REAL> newcoord(3,0.);
    //    TPZVec <REAL> gcoord1(3,0.);
    //    TPZVec <REAL> gcoord2(3,0.);
    //    double theta=M_PI/4;
    //    gcoord1[0]=-1;
    //    gcoord2[0]=hx-1;
    //    gcoord1[1]=-1.;
    //    gcoord2[1]=hy-1;
    //    gcoord1[2]=0-1;
    //   gcoord2[2]=0-1;
    //Inicialização dos nós:
    
    // TPZManVector<int> nelem(2,1);
    //  nelem[0]= nx;
    // nelem[1] = ny;
    
    // TPZGenGrid gengrid(nelem,gcoord1,gcoord2);
    
    //    gengrid.SetElementType(ETriangle);
    
    REAL distort = 0.7;
    //    gengrid.SetDistortion(distort);
    
    //gengrid.Read(gmesh);
    //gengrid.SetBC(gmesh,4,matBCbott);
    //gengrid.SetBC(gmesh,5,matBCright);
    //gengrid.SetBC(gmesh,6,matBCtop);
    //gengrid.SetBC(gmesh,7,matBCleft);
    
    
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
    //  TPZVec<long> pointtopology(1);
    //  pointtopology[0] = 1;
    //
    //  gmesh->CreateGeoElement(EPoint,pointtopology,matPoint,id);
    
    
    //Vetor auxiliar para armazenar as conecções entre elementos:
    
    TPZVec <long> connect(4,0);
    
    
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
    
    {
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
    }
    
    //Impressão da malha geométrica:
    
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;
    
    
    
}

TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}


TPZCompMesh *CMesh_S(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim);//Insere dimensão do modelo
    
    
    //Definição do espaço de aprximação:
    
    cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
    
    
    //Criando material cujo nSTATE = 2:
    
    REAL E = 206.815026; //* @param E elasticity modulus
    REAL nu=0.; //* @param nu poisson coefficient
    REAL fx=0.;//* @param fx forcing function \f$ -x = fx \f$
    REAL fy=0.;//* @param fx forcing function \f$ -x = fx \f$
    int plain=0.; //* @param plainstress = 1 \f$ indicates use of plainstress
    
    TPZMaterial * material = new TPZMixedElasticityMaterial(matID,E,nu,fx,fy,plain,dim);
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
    //TPZMaterial * material2 = new TPZMixedElasticityMaterial(3,E,nu,fx,fy,plain,dim);
    //cmesh->InsertMaterialObject(material2); //Insere material na malha
    
    
    //Dimensões do material (para H1 e descontinuo):
    //TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
    //material->SetMaterial(xkin, xcin, xfin);
    
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZFMatrix<STATE> val2s(2,1,0.);
    val2s(0,0) = 10.0; // vx -> 0
    val2s(1,0) = 0.0; // vy -> 0
    
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
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

TPZCompMesh *CMesh_U(TPZGeoMesh *gmesh, int pOrder)
{
    
    // @omar::
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    
    
    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    //Criando material cujo nSTATE = 2:
    
    REAL E = 206.815026; //* @param E elasticity modulus
    REAL nu=0.; //* @param nu poisson coefficient
    REAL fx=0.;//* @param fx forcing function \f$ -x = fx \f$
    REAL fy=-20.;//* @param fx forcing function \f$ -x = fx \f$
    int plain=0.; //* @param plainstress = 1 \f$ indicates use of plainstress
    
    TPZMaterial * material = new TPZMixedElasticityMaterial(matID,E,nu,fx,fy,plain,dim);
    
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
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    std::set<int> materialids;
    materialids.insert(matID);
    //materialids.insert(3);
    cmesh->AutoBuild(materialids);
    cmesh->LoadReferences();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild();
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    //    cmesh->AdjustBoundaryElements();
    //    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

TPZCompMesh *CMesh_P(TPZGeoMesh *gmesh, int pOrder)
{
    
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
    
    TPZMaterial *material = new TPZMatPoisson3d(matID,dim);//criando material que implementa a formulacao fraca do problema modelo
    
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
    cmesh->AutoBuild(materialids);
    cmesh->LoadReferences();
    //    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    //    cmesh->AutoBuild();
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    long nelem = cmesh->NElements();
    for (long el=0; el<nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(cel);
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



TPZCompMesh *CMesh_Girk(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    int bc_inte_order = 10;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    //    TElasticityExample1 example;
    
    // Criando material:
    
    
    REAL E = 20.59; //* @param E elasticity modulus
    REAL nu=0.; //* @param nu poisson coefficient
    
    
    REAL fx=0.;//* @param fx forcing function \f$ -x = fx \f$
    REAL fy = -32.69;//* @param fx forcing function \f$ -x = fx \f$
    int plain=1.; //* @param plainstress = 1 \f$ indicates use of plainstress
    
    TPZMixedElasticityMaterial * material1 = new TPZMixedElasticityMaterial(1,E,nu,fx,fy,plain,dim);
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
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    REAL x;
    val2(0,0) = 0; // vx -> 0 //val2 represent norm stress;
    val2(1,0) = 0; // vy -> 0
    val1(0,0) = 0;
    val1(1,1) = material1->gBigNumber;
    
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
    
    val2(1,0) = 43.533;
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


TPZCompMesh *CMesh_AxiS(TPZGeoMesh *gmesh, int pOrder, TElasticityExample1 &example)
{
    
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    //    TElasticityExample1 example;
    
    // Criando material:
    
    example.fProblemType = TElasticityExample1::EThiago;
    example.fStressState   = TElasticityExample1::EPlaneStrain;
    REAL E = 206.815026; //* @param E elasticity modulus
    REAL nu=0.; //* @param nu poisson coefStressStatfStressStateficient
    
    
    REAL fx=0.;//* @param fx forcing function \f$ -x = fx \f$
    REAL fy=-20.;//* @param fx forcing function \f$ -x = fx \f$
    int plain=0.; //* @param plainstress = 1 \f$ indicates use of plainstress
    
    TPZMixedElasticityMaterial * material = new TPZMixedElasticityMaterial(matID,E,nu,fx,fy,plain,dim);
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
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    val2(0,0) = 0; // vx -> 0 //val2 represent norm stress;
    val2(1,0) = 0; // vy -> 0
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


TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    //    TElasticityExample1 example;
    
    // Criando material:
    
    // example is initialized in the calling method
    //    example.fProblemType = TElasticityExample1::EThiago;
    //    example.fStressState   = TElasticityExample1::EPlaneStrain;
    
    REAL E = TElasticityExample1::fElast;// 206.815026; //* @param E elasticity modulus
    REAL nu = TElasticityExample1::fNu; // 0.;//0.30400395; //* @param nu poisson coefficient
    
    
    REAL fx=0.;//* @param fx forcing function \f$ -x = fx \f$
    REAL fy=0.;//* @param fx forcing function \f$ -x = fx \f$
    int plain=0; //* @param plainstress = 1 \f$ indicates use of plainstress
    if (TElasticityExample1::fStressState == TElasticityExample1::EPlaneStrain) {
        plain = 0;
    }
    else if(TElasticityExample1::fStressState == TElasticityExample1::EPlaneStress)
    {
        plain = 1;
    }
    else if(TElasticityExample1::fStressState == TElasticityExample1::EAxiSymmetric)
    {
        plain = 0;
    }
    else
    {
        DebugStop();
    }
    
    TPZMixedElasticityMaterial * material = new TPZMixedElasticityMaterial(matID,E,nu,fx,fy,plain,dim);
    
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
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    val2(0,0) = 0; // vx -> 0 //val2 represent norm stress;
    val2(1,0) = 0; // vy -> 0
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


void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    long nel = gmesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->MaterialId() != matfrom) {
            continue;
        }
        
        int nsides= gel->NSides();
        
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        if (celstack.size() != 2) {
            DebugStop();
        }
        gel->SetMaterialId(mattarget);
        long index;
        new TPZMultiphysicsInterfaceElement(cmesh,gel,index,celstack[1],celstack[0]);
    }
    
}


void CreateCondensedElements(TPZCompMesh *cmesh)
{
    long nel = cmesh->NElements();
    cmesh->ComputeNodElCon();
    for(long el = 0; el<nel ; el++)
    {
        TPZCompEl *cel = cmesh->Element(el);
        TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(!mphys) continue;
        TPZGeoEl *gel = mphys->Reference();
        if(gel->Dimension() != cmesh->Dimension()) continue;
        
        long newconnectindex = cmesh->AllocateNewConnect(4,1,1);
        cmesh->ConnectVec()[newconnectindex].SetLagrangeMultiplier(2);
        cmesh->ConnectVec()[newconnectindex].IncrementElConnected();
        cmesh->ConnectVec()[newconnectindex].IncrementElConnected();
        REAL delx = abs(gel->Node(1).Coord(0)-gel->Node(0).Coord(0));
        REAL dely = abs(gel->Node(1).Coord(1)-gel->Node(0).Coord(1));
        int constridf = 0;
        if(delx > dely)
        {
            constridf = 1;
        }
        int numel = mphys->NumberOfCompElementsInsideThisCompEl();
        if(numel != 3) DebugStop();
        TPZManVector<int> numconnects(numel,0);
        for(int i=0; i<numel; i++) numconnects[i] = mphys->Element(i)->NConnects();
        TPZManVector<TPZCompElLagrange::TLagrange, 4> EquationDelay(4);
        EquationDelay[0].fConnect[0] = mphys->ConnectIndex(numconnects[0]);
        EquationDelay[0].fIdf[0] = 0;
        EquationDelay[0].fConnect[1] = newconnectindex;
        EquationDelay[0].fIdf[1] = 0;
        EquationDelay[1].fConnect[0] = mphys->ConnectIndex(numconnects[0]);
        EquationDelay[1].fIdf[0] = 1;
        EquationDelay[1].fConnect[1] = newconnectindex;
        EquationDelay[1].fIdf[1] = 1;
        EquationDelay[2].fConnect[0] = mphys->ConnectIndex(numconnects[0]+1);
        EquationDelay[2].fIdf[0] = constridf;
        EquationDelay[2].fConnect[1] = newconnectindex;
        EquationDelay[2].fIdf[1] = 2;
        EquationDelay[3].fConnect[0] = mphys->ConnectIndex(numconnects[0]+numconnects[1]+numconnects[2]-1);
        EquationDelay[3].fIdf[0] = 0;
        EquationDelay[3].fConnect[1] = newconnectindex;
        EquationDelay[3].fIdf[1] = 3;
        long elindex;
        new TPZCompElLagrange(*cmesh,EquationDelay,elindex);
        long groupindex;
        TPZElementGroup *elgr = new TPZElementGroup(*cmesh,groupindex);
        elgr->AddElement(cel);
        elgr->AddElement(cmesh->Element(elindex));
        TPZCondensedCompEl *condensed = new TPZCondensedCompEl(elgr);
    }
    cmesh->ExpandSolution();
}

void CreateCondensedElements2(TPZCompMesh *cmesh)
{
    long nel = cmesh->NElements();
    cmesh->ComputeNodElCon();
    for(long el = 0; el<nel ; el++)
    {
        TPZCompEl *cel = cmesh->Element(el);
        TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(!mphys) continue;
        TPZGeoEl *gel = mphys->Reference();
        if(gel->Dimension() != cmesh->Dimension()) continue;
        
        int numel = mphys->NumberOfCompElementsInsideThisCompEl();
        if(numel != 3) DebugStop();
        TPZManVector<int> numconnects(numel,0);
        for(int i=0; i<numel; i++) numconnects[i] = mphys->Element(i)->NConnects();
        cel->Connect(numconnects[0]).IncrementElConnected();
        cel->Connect(numconnects[0]+1).IncrementElConnected();
        cel->Connect(numconnects[0]+numconnects[1]+gel->NCornerNodes()-1).IncrementElConnected();
        TPZCondensedCompEl *condensed = new TPZCondensedCompEl(cel);
    }
}

/// change the order of the internal connect to the given order
void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder)
{
    long nel = cmesh->NElements();
    for(long el = 0; el<nel; el++)
    {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() != cmesh->Dimension())
        {
            continue;
        }
        int nc = cel->NConnects();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel) DebugStop();
        
        intel->ForceSideOrder(gel->NSides()-1,pOrder);
    }
    cmesh->ExpandSolution();
}

// Compute the area


REAL AxiArea(TPZGeoMesh * gmesh, std::set<int> matids)
{
    long nel = gmesh->NElements();
    REAL result = 0.;
    for(long el=0; el<nel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        int matid = gel->MaterialId();
        if(matids.find(matid) == matids.end()) continue;
        int nsides = gel->NSides();
        int order = 1;
        TPZIntPoints *intpoints = gel->CreateSideIntegrationRule(nsides-1,order);
        int np = intpoints->NPoints();
        TPZManVector<REAL,3> point(2,0.), x(3,0.);
        REAL weight;
        
        /** @brief Compute a decomposition of the gradient of the mapping function, as a rotation matrix (Jacobian) and orthonormal basis (axes)  **/
        //        void Jacobian(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const;
        
        
        TPZFNMatrix<6,REAL> jac(2,2),jacinv(2,2),axes(2,3);
        REAL detjac;
        //         TPZManVector<REAL,3> co(3);
        //         std::cout << " gel index "  << el << std::endl;
        //         for(int n=0; n<3; n++)
        //         {
        //             gel->NodePtr(n)->GetCoordinates(co);
        //             std::cout << co << " \n ";
        //         }
        REAL elarea = 0.;
        for(int ip = 0; ip < np; ip++)
        {
            //        virtual void Point(int i, TPZVec<REAL> &pos, REAL &w) const = 0;
            intpoints->Point(ip,point,weight);
            gel->Jacobian(point, jac, axes, detjac, jacinv);
            gel->X(point,x);
            elarea += abs(detjac)*weight*x[0];
        }
        //         std::cout << "elarea " << elarea << std::endl;
        result += elarea;
        delete intpoints;
    }
    return result;
}

// integrate r sig.n on the bottom
STATE IntegrateBottom(TPZCompMesh *cmesh, int targetmatid)
{
    long nel = cmesh->NElements();
    STATE integSigy = 0.;
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        int matid = gel->MaterialId();
        if (targetmatid != matid) {
            continue;
        }
        if(gel->Dimension() != 1) DebugStop();
        TPZGeoElSide gelside(gel,2);
        TPZGeoElSide neighbour = gelside.Neighbour();
        int varindex = neighbour.Element()->Reference()->Material()->VariableIndex("SigmaY");
        if(neighbour.Element()->Dimension() != 2 || neighbour.Element()->Reference() == 0) DebugStop();
        TPZTransform<REAL> tr(1);
        gelside.SideTransform3(neighbour, tr);
        tr = neighbour.Element()->SideToSideTransform(neighbour.Side(), neighbour.Element()->NSides()-1).Multiply(tr);
        TPZIntPoints *rule = gel->CreateSideIntegrationRule(2, 7);
        TPZCompEl *neighcel = neighbour.Element()->Reference();
        int np = rule->NPoints();
        for (int ip = 0; ip<np; ip++) {
            TPZManVector<REAL,3> locpoint(1), volpoint(2);
            REAL weight;
            rule->Point(ip, locpoint, weight);
            tr.Apply(locpoint, volpoint);
            TPZManVector<STATE,3> sol(1),x(3);
            neighcel->Solution(volpoint, varindex, sol);
            gel->X(locpoint, x);
            integSigy += sol[0]*weight*x[0];
        }
        
        delete rule;
    }
    return integSigy;
}

