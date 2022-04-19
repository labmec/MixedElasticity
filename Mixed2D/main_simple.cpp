#include <fstream>

#include "pzlog.h"
#include <iostream>
#include <string>

#include <cmath>
#include <set>
#include "common_files.h"
#include "TPZAnalyticSolution.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGeoMeshTools.h"
#include "pzcheckgeom.h"
#include "TPZMHMeshControl.h"

void InsertMaterialObjects(TPZCompMesh &cmeshref);



#ifdef PZ_LOG
static TPZLogger logger("testmhm");
#endif

TPZAnalyticSolution *gAnalytic = 0;

int main(int argc, char *argv[]) {
    //    TPZMaterial::gBigNumber = 1.e16;
    
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    EConfigDarcy conf = EPhil;
    int initial_p = 1;
    int final_p = 1;
    int initial_h = 1;
    int final_h = 1;
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
            conf = EConfigDarcy(atoi(argv[1]));
    };
    int n_ref_p = final_p - initial_p + 1;
    int n_ref_h = final_h - initial_h + 1;
    
#ifdef PZ_USING_MKL
    mkl_set_dynamic(0); // disable automatic adjustment of the number of threads
    mkl_set_num_threads(numthreads);
#endif
    
    std::stringstream rootname;
    double hx = 2, hy = 2; //Dimensões em x e y do domínio
    double x0 = -1;
    double y0 = -1;
    
    //Problem data:
    switch (conf) {
        case EPhil:
            if(dim == 2 || dim == 3)
            {
                TLaplaceExample1 *laplace = new TLaplaceExample1;
                laplace->fExact = TLaplaceExample1::EX;
                gAnalytic = laplace;
            }
            else
                DebugStop();
            hx = 1;
            hy = 1;
            x0 = 0;
            y0 = 0;
            
            rootname << ConfigRootnameDarcy[conf] + "_EX";
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
            int nelx = h_level, nely = h_level; //Number of elements in x and y directions
            std::cout << "********* " << "Number of h refinements: " << href << " (" << nelx << "x" << nely << " elements). p order: " << pref + 1 << ". *********" << std::endl;
            unsigned int nx = nelx + 1, ny = nely + 1; //Number of nodes in x and y directions
            unsigned int stressPOrder = pref + 1; //Polynomial order of the approximation
            int stressInternalPOrder = stressPOrder; //k
            int displacementPOrder = elementType == ETriangular ? stressInternalPOrder - 1 : stressInternalPOrder;
            int rotationPOrder = displacementPOrder;
            TPZGeoMesh *gmesh = 0;
            if(dim == 2)
            {
                TPZManVector<REAL,3> minX = {-1.,-1.,0.};
                TPZManVector<REAL,3> maxX = {1.,1.,0.};
                TPZManVector<int,5> matids = {1,-1,-1,-1,-1};
                TPZManVector<int> nDivs = {nelx,nelx};
                MMeshType meshType = MMeshType::EQuadrilateral;
                rootname << "_quad_";
                bool createBoundEls = true;
                gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                        matids, nDivs, meshType, createBoundEls);
//                TPZCheckGeom check(gmesh);
//                check.UniformRefine(href);
            }
            else if(dim == 3)
            {
                TPZManVector<REAL,3> minX = {-1.,-1.,-1.};
                TPZManVector<REAL,3> maxX = {1.,1.,1.};
                TPZManVector<int,5> matids = {1,-1,-1,-1,-1,-1};
                TPZManVector<int> nDivs = {1,1,1};
                MMeshType meshType = MMeshType::EHexahedral;
                rootname << "_hexa_";
                bool createBoundEls = true;
                gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                        matids, nDivs, meshType, createBoundEls);
                TPZCheckGeom check(gmesh);
                check.UniformRefine(href);
//                gmesh = CreateGMesh3D(nelx, nely, hx, hy, x0, y0, elementType); //Creates the geometric mesh
                
            }
#ifdef PZDEBUG
            std::ofstream fileg("MalhaGeo.txt"); //Prints the geometric mesh in txt format
            std::ofstream filegvtk("MalhaGeo.vtk"); //Prints the geometric mesh in vtk format
            gmesh->Print(fileg);
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk, true);
#endif
            TPZAutoPointer<TPZMultiphysicsCompMesh> cmesh_m_HDiv;
            
            TPZVec<int64_t> coarseindices(gmesh->NElements());
            int64_t count = 0;
            int64_t nel = gmesh->NElements();
            for (int64_t el=0; el<nel; el++) {
                TPZGeoEl *gel = gmesh->Element(el);
                if(gel && gel->Dimension() == gmesh->Dimension() && gel->Level()==0)
                {
                    coarseindices[count++] = el;
                }
            }
            coarseindices.resize(count);
            TPZMHMeshControl control(gmesh);
            control.DefinePartitionbyCoarseIndices(coarseindices);
            control.fMaterialIds = {1};
            control.fMaterialBCIds = {-1};
            if(dim == 2)
            {
                control.SetProblemType(TPZMHMeshControl::EScalar);
            } else if(dim == 3)
            {
                control.SetProblemType(TPZMHMeshControl::EScalar);
            }
            InsertMaterialObjects(control.CMesh());
            control.SetInternalPOrder(stressInternalPOrder);
            control.SetSkeletonPOrder(pref);
            bool substruct = true;
            control.BuildComputationalMesh(substruct);
            cmesh_m_HDiv = control.CMesh();
            
            SolveProblem(cmesh_m_HDiv, rootname);
            
            ComputeError(cmesh_m_HDiv, rootname, href, pref, control);

        }
    }
    return 0;
}

#include "DarcyFlow/TPZHybridDarcyFlow.h"

void InsertMaterialObjects(TPZCompMesh &cmeshref)
{
    TPZCompMesh *cmesh = &cmeshref;
    int dim = cmesh->Dimension();

    TPZHybridDarcyFlow *material = new TPZHybridDarcyFlow(1,dim);
    //Condições de contorno:

    TLaplaceExample1 *laplace = dynamic_cast<TLaplaceExample1 *>(gAnalytic);
    if(!laplace) DebugStop();
    
    material->SetForcingFunction(laplace->ForceFunc(), 3);
    
    material->SetExactSol(laplace->ExactSolution(),3);
    
    cmesh->InsertMaterialObject(material);
    TPZFMatrix<REAL> val1(dim, dim, 0.);
    TPZManVector<REAL,3> val2(dim, 0.);


    {
        auto * BCond1 = material->CreateBC(material, -1, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
        BCond1->SetForcingFunctionBC(gAnalytic->ExactSolution());
        TPZMaterial *bc = dynamic_cast<TPZMaterial *>(BCond1);
        //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
        //BCond1->SetForcingFunction(solucao_exact,bc_inte_order);
        cmesh->InsertMaterialObject(bc); //Insere material na malha
    }

}
