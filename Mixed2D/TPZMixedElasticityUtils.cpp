#include "TPZMixedElasticityUtils.h"
#include "TPZGenGrid2D.h"
#include "TPZExtendGridDimension.h"
#include "TPZGeoMeshTools.h"
#include "pzcheckgeom.h"

TPZGeoMesh* TPZMixedElasticityUtils::CreateGMesh(int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType, int matIdDomain, std::set<int> &matIdBC) {
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
            gengrid.SetElementType(MMeshType::EQuadrilateral);
            gengrid.SetDistortion(0.25);
            break;
    }

    gengrid.Read(gmesh, matIdDomain);

    int cont = 4;
    for (auto matId : matIdBC)
    {
        gengrid.SetBC(gmesh, cont, matId);
        cont++;
    }
    // gengrid.SetBC(gmesh, 5, matBCright);
    // gengrid.SetBC(gmesh, 6, matBCtop);
    // gengrid.SetBC(gmesh, 7, matBCleft);

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

TPZGeoMesh* TPZMixedElasticityUtils::CreateGMesh3D(int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType, int matIdDomain, std::set<int> &matIdBC) {
//Creating geometric mesh, nodes and elements.
//Including nodes and elements in the mesh object:
    TPZGeoMesh *gmesh2D = CreateGMesh(nelx, nely, hx, hy, x0, y0, meshType, matIdDomain, matIdBC);
    TPZExtendGridDimension extend(gmesh2D, hx);
    TPZGeoMesh *gmesh3D = extend.ExtendedMesh(nelx,-5,-6);
//    delete gmesh2D;
    return gmesh3D;
}