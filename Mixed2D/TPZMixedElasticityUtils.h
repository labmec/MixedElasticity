#ifndef MIXED_ELASTICITY_UTILS
#define MIXED_ELASTICITY_UTILS

#include <iostream>
#include <set>
#include "TPZLinearAnalysis.h"
#include <pzgmesh.h> 

class TPZMixedElasticityUtils{

public:

    TPZMixedElasticityUtils() = default;

    // local enum for mesh types @ToDo these names might lead to confusion. We should consider changing.
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
    TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType, int matIdDomain, std::set<int> &matIdBC);

    /**
    * @brief Funcao para criar a malha geometrica do problema a ser simulado
    * @note A malha sera tridimensional formada por nel elementos de tamanho elsize
    * @param nelx, nely number of elements in the x and y direction (the number of elements in the z direction = nelx
    * @param hx, hy size of the domain in x and y (size in z = hx)
    * @param x0, y0 bottom left point coordinate (bottom left in z = 0)
    * @param meshType = triangle, quadrilateral or trapeze
    */
    TPZGeoMesh *CreateGMesh3D(int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType, int matIdDomain, std::set<int> &matIdBC);



};

#endif