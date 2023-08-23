//
//  VoronoiAnalise.hpp
//  mhm-elas-3d
//
//  Created by Philippe Devloo on 21/08/23.
//

#ifndef VoronoiAnalise_hpp
#define VoronoiAnalise_hpp

#include <stdio.h>
#include <utility>
#include "pzmanvector.h"
#include "pzfmatrix.h"


struct TPZFaceDefinition {

    // indices of the left and right domain
    std::pair<int64_t,int64_t> fLeftRightDomain;
    // normal between the domains
    TPZManVector<REAL,3> fNormal;
    // center of the surface linking the domains
    TPZManVector<REAL,3> fCenter;
    // two direct to define the tangent tractions
    TPZFNMatrix<6,REAL> fAxes;
    // scale factor
    REAL fScale;
    // order of approximation
    int fPOrder;
    // connect index associated with the face
    int64_t fConnectIndex;
    // list of aligned face elements between both domains
    std::list<TPZCompEl *> fElementList;
    
    // Total area of the interface
    REAL fTotalArea = 0.;
    
    // return the number of shape functions
    int NShape();
    // compute the value of the shape functions at a point in the cartesian space
    void Shape(const TPZVec<REAL> &x, TPZFMatrix<REAL> &phi);
    
    // project and restrain the shape functions
    void Project(TPZCompEl *cel);
    
    // Initialize geometric datastructure
    void InitializeDataStructure(const int pord);
    
    /// Compute the center of of the polygon and its area
    void ComputeCenter();
    
    
};
#endif /* VoronoiAnalise_hpp */
