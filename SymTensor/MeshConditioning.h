//
//  MeshConditioning.hpp
//  SymTensor
//
//  Created by Philippe Devloo on 13/05/25.
//

#ifndef MeshConditioning_hpp
#define MeshConditioning_hpp

#include <stdio.h>
#include "TPZRefPatternDataBase.h"
#include "TPZRefPattern.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"

#include "TPZMultiphysicsCompMesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzgeoelbc.h"
#include "TPZApproxCreator.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZMixedSymElasticityND.h"
#include "TPZInterfaceSymTensor.h"

#include "pzelementgroup.h"
#include "pzcondensedcompel.h"

#include "pzfstrmatrix.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"  //symmetric sparse matrix storage
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZAnalyticSolution.h"

extern TElasticity2DAnalytic gAnalytic;


enum matids {matID = 1, matBCD = 2, matBCN = 3, matHEXAGON = 4, matWrap = 5, matIntFacePos = 6, matIntFaceNeg = 7, matLagrDisp = 8, matIntFacePosLoworder = 9, matIntFaceNegLoworder = 10, matLagrDispLoworder = 11, matTensionLoworder = 12};

/// uniformly refine the mesh
void UniformRefine(TPZGeoMesh *gmesh, int nref);

/// divide the volumetric elements with a specific refinement pattern
void CreateJohnsonMercier(TPZGeoMesh *gmesh);

/// Add geometric boundary elements
void AddWrapElements(TPZGeoMesh *gmesh);

/// Add geometric boundary elements to represent triple hybridization
void AddWrapElementsTripleHybrid(TPZGeoMesh *gmesh);

/// Add matHexagon elements indicating the contour of the Johnson Mercier elements
void AddMatHexagon(TPZGeoMesh *gmesh);

/// create a discontinuous mesh with continuity of center nodes
TPZCompMesh *CreateTensorSpace(TPZGeoMesh *gmesh, int porder, int porderlow = -1);

/// create the displacement space and lagrange multipliers
TPZCompMesh *CreateDisplacementSpace(TPZGeoMesh *gmesh, int porder, int porderlow = -1);

TPZMultiphysicsCompMesh *GenerateMultiphysicsMesh(TPZGeoMesh *gmesh, int porder, int porderlow = -1);

void AddInterfaceElements(TPZMultiphysicsCompMesh *mfmesh, int porderlow = -1);

std::set<int64_t> InactiveInterfaceConnects(TPZMultiphysicsCompMesh *mfmesh);

void GroupElements(TPZMultiphysicsCompMesh *mfmesh, int porderlow = -1);

void AddBCHoneycomb(TPZGeoMesh *gmesh, int bcidD, int bcidN);


#endif /* MeshConditioning_hpp */
