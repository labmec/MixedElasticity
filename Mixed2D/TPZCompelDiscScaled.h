//
//  TPZCompelDiscScaled.h
//  MixedElasticityGirkmann
//
//  Created by Philippe Devloo on 03/05/18.
//

#ifndef TPZCompelDiscScaled_h
#define TPZCompelDiscScaled_h

#include <stdio.h>
#include <iostream>
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzcompel.h"
#include "pzinterpolationspace.h"
#include "TPZCompelDisc.h"
#include "pzgeoel.h"
#include "pzreal.h"
#include "TPZShapeDisc.h"
#include "tpzautopointer.h"
#include "pzquad.h"
#include "pzfunction.h"
#include "pzfmatrix.h"

class TPZCompElDiscScaled : public TPZCompElDisc {

public:
    
//    void SetScale();
//    
//    virtual void ComputeShape(TPZVec<REAL> &intpoint,TPZMaterialData &data);

    virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphidxi);
   
    void ShapeX(TPZVec<REAL> &X, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);

};

#endif /* TPZCompelDiscScaled_h */

