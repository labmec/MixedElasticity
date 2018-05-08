//
//  TPZCompElDiscScaled.cpp
//  MixedElasticityGirkmann
//
//  Created by Philippe Devloo on 03/05/18.
//

#include "TPZCompElDiscScaled.h"
#include <math.h>


void TPZCompElDiscScaled::ComputeScale() {
    TPZGeoEl * ref = this->Reference();
    fScale = ref->CharacteristicSize();
}

void TPZCompElDiscScaled::Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
    phi*=(1/fScale);
    dphi*(1/fScale);
}




