//
//  TPZCompelDiscScaled.cpp
//  MixedElasticityGirkmann
//
//  Created by Philippe Devloo on 03/05/18.
//

#include "TPZCompelDiscScaled.h"
#include "pzinterpolationspace.h"
#include "TPZCompElDisc.h"
#include "pzgeoel.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include <math.h>
#include "pzmaterialdata.h"
#include "tpzautopointer.h"


TPZCompElDiscScaled::~TPZCompElDiscScaled() {
    TPZGeoEl * ref = this->Reference();
    if (ref){
        if(ref->Reference() == this) ref->ResetReference();
    }//if (ref)
}

void TPZCompElDiscScaled::SetScale() {
    TPZGeoEl * ref = this->Reference();
    fScale = ref->CharacteristicSize();
}


void TPZCompElDiscScaled::Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
    phi=phi*(1/fScale);
    dphi=dphi*(1/fScale);
}




