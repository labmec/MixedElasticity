//
//  TPZCompelDiscScaled.cpp
//  MixedElasticityGirkmann
//
//  Created by Philippe Devloo on 03/05/18.
//

#include "TPZCompelDiscScaled.h"
#include "pzinterpolationspace.h"
#include "TPZCompelDisc.h"
#include "pztransfer.h"
#include "pzelmat.h"
#include "pzmatrix.h"
#include "pzelmat.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzerror.h"
#include "pzconnect.h"
#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzmanvector.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
#include "pzgraphel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "pzgraphel1d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pztrigraph.h"
#include "tpzgraphelt2dmapped.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel.h"

#include <sstream>
#include <cmath>
#include <algorithm>

#include "time.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include <math.h>
#include <stdio.h>
#include "pzmaterialdata.h"
#include "tpzautopointer.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcompeldisc"));
#endif

//void TPZCompElDiscScaled::ComputeShape(TPZVec<REAL> &intpoint,TPZMaterialData &data){
//    this->ComputeShape(intpoint, data.x, data.jacobian, data.axes,data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
//}

void TPZCompElDiscScaled::Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
    
    if(fUseQsiEta==true){
        this->ShapeX(qsi,phi,dphi);
    }else{
        
        TPZManVector<REAL,4> x(3);
        this->Reference()->X(qsi,x);
        this->ShapeX(x,phi,dphi);
    }
}


void TPZCompElDiscScaled::ShapeX(TPZVec<REAL> &X, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi){
    const int Degree = this->Degree();
    if(Degree < 0) return;
    const int dim = this->Dimension();
    
    TPZGeoEl *ref = this->Reference();
    REAL scale = ref->CharacteristicSize();
    this->SetConstC(scale);
    
    pzshape::TPZShapeDisc::Shape(dim, fConstC,fCenterPoint,X,Degree,phi,dphi,fShapefunctionType);
    //now appending external shape functions
    this->AppendExternalShapeFunctions(X,phi,dphi);
}//method



