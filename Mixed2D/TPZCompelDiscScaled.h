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
#include "TPZCompElDisc.h"
#include "pzgeoel.h"
#include "TPZShapeDisc.h"


class TPZCompElDiscScaled : public TPZCompElDisc {

protected:
    
    REAL fScale = 1.;
    
public:
    
    /** @brief Default constructor */

    TPZCompElDiscScaled() : TPZCompElDisc(){ }
    
    /** @brief Constructor of the discontinuous element associated with geometric element */
    TPZCompElDiscScaled(TPZCompMesh &mesh,TPZGeoEl *ref,int64_t &index) : TPZCompElDisc(mesh,ref,index){ } ;//original

    /** @brief Copy constructor */
    TPZCompElDiscScaled(TPZCompMesh &mesh, const TPZCompElDiscScaled &copy) : TPZCompElDisc(mesh,copy), fScale(copy.fScale){ };

    /** @brief Default destructor */
    ~TPZCompElDiscScaled();
    
    /** @brief Return the shape functions scale */
    REAL Scale() const {return fScale;}

    /** @brief Set a constant value to shape functions scale */
    void SetScale(REAL c){fScale = c;}

    /** @brief Set the shape functions scale related to reference characteristic size*/
    void SetScale();
    
    /** @brief Compute shape functions divided by scale value*/
    virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphidxi);
   


};

#endif /* TPZCompelDiscScaled_h */

