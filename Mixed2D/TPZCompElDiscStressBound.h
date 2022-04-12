//
//  TPZCompElDiscStressBound.h
//
//  Created by Jeferson Fernandes on 12/04/22.
//

#ifndef TPZCompElDiscStressBoundBound_h
#define TPZCompElDiscStressBoundBound_h

#include <stdio.h>
#include <iostream>
#include "TPZCompElDisc.h"
#include "pzgeoel.h"
#include "pzinterpolationspace.h"
#include "TPZMaterialData.h"


template<class StressDef>
class TPZCompElDiscStressBoundBound : public TPZCompElDisc {

public:

    /** @brief Default constructor */
    TPZCompElDiscStressBound() : TPZCompElDisc() {
    }
    
    /** @brief Constructor of the discontinuous element associated with geometric element */
    TPZCompElDiscStressBound(TPZCompMesh &mesh,TPZGeoEl *ref);

    /** @brief Copy constructor */
    TPZCompElDiscStressBound(TPZCompMesh &mesh, const TPZCompElDiscStressBound &copy);

    virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override {
		return new TPZCompElDiscStressBound(mesh,*this);
	}

    /**
     * @brief Initialize a material data and its attributes based on element dimension, number
     * of state variables and material definitions
     */
    virtual void InitMaterialData(TPZMaterialData &data) override;

    /** 
	 * @brief Compute shape functions based on master element in the classical FEM manne. 
	 * @param[in] intpoint point in master element coordinates 
	 * @param[in] data stores all input data
	 */
    virtual void ComputeShape(TPZVec<REAL> &intpoint,TPZMaterialData &data) override;


    int ClassId() const override;


    /** @brief Returns the shapes number of the element */
	virtual int NShapeF() const override;

    void SetStressDef(TPZAutoPointer<StressDef > stress){
        fStress = stress;
    };

protected:
	/** @brief A pz function to allow the inclusion of extra shape functions which are defined externally. 
     * The arguments in this case are: 
    */
    TPZAutoPointer<StressDef> fStress;


};

#endif