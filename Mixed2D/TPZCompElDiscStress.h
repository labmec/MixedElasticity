//
//  TPZCompElDiscStress.h
//
//  Created by Jeferson Fernandes on 08/04/22.
//

#ifndef TPZCompElDiscStress_h
#define TPZCompElDiscStress_h

#include <stdio.h>
#include <iostream>
#include "TPZCompElDisc.h"
#include "pzgeoel.h"
#include "pzinterpolationspace.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"


template<class StressDef>
class TPZCompElDiscStress : public TPZCompElDisc {

public:

    /** @brief Default constructor */
    TPZCompElDiscStress() : TPZCompElDisc() {
    }
    
    /** @brief Constructor of the discontinuous element associated with geometric element */
    TPZCompElDiscStress(TPZCompMesh &mesh,TPZGeoEl *ref);

    /** @brief Copy constructor */
    TPZCompElDiscStress(TPZCompMesh &mesh, const TPZCompElDiscStress &copy);

    virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override {
		return new TPZCompElDiscStress(mesh,*this);
	}

    /** @brief Default destructor */
	virtual ~TPZCompElDiscStress(){

    };

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

    void SetStressDef(StressDef* stress){
        fStress = stress;
    };
    
    virtual void ReallyComputeSolution(TPZMaterialDataT<STATE> &data) override;
    
    /** @brief Returns the number of shapefunctions associated with a connect*/
	virtual int NConnectShapeF(int inod, int order) const override;

    /** @brief Returns the max order of interpolation. */
	virtual int MaxOrder() override;
    
protected:
	/** @brief A pz function to allow the inclusion of extra shape functions which are defined externally. 
     * The arguments in this case are: 
    */
    StressDef* fStress;


};

#endif