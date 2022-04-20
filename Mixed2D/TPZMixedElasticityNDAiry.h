/**
 * @file
 * @brief Contains the TPZMixedElasticityNDAiry class which implements an N-dimensional elastic material enriched
 * with analytical Airy stress function.
 */

#ifndef MIXEDELASMATNDHPPAIRY
#define MIXEDELASMATNDHPPAIRY

#include <iostream>

#include "TPZMixedElasticityND.h"
#include "TPZMaterial.h"
#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"

/**
 * @ingroup material
 * @brief This class implements an N-dimensional elastic material
 */
class TPZMixedElasticityNDAiry : public TPZMixedElasticityND  {

// type alias to improve constructor readability
using TBase = TPZMixedElasticityND;

public:
    /** @brief Default constructor */
    TPZMixedElasticityNDAiry();
    /** 
     * @brief Creates an elastic material with:
     * @param id material id
     * @param E elasticity modulus
     * @param nu poisson coefficient
     * @param fx forcing function \f$ -x = fx \f$ 
     * @param fy forcing function \f$ -y = fy \f$
     * @param planestress \f$ planestress = 1 \f$ indicates use of planestress
     * @param Dimension of the problem
     */
    TPZMixedElasticityNDAiry(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress = 1, int Dimension = 1);

    /** @brief Copies the data of one TPZMixedElasticityND object to another */
    TPZMixedElasticityNDAiry(const TPZMixedElasticityNDAiry &copy);

    /** @brief Default destructor */
    virtual ~TPZMixedElasticityNDAiry();

    /** @brief Creates a new material from the current object   ??*/
    virtual TPZMaterial * NewMaterial()  const override {
        return new TPZMixedElasticityNDAiry(*this);
    }

    /** @brief Returns the material name*/
    virtual std::string Name() const override {
        return "TPZMixedElasticityNDAiry";
    }

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /** @brief Calculates the element stiffness matrix using 4 spaces - Stress tensor, Airy stress functions, displacement, and skew-symmetric tensor (for weak symmetry) */
    virtual void Contribute_4spaces(const TPZVec<TPZMaterialDataT<STATE>> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
                      
    /**
     * @brief Returns the solution associated with the var index based on the
     * finite element approximation at a point
     * @param [in] datavec material data associated with a given integration point
     * @param [in] var index of the variable to be calculated
     * @param [out] solOut vector to store the solution
     */
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut) override;


    /**
     * @brief Calculates the approximation error at a point
     * @param [in] data material data of the integration point
     * @param [out] errors calculated errors
     */
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;

    virtual int ClassId() const override;

    virtual int VariableIndex(const std::string &name) const override;

};

#endif
