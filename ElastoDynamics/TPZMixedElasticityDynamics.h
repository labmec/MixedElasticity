#pragma once

#include <iostream>

#include "TPZMaterial.h"
#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "Material/Elasticity/TPZMixedElasticityND.h"

class TPZMixedElasticityDynamics : public TPZMixedElasticityND {

public:
    enum EContributeType
    {
        EStatic = 0,
        EDynamic = 1,
        ERhs = 2,
        EAll = 3
    };

    TPZMixedElasticityDynamics();
    
    TPZMixedElasticityDynamics(int id, REAL E, REAL nu, REAL rho, REAL deltat, int dimension = 2, int planestress = 1);

    TPZMixedElasticityDynamics(int id, int dimension = 2);

    TPZMixedElasticityDynamics(const TPZMixedElasticityDynamics &copy);

    virtual ~TPZMixedElasticityDynamics();

    virtual int VariableIndex(const std::string &name) const override;

    virtual int NSolutionVariables(int var) const override;

    virtual int NEvalErrors() const override;

    virtual int NStateVariables() const override;
    
    /** @name Contribute methods */
    /** @{ */

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

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

    /** @brief Calculates the element stiffness matrix using 3 spaces - Stress tensor, displacement, and skew-symmetric tensor (for weak symmetry) */
    virtual void Contribute_3spaces(const TPZVec<TPZMaterialDataT<STATE>> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    /** @brief Calculates the element stiffness matrix using 3 spaces - Stress tensor, displacement, and skew-symmetric tensor (for weak symmetry) */
    virtual void Contribute_5spaces(const TPZVec<TPZMaterialDataT<STATE>> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    virtual void ContributeRigidBodyMode(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, int firstRBEq, int firstRBspace);

    /*
     * @brief Fill requirements for volumetric contribute
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

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

protected:

    REAL fdt; // time step for dynamics
    REAL fRho; // density of the material
    EContributeType fContributeType = EAll; // type of contribution to be assembled (static or dynamic) see documentation
};
