/**
 * @file
 * @brief Contains the TPZMixedElasticityMaterialLocal class which implements a two dimensional elastic material in plane stress or strain.
 */

#ifndef MIXEDELASMATHPP
#define MIXEDELASMATHPP

#include <iostream>

#include "TPZMaterial.h"
#include "TPZMaterialDataT.h"
#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"


/**
 * @ingroup material
 * @brief This class implements a two dimensional elastic material in plane stress or strain
 */
class TPZMixedElasticityMaterialLocal : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
TPZMatErrorCombinedSpaces<STATE>> {

// type alias to improve constructor readability
using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
    TPZMatErrorCombinedSpaces<STATE> >;


public:

    enum MVoigt {
        Exx, Eyy, Exy, Eyx
    };

    /** @brief Default constructor */
    TPZMixedElasticityMaterialLocal();
    /** 
     * @brief Creates an elastic material with:
     * @param id material id
     * @param E elasticity modulus
     * @param nu poisson coefficient
     * @param fx forcing function \f$ -x = fx \f$ 
     * @param fy forcing function \f$ -y = fy \f$
     * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
     */
    TPZMixedElasticityMaterialLocal(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress = 1, int fDimension = 1);

    /// dimension of the material

    TPZMixedElasticityMaterialLocal(int id);

    /** @brief Copies the data of one TPZMixedElasticityMaterialLocal object to another */
    TPZMixedElasticityMaterialLocal(const TPZMixedElasticityMaterialLocal &copy);


    int VariableIndex(const std::string &name);


    virtual int NSolutionVariables(int var) const override;

    /** index of Stress */
    int SIndex() {
        return 0;
    }

    /** index of Displacement */
    int UIndex() {
        return 1;
    }

    /** index of Rotation */
    int PIndex() {
        return 2;
    }

    virtual int NEvalErrors() const override;
    
    void ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi);

    void ComputeDeformationVector(TPZVec<STATE> &PhiStress, TPZVec<STATE> &APhiStress);

    void ComputeStressVector(TPZVec<STATE> &Deformation, TPZVec<STATE> &Stress);

    void ElasticityModulusTensor(TPZFMatrix<STATE> &MatrixElast);

    /** @brief Creates a new material from the current object   ??*/
    virtual TPZMaterial * NewMaterial() {
        return new TPZMixedElasticityMaterialLocal(*this);
    }

    /** @brief Default destructor */
    virtual ~TPZMixedElasticityMaterialLocal();

    /**
     * @brief Set parameters of elastic material:
     * @param First  Lame Parameter Lambda
     * @param Second Lame Parameter Mu -> G
     * @param fx forcing function \f$ -x = 0 \f$
     * @param fy forcing function \f$ -y = 0 \f$
     */
    void SetElasticParameters(REAL Eyoung, REAL nu) {
        this->SetElasticity(Eyoung, nu);
    }

    /** @brief Set elasticity parameters */
    void SetElasticity(REAL E, REAL nu) {
        fE = E; // Young modulus
        fnu = nu; // poisson coefficient
        flambda = (E * nu) / ((1 + nu)*(1 - 2 * nu));
        fmu = E / (2 * (1 + nu));

    }

    /** @brief Set elasticity parameters */
    void SetLameParameters(REAL Lambda, REAL mu) {
        fE = (mu * (3.0 * Lambda + 2.0 * mu)) / (Lambda + mu);
        fnu = (Lambda) / (2 * (Lambda + mu));

        flambda = Lambda;
        fmu = mu;
    }

    /// set the material configuration to AxisSymmetric

    void SetAxisSymmetric() {
        fAxisSymmetric = 1;
    }
    /// Set the material configuration to plane strain

    void SetPlaneStrain() {
        fPlaneStress = 0;
    }

    /// Set the material configuration to plane stress

    void SetPlaneStress() {
        fPlaneStress = 1;
    }

    /** @brief Set forcing function */
    void SetBodyForce(REAL fx, REAL fy) {
        fForce[0] = fx;
        fForce[1] = fy;
        fForce[2] = 0.;
    }

    /** @brief Returns the model dimension */
    int Dimension() const override{
        return 2;
    }

    /** @brief Returns the number of state variables associated with the material */
	virtual int NStateVariables() const override{ return 2; }

    /** @brief Print the material data*/
    virtual void Print(std::ostream & out = std::cout);

    /** @brief Returns the material name*/
    std::string Name() {
        return "TPZMixedElasticityMaterialLocal";
    }

    /** @brief Returns the number of components which form the flux function */
    virtual short NumberOfFluxes() {
        return 3;
    }

    /** @brief Returns the number of components which form the flux function */
    virtual int NFluxes() {
        return 3;
    }

    /** @name Contribute */
    /** @{ */
    /**
     * @brief It computes a contribution to the stiffness matrix
     * and load vector at one integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                            REAL weight,TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;
    /** @name ContributeBC
        @ingroup Contribute*/
    /**@{*/
    /**
     * @brief It computes a contribution to the stiffness matrix
     * and load vector at one BC integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                              REAL weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef,
                              TPZBndCondT<STATE> &bc) override;
    /**@}*/
    /**@}*/

    /**
     * @brief Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered
     * as not necessary.
     */
    virtual void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;


    STATE GetLambda() const {
        STATE lambda = (fnu * fE) / ((1. + fnu)*(1. - 2. * fnu));
        return lambda;
    }

    STATE GetMU() const {
        STATE mu = fE / (2. * (1. + fnu));
        return mu;
    }


    /** inner product of two tensors. See Gurtin (2003), p. 5. */
    STATE Inner(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T);

    /// Transform a tensor to a Voigt notation
    static void ToVoigt(TPZFMatrix<STATE> &S, TPZVec<STATE> &Svoigt);

    /// Transform a Voigt notation to a tensor
    static void FromVoigt(TPZVec<STATE> &Svoigt, TPZFMatrix<STATE> &S);

    /** inner product of two vectors. See Gurtin (2003), p. 5. */
    template<class TVar>
    TVar InnerVec(const TPZVec<TVar> &S, const TPZVec<TVar> &T);

    /** trace of the tensor GradU = Div(U)*/
    STATE Tr(TPZFMatrix<REAL> &GradU);

    /** transpose of the tensor GradU = Div(U)*/
    STATE Transpose(TPZFMatrix<REAL> &GradU);

    /** Fill the vector of gradient for each phi */
    void FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<STATE> > &GradPhi);

    /// transform a H1 data structure to a vector data structure
    void FillVecShapeIndex(TPZMaterialData &data);

public:
    /*
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


    /** @brief Returns the elasticity modulus E */
    REAL E() {
        return fE;
    }

    /** @brief Returns the poison coefficient modulus E */
    REAL Nu() {
        return fnu;
    }

    virtual int ClassId() const override;

    virtual void Read(TPZStream &buf, void *context) override;

    virtual void Write(TPZStream &buf, int withclassid) const override;



protected:
    /** @brief Elasticity modulus */
    REAL fE;

    /** @brief Poison coefficient */
    REAL fnu;

    /** @brief first Lame Parameter */
    REAL flambda;

    /** @brief Second Lame Parameter */
    REAL fmu;

    /** @brief Forcing vector */
    TPZManVector<REAL, 3> fForce = TPZManVector<REAL, 3>(2, 0.);

    /** @brief Uses plain stress */
    int fPlaneStress = 1;

    /// dimension of the material
    int fDimension = 2;

    // Matrix A
    TPZFMatrix<STATE> fMatrixA;

    /// flag indicates axix-AxisSymmetric 
    bool fAxisSymmetric = false;

};

#endif
