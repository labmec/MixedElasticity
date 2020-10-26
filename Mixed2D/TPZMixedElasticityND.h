/**
 * @file
 * @brief Contains the TPZMixedElasticityND class which implements an N-dimensional elastic material.
 */

#ifndef MIXEDELASMATNDHPP
#define MIXEDELASMATNDHPP

#include <iostream>

#include "TPZMaterial.h"
#include "pzdiscgal.h"

/**
 * @ingroup material
 * @brief This class implements an N-dimensional elastic material
 */
class TPZMixedElasticityND : public TPZDiscontinuousGalerkin {

    struct TElasticityAtPoint
    {
        /** @brief Elasticity modulus */
        REAL fE;
        
        /** @brief Poison coefficient */
        REAL fnu;
        
        /** @brief first Lame Parameter */
        REAL flambda;
        
        /** @brief Second Lame Parameter */
        REAL fmu;

        TElasticityAtPoint(REAL E, REAL nu) : fE(E), fnu(nu)
        {
            flambda = (E * nu) / ((1 + nu)*(1 - 2 * nu));
            fmu = E / (2 * (1 + nu));
        }
    };



public:

    enum MVoigt {
        Exx, Eyy, Exy, Eyx, Exz, Ezx, Eyz, Ezy, Ezz
    };

    /** @brief Default constructor */
    TPZMixedElasticityND();
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
    TPZMixedElasticityND(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress = 1, int Dimension = 1);

    /// dimension of the material

    TPZMixedElasticityND(int id);

    /** @brief Copies the data of one TPZMixedElasticityND object to another */
    TPZMixedElasticityND(const TPZMixedElasticityND &copy);


    virtual int VariableIndex(const std::string &name) override;


    virtual int NSolutionVariables(int var) override;

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

    virtual int NEvalErrors() override;
    
    void ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi);

    void ComputeDeformationVector(TPZVec<STATE> &PhiStress, TPZVec<STATE> &APhiStress, TElasticityAtPoint &elastb);

    void ComputeStressVector(TPZVec<STATE> &Deformation, TPZVec<STATE> &Stress,  TElasticityAtPoint &elast);

    void ElasticityModulusTensor(TPZFMatrix<STATE> &MatrixElast, TElasticityAtPoint &elast);

    /** @brief Creates a new material from the current object   ??*/
    virtual TPZMaterial * NewMaterial()  override {
        return new TPZMixedElasticityND(*this);
    }

    /** @brief Default destructor */
    virtual ~TPZMixedElasticityND();

    /**
     * @brief Set parameters of elastic material:
     * @param Eyoung Young's elasticity modulus
     * @param nu Poisson's coefficient
     */
    void SetElasticParameters(REAL Eyoung, REAL nu) {
        this->SetElasticity(Eyoung, nu);
    }

    /** @brief Set elasticity parameters */
    void SetElasticity(REAL E, REAL nu) {
        fE_const = E; // Young modulus
        fnu_const = nu; // poisson coefficient
        flambda_const = (E * nu) / ((1 + nu)*(1 - 2 * nu));
        fmu_const = E / (2 * (1 + nu));

    }

    /// Set a variable elasticity and poisson coefficient
    void SetElasticityFunction(TPZAutoPointer<TPZFunction<STATE> > func)
    {
        fElasticity = func;
    }
    
    /** @brief Set elasticity parameters in the form of Lame's parameters*/
    void SetLameParameters(REAL Lambda, REAL mu) {
        fE_const = (mu * (3.0 * Lambda + 2.0 * mu)) / (Lambda + mu);
        fnu_const = (Lambda) / (2 * (Lambda + mu));

        flambda_const = Lambda;
        fmu_const = mu;
    }

    /// set the material configuration to AxisSymmetric

    void SetAxisSymmetric() {
        if(fDimension == 3) DebugStop();
        fAxisSymmetric = 1;
    }
    /// Set the material configuration to plane strain

    void SetPlaneStrain() {
        if(fDimension != 2) DebugStop();
        fPlaneStress = 0;
        // fPlaneStrain = 1;
    }

    /// Set the material configuration to plane stress

    void SetPlaneStress() {
        if(fDimension != 2) DebugStop();
        fPlaneStress = 1;
        // fPlaneStrain = 0;
    }

    /** @brief Set forcing function */
    void SetBodyForce(REAL fx, REAL fy, REAL fz = 0.) {
        fForce[0] = fx;
        fForce[1] = fy;
        fForce[2] = fz;
    }

    /** @brief Returns the model dimension */
    int Dimension() const  override {
        return fDimension;
    }

    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() const override;

    /** @brief Print the material data*/
    virtual void Print(std::ostream & out = std::cout) override;

    /** @brief Returns the material name*/
    virtual std::string Name()  override {
        return "TPZMixedElasticityND";
    }

    /** @brief Returns the number of components which form the flux function */
    virtual short NumberOfFluxes() {
        DebugStop();
        return 3;
    }

    /** @brief Returns the number of components which form the flux function */
    virtual int NFluxes()  override {
        // DebugStop();
        return 3;
    }

    /** @name Contribute methods */
    /** @{ */

    /** @brief Calculates the element stiffness matrix */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    /** @brief Calculates the element stiffness matrix - simulate compaction as aditional variable */
    virtual void Contribute(TPZVec<TPZMaterialData> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    /** @brief Calculates the element stiffness matrix using 3 spaces - Stress tensor, displacement, and skew-symmetric tensor (for weak symmetry) */
    virtual void Contribute_3spaces(TPZVec<TPZMaterialData> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /** @brief Calculates the element stiffness matrix using 5 spaces - 3 from Contribute_3spaces() + Rigid body motions, and distributed forces */
    virtual void Contribute_5spaces(TPZVec<TPZMaterialData> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);

    /** @brief Calculates the element stiffness matrix - simulate compaction as aditional variable */
    //    virtual void Contribute(TPZVec<TPZMaterialData> &data, REAL weight, TPZFMatrix<STATE> &ef)
    //    {
    //        DebugStop();
    //    }

    /** @brief Calculates the element stiffness matrix */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override {
        TPZDiscontinuousGalerkin::Contribute(data, weight, ef);
    }


    /** @brief Applies the element boundary conditions Mixed */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;


    /** @brief Applies the element boundary conditions */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight,
            TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;

    /** @brief Applies the element boundary conditions */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight,
            TPZFMatrix<STATE> &ef, TPZBndCond &bc) override {
        TPZDiscontinuousGalerkin::ContributeBC(data, weight, ef, bc);
    }

    //virtual void FillDataRequirements(TPZMaterialData &data);
    virtual void FillDataRequirements(TPZMaterialData &data) override;
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;

    virtual void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data) override;

    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)  override {
        PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
    }

    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)  override {
        PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
    }

    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)  override {
        PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
    }

    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &left, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)  override {
        PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
    }

    /** inner product of two tensors. See Gurtin (2003), p. 5. */
    STATE Inner(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T);

    /// Transform a tensor to a Voigt notation
    void ToVoigt(TPZFMatrix<STATE> &S, TPZVec<STATE> &Svoigt);

    /// Transform a Voigt notation to a tensor
    void FromVoigt(TPZVec<STATE> &Svoigt, TPZFMatrix<STATE> &S);

    /** inner product of two vectors. See Gurtin (2003), p. 5. */
    template<class TVar>
    TVar InnerVec(const TPZVec<TVar> &S, const TPZVec<TVar> &T);

    /** trace of the tensor GradU = Div(U)*/
    STATE Tr(TPZFMatrix<REAL> &GradU);

    /** transpose of the tensor GradU */
    STATE Transpose(TPZFMatrix<REAL> &GradU);

    /** Fill the vector of gradient for each phi */
    void FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<STATE> > &GradPhi);

    /// transform a H1 data structure to a vector data structure
    void FillVecShapeIndex(TPZMaterialData &data);

public:

    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;

    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZVec<TPZMaterialData> &data, int var, TPZVec<STATE> &Solout) override;

    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<STATE> &Solout)  {
        TPZDiscontinuousGalerkin::SolutionDisc(data, dataleft, dataright, var, Solout) ;
    }

    /** @brief Computes the value of the flux function to be used by ZZ error estimator */
    virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux) override;

    /** 
     * @brief Computes the error due to the difference between the interpolated flux \n
     * and the flux computed based on the derivative of the solution
     */
    void Errors(TPZVec<REAL> &x, TPZVec<STATE> &u,
            TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
            TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &values) override; //Cedric

    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors) override;

    virtual int ClassId() const override;

    virtual void Read(TPZStream &buf, void *context) override;

    virtual void Write(TPZStream &buf, int withclassid) const override;



protected:
    /** @brief Forcing vector */
    TPZManVector<REAL, 3> fForce = TPZManVector<REAL, 3>(3, 0.);

    /** @brief Plane stress flag */
    int fPlaneStress = 1;
    
    /** @brief Plane strain flag */
    // int fPlaneStrain = 0;
    
    /// dimension of the material
    int fDimension = 2;

    /** Elasticity function */
    TPZAutoPointer<TPZFunction<STATE> > fElasticity;
    
    /** @brief Elasticity modulus */
    REAL fE_const;
    
    /** @brief Poison coefficient */
    REAL fnu_const;
    
    /** @brief first Lame Parameter */
    REAL flambda_const;
    
    /** @brief Second Lame Parameter */
    REAL fmu_const;
    

    // Matrix A
    TPZFMatrix<STATE> fMatrixA;

    /// flag indicates axix-AxisSymmetric 
    bool fAxisSymmetric = false;

};

#endif
