#include "TPZMixedElasticityDynamics.h"
#include "pzelmat.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "TPZMaterialDataT.h"
#include <math.h>

TPZMixedElasticityDynamics::TPZMixedElasticityDynamics() {}

TPZMixedElasticityDynamics::TPZMixedElasticityDynamics(int id, REAL E, REAL nu, REAL rho, REAL deltat, int dimension, int planestress) : TPZMixedElasticityND(id, E, nu, 0., 0., planestress, dimension), fRho(rho), fdt(deltat) {
    this->SetElasticity(E, nu);
    fForce[0] = 0.; // X component of the body force
    fForce[1] = 0.; // Y component of the body force
    fForce[2] = 0.; // Z component of the body force
}

TPZMixedElasticityDynamics::TPZMixedElasticityDynamics(int id, int dimension) : TPZMixedElasticityND(id, dimension), fRho(0.), fdt(0.) {
    this->SetElasticity(1., 0.);
    fForce[0] = 0.; // X component of the body force
    fForce[1] = 0.; // Y component of the body force
    fForce[2] = 0.; // Z component of the body force
}

TPZMixedElasticityDynamics::~TPZMixedElasticityDynamics() {}

TPZMixedElasticityDynamics::TPZMixedElasticityDynamics(const TPZMixedElasticityDynamics &copy) : 
TPZMixedElasticityND(copy),
fRho(copy.fRho),
fdt(copy.fdt) {}

int TPZMixedElasticityDynamics::VariableIndex(const std::string &name) const {

    if (!strcmp("Displacement", name.c_str())) return 1;
    if (!strcmp("Stress", name.c_str())) return 2;
    if (!strcmp("Flux", name.c_str())) return 2;
    if (!strcmp("Rotation", name.c_str())) return 3;
    if (!strcmp("Strain", name.c_str())) return 4;
    if (!strcmp("ExactDisplacement", name.c_str())) return 5;
    if (!strcmp("ExactStress", name.c_str())) return 6;
    
    std::cout << "\n\nVar index not implemented\n\n";
    DebugStop();
    
    return 0;
}

int TPZMixedElasticityDynamics::NSolutionVariables(int var) const {
    int nstate = fDimension;
    switch (var) {
        case 1:
            return nstate;
        case 2:
        case 3:
        case 4:
            return nstate*nstate;
        case 5:
            return nstate;
        case 6:
            return nstate*nstate;
        default:
            std::cout << "\n\nVar index not implemented!!!\n\n";
            DebugStop();
            return 0;
    }
}

int TPZMixedElasticityDynamics::ClassId() const {
    return Hash("TPZMixedElasticityDynamics") ^ TPZMixedElasticityND::ClassId() << 1;
}

int TPZMixedElasticityDynamics::NEvalErrors() const {
    return 7;
}

void TPZMixedElasticityDynamics::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {
    //values[0] = 0.;
    //TPZManVector<REAL, 4> SigmaV(4, 0.), sigma_exactV(4, 0.), eps_exactV(4, 0.), EPSZV(4, 0.);
    TPZVec<STATE> u_exact(fDimension,0.);
    TPZFMatrix<STATE> du_exact(fDimension,fDimension,0.);
    TPZManVector<STATE> divsigma(fDimension,0.);
    if (this->fExactSol) {
        this->fExactSol(data[0].x, u_exact, du_exact);
    }
    if (this->fForcingFunction) {
        this->fForcingFunction(data[0].x, divsigma);
    }

    int nstate = fDimension;
    int matdim = nstate*nstate;
	TPZManVector<STATE, 9> SigmaV(matdim, 0.), sigma_exactV(matdim, 0.), eps_exactV(matdim, 0.), EPSZV(matdim, 0.);
	TPZFNMatrix<9, STATE> sigma(nstate, nstate, 0.), eps(nstate, nstate, 0.), grad(nstate, nstate, 0.);
    TPZFNMatrix<9, STATE> eps_exact(nstate, nstate, 0.);
    TPZManVector<REAL, 3> x = data[0].x;
    int dim = Dimension();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            sigma(i, j) = data[0].sol[0][j + i * 3];
            if (fAxisSymmetric) {
                sigma(i, j) /= x[0];
            }
        }
    }
    ToVoigt(sigma, SigmaV);
    TPZManVector<STATE, 3> divSigma(dim,0.);
    for (int i = 0; i < dim; i++) {
        divSigma[i] = data[0].divsol[0][i];
    }

    TPZManVector<STATE, 3> disp(dim);
    for (int i = 0; i < dim; i++) {
        disp[i] = data[1].sol[0][i];
    }

#ifdef LOG4CXX
    if (logdata->isDebugEnabled())
 {
        std::stringstream sout;
        sout << "DISP*************************** = " << disp << std::endl;
        sigma.Print("sigma************************************ = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(logdata, sout.str())
    }
#endif


    eps_exact(0, 0) = du_exact(0, 0);
    eps_exact(1, 0) = eps_exact(0, 1) = 0.5 * (du_exact(0, 1) + du_exact(1, 0));
    eps_exact(1, 1) = du_exact(1, 1);
    if(dim == 3)
    {
        eps_exact(2, 2) = du_exact(2, 2);
        eps_exact(2, 0) = eps_exact(0, 2) = 0.5 * (du_exact(0, 2) + du_exact(2, 0));
        eps_exact(1, 2) = eps_exact(2, 1) = 0.5 * (du_exact(2, 1) + du_exact(1, 2));
    }
    ToVoigt(eps_exact, eps_exactV);
    
    TElasticityAtPoint elast(fE_const,fnu_const);
    if(fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
		TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity(x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }
    
    ComputeStressVector(eps_exactV, sigma_exactV, elast);

    TPZManVector<STATE,3> rotation(1,0.), rotationExact(1,0.);
    if(dim == 3) {
        rotation.Resize(fDimension,0.);
        rotationExact.Resize(fDimension,0.);
    }
    rotation[0] = data[2].sol[0][0];
    rotationExact[0] = 0.5*(du_exact(1, 0)-du_exact(0, 1));
    if(dim == 3)
    {
        rotation[1] = data[2].sol[0][1];
        rotationExact[1] = 0.5*(du_exact(2, 0)-du_exact(0, 2));
        rotation[2] = data[2].sol[0][2];
        rotationExact[2] = 0.5*(du_exact(2, 1)-du_exact(1, 2));

    }
    // std::cout << "Rotation = "<< rotation << std::endl;
    // std::cout << "rotationExact = "<< rotationExact << std::endl << std::endl;

#ifdef PZDEBUG
    {
        TPZManVector<STATE, 9> eps_again(matdim);
        ComputeDeformationVector(sigma_exactV, eps_again, elast);
        for (int i = 0; i < matdim; i++) {
            if (abs(eps_again[i] - eps_exactV[i]) > 1.e-8) {
                DebugStop();
            }
        }
    }
#endif
    
    //L2_Error: for the displacement (u)
    errors[3] = 0.;
    for (int idf = 0; idf < dim; idf++) {
        errors[3] += (disp[idf] - u_exact[idf])*(disp[idf] - u_exact[idf]);
    }
    //Energe norm
    //TPZManVector<REAL,4> SIGMA(4,0.) , EPSZ(4,0.);

    
    ToVoigt(sigma, SigmaV);

    ComputeDeformationVector(SigmaV, EPSZV, elast);

    errors[0] = 0.;
    errors[1] = 0.;
    for (int i = 0; i < matdim; i++) {
        //L2_Error: for the stress tensor (sigma)
        errors[0] += (SigmaV[i] - sigma_exactV[i])*(SigmaV[i] - sigma_exactV[i]);
        
        //Energy_Error: for the stress tensor (sigma)
        errors[1] += (SigmaV[i] - sigma_exactV[i])*(EPSZV[i] - eps_exactV[i]);
    }
    if (errors[1] < 0.) {
        std::cout << "I should stop \n";
    }

    TPZManVector<STATE,3> divSigmaExact(fDimension,0.);
    if(HasForcingFunction())
    {
        fForcingFunction(data[0].x, divSigmaExact);
        for(int i=0; i<fDimension; i++) divSigmaExact[i] *= -1.;
    }
    errors[2] = 0.;
    errors[4] = 0.;
    for (int idf = 0; idf < fDimension; idf++)
    {
        //L2_Error: divergent of the stress tensor (div(sigma))
        
        errors[2] += (divSigma[idf]-divSigmaExact[idf])*(divSigma[idf]-divSigmaExact[idf]);
        
    }
    int nrot = 1;
    if(fDimension == 3) nrot = 3;
    for(int ir = 0; ir<nrot; ir++)
    {
        //L2_Error: for the rotation (q)
        errors[4] += (rotation[ir]-rotationExact[ir])*(rotation[ir]-rotationExact[ir]);
    }
    //L2_Error: for the symmetry measure (asym)
    errors[5] = pow(SigmaV[Exy]-SigmaV[Eyx],2);
    if(fDimension == 3)
    {
        errors[5] += pow(SigmaV[Exz]-SigmaV[Ezx],2);
        errors[5] += pow(SigmaV[Ezy]-SigmaV[Eyz],2);

    }
    
    //Energy_Norm: for the exact displacement solution
    errors[6] = 0.;
    for(int i=0; i<matdim; i++)
    {
        errors[6] += eps_exactV[i]*sigma_exactV[i];
    }
}

void TPZMixedElasticityDynamics::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    int nref = datavec.size();
    for (int i = 0; i < nref; i++) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsNeighborSol = false;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = false;
    }
}

void TPZMixedElasticityDynamics::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int nref = datavec.size();
    for (int i = 0; i < nref; i++) {
        datavec[i].fNeedsNormal = true;
    }
}

int TPZMixedElasticityDynamics::NStateVariables() const {
    return fDimension;
}

void TPZMixedElasticityDynamics::Solution(const TPZVec<TPZMaterialDataT<STATE>> &data, int var, TPZVec<STATE> &Solout) {
#ifdef PZDEBUG
    // if (data.size() != 3) {
    //     DebugStop();
    // }
#endif
    TPZManVector<REAL, 3> x = data[0].x;
    REAL R = x[0];
    if (R < 1.e-6) {
        R = 1.e-6;
    }
    TPZFNMatrix<9, STATE> sigma(3, 3, 0.), sigmah(3, 3, 0.), eps(3, 3, 0.);
    int dim = Dimension();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < 3; j++) {
            if (fAxisSymmetric) {
                sigma(i, j) = data[0].sol[0][j + i * 3] / R;
            } else {
                sigma(i, j) = data[0].sol[0][j + i * 3];
            }
        }
    }
    TPZManVector<STATE, 3> disp(dim);
    if ( var == 9 && data[1].sol[0].size() != 0){
        for (int i = 0; i < dim; i++) {
            disp[i] = data[1].sol[0][i];
        }
    }
    
    TPZFNMatrix<9, STATE> antisym(dim, dim, 0.);
    antisym(0, 1) = data[2].sol[0][0];
    antisym(1, 0) = -antisym(0, 1);
    if(dim == 3)
    {
        antisym(0, 2) = data[2].sol[0][1];
        antisym(2, 0) = -antisym(0, 2);
        antisym(1, 2) = data[2].sol[0][2];
        antisym(2, 1) = -antisym(1, 2);
    }

    TElasticityAtPoint elast(fE_const,fnu_const);
    if(fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
		TPZManVector<STATE, 3> result(dim);
        TPZFNMatrix<9,STATE> Dres(0,0);
        fElasticity(x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }
    
    REAL mu = elast.fmu;
    REAL E = elast.fE;
    REAL nu = elast.fnu;
    
    if(var == 28)
    {
        Solout[0] = E;
        return ;
    }
    if(var == 29)
    {
        Solout[0] = elast.fnu;
        return;
    }
    

    REAL Pressure;

    //TPZManVector<REAL, 4> SIGMA(4, 0.), EPSZ(4, 0.);
	TPZManVector<STATE, 9> SIGMA(dim*dim, 0.), EPSZ(dim*dim, 0.);

    ToVoigt(sigma, SIGMA);

    ComputeDeformationVector(SIGMA, EPSZ, elast);

    FromVoigt(EPSZ, eps);

    if(dim ==2)
    {
        if (this->fPlaneStress == 1) {
            eps(2, 2) = -1. / E * nu * (sigma(0, 0) + sigma(1, 1));
        } else {
            sigma(2, 2) = nu * (sigma(0, 0) + sigma(1, 1));
            eps(2, 2) = 0;
        }
    }
    Pressure = -1/3. * (sigma(0, 0) + sigma(1, 1) + sigma(2,2));
    sigmah = sigma;
    sigmah(0, 0) += Pressure;
    sigmah(1, 1) += Pressure;
    sigmah(2, 2) += Pressure;
    // Displacement
    if (var == 9) {
        Solout[2] = 0.;
        for(int idf = 0; idf < dim; idf++)
        {
            Solout[idf] = disp[idf];
        }
        return;
    }
    // Pressure
    if (var == 1) {
        Solout[0] = Pressure;
        return;
    }
    // Sigmaz
    if (var == 12) {
        Solout[0] = sigma(2, 2);

        return;
    }
    // SigmaY
    if (var == 6) {
        Solout[0] = sigma(1, 1);
        return;
    }

    // SigmaX                
    if (var == 5) {
        Solout[0] = sigma(0, 0);
        return;
    }
    //  TauXY
    if (var == 8) {
        Solout[0] = 0.5 * (sigma(0, 1) + sigma(1, 0));
        return;
    }
    //  TauXZ
    if (var == 30) {
        Solout[0] = 0.5 * (sigma(0, 2) + sigma(2, 0));
        return;
    }
    //  TauYZ
    if (var == 31) {
        Solout[0] = 0.5 * (sigma(1, 2) + sigma(2, 1));
        return;
    }

    // Exact displacement                
    if (var == 33) {
        TPZVec<STATE> u_exact(fDimension,0.);
        TPZFMatrix<STATE> du_exact(fDimension,fDimension,0.);
        if (this->fExactSol) {
            this->fExactSol(data[0].x, u_exact, du_exact);
        }
        for (int idf = 0; idf < dim; idf++) {
            Solout[idf] = u_exact[idf];
        }
        return;
    }
    // Exact stress                
    if (var == 34) {
        
        TPZVec<STATE> u_exact(fDimension,0.);
        TPZFMatrix<STATE> du_exact(fDimension,fDimension,0.);
        if (this->fExactSol) {
            this->fExactSol(data[0].x, u_exact, du_exact);
        }
        // std::cout << "duexact = " << du_exact << std::endl;
        // For 2D only
        Solout[0] = du_exact(0,0);//Sigma x
        Solout[1] = du_exact(1,1);//Sigma y
        Solout[2] = (du_exact(1,0)+du_exact(0,1))/2;//Tau xy


        return;
    }




    //Strain
    if (var == 11) {
        Solout[Exx] = eps(0, 0);
        Solout[Eyx] = eps(1, 0);
        Solout[Exy] = eps(0, 1);
        Solout[Eyy] = eps(1, 1);
        if(dim == 3)
        {
            Solout[Exz] = eps(0, 2);
            Solout[Eyz] = eps(1, 2);
            Solout[Ezx] = eps(2, 0);
            Solout[Ezy] = eps(2, 1);
            Solout[Ezz] = eps(2, 2);
        }
        return;
    }

    //Stress
    if (var == 10) {
        Solout[Exx] = sigma(0, 0);
        Solout[Eyx] = sigma(1, 0);
        Solout[Exy] = sigma(0, 1);
        Solout[Eyy] = sigma(1, 1);
        if(dim == 3)
        {
            Solout[Exz] = sigma(0, 2);
            Solout[Eyz] = sigma(1, 2);
            Solout[Ezx] = sigma(2, 0);
            Solout[Ezy] = sigma(2, 1);
            Solout[Ezz] = sigma(2, 2);
        }
        return;
    }

    //I1
    if (var == 21) {
        Solout[0] = sigma(0, 0) + sigma(1, 1) + sigma(2, 2);
        return;
    }
    //J2
    if (var == 20) {
        Solout[0] = sigmah(1, 1) * sigmah(0, 0) + sigmah(1, 1) * sigma(2, 2)
         + sigmah(0,0) * sigmah(2,2) - sigmah(1, 0) * sigmah(1, 0) - sigmah(0, 1) * sigmah(0, 1);
        return;
    }
    //NormalStrain?

    //PrincipalStrain1
    if (var == 3) {
        Solout[0] = Pressure + sqrt(0.25 * (sigma(0, 0) - sigma(1, 1))*(sigma(0, 0) - sigma(1, 1)) + sigma(1, 2));
        return;
    }
    //Rotation
    if (var == 27) {
        Solout[0] = antisym(0, 1);
        if(dim == 3)
        {
            Solout[1] = antisym(0, 2);
            Solout[2] = antisym(1, 2);
        }
        return;
    }
}

void TPZMixedElasticityDynamics::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) {
    int nspaces = datavec.size();
    switch(nspaces){
        case 3:
            Contribute_3spaces(datavec, weight, ek, ef);
            break;
        case 5:
        case 7:
            Contribute_5spaces(datavec, weight, ek, ef);
            break;
        default:
            DebugStop();
    }    

}

/** @brief Calculates the element stiffness matrix using 3 spaces - Stress tensor, displacement, and skew-symmetric tensor (for weak symmetry) */
void TPZMixedElasticityDynamics::Contribute_3spaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){

    int nshapeS, nshapeU, nshapeP;
    nshapeS = datavec[0].fDeformedDirections.Cols();
    nshapeU = datavec[1].phi.Rows();
    nshapeP = datavec[2].phi.Rows();
    const int firstequation_S = 0;
    const int firstequation_U = firstequation_S + nshapeS*fDimension;
    const int firstequation_P = firstequation_U + nshapeU*fDimension;

    // Number of voight terms
    int voigtdim = fDimension*fDimension;

    TElasticityAtPoint elast(fE_const, fnu_const);
    if (fElasticity) {
        TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4, STATE> Dres(0, 0);
        fElasticity(datavec[0].x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        elast = TElasticityAtPoint(E, nu);
    }

    // number of asymetric tensors for each shape function
    int nrotations = (fDimension == 3) ? 3 : 1;

    //K11 - (Matrix A * stress tensor) x test-function stress tensor
    TPZFNMatrix<200, STATE> PhiSVoight(voigtdim, nshapeS*fDimension, 0.);

    for (int i = 0; i < nshapeS; i++) {
        for (int k = 0; k < fDimension; k++) {
            TPZFNMatrix<9,STATE> phiCrossX(fDimension, fDimension, 0.);
            for (int e = 0; e < fDimension; e++) {
                phiCrossX(k, e) = datavec[0].fDeformedDirections(e, i);
            }

            TPZManVector<STATE, 9>  phiCrossXVoight(voigtdim, 0.0);
            ToVoigt(phiCrossX, phiCrossXVoight);
            for (int l = 0; l < voigtdim; l++) {
                PhiSVoight(l, i*fDimension + k) = phiCrossXVoight[l];
            }
        }
    }

    TPZFNMatrix<200, STATE> APhiSVoight(voigtdim, nshapeS*fDimension, 0.);
    TPZFNMatrix<81, STATE> MatrixElast(voigtdim, voigtdim, 0.);
    ElasticityModulusTensor(MatrixElast, elast);
    MatrixElast.Multiply(PhiSVoight, APhiSVoight);

    REAL factor = weight;
    ek.AddContribution(firstequation_S, firstequation_S, PhiSVoight, true, APhiSVoight, false, weight);

    //K21 and K12 - divergent of test-function stress tensor * displacement vector
    TPZFNMatrix<200, STATE> PhiU(fDimension, nshapeU*fDimension, 0.);
    for (int i = 0; i < nshapeU; i++) {
        for (int k = 0; k < fDimension; k++) {
            PhiU(k, i*fDimension + k) = datavec[1].phi(i, 0);
        }
    }

    TPZFNMatrix<200, STATE> DivPhiS(fDimension, nshapeS*fDimension, 0.);
    for (int i = 0; i < nshapeS; i++) {
        for (int k = 0; k < fDimension; k++) {
            DivPhiS(k, i*fDimension + k) = datavec[0].divphi(i);
        }
    }

    ek.AddContribution(firstequation_S, firstequation_U, DivPhiS, true, PhiU, false, factor);
    ek.AddContribution(firstequation_U, firstequation_S, PhiU, true, DivPhiS, false, factor);

    //K22 - test function displacement vector x displacement vector (mass matrix)
    factor *= -fRho;
    ek.AddContribution(firstequation_U, firstequation_U, PhiU, true, PhiU, false, factor);

    //K31 and K13 - test-function stress tensor x rotation tensor p
    TPZFNMatrix<200, STATE> PhiPVoight(voigtdim, nshapeP*nrotations, 0.);
    for (int i = 0; i < nshapeP; i++) {
        int cont = 0; // counter for the number of rotations
        for (int k = 0; k < fDimension-1; k++) {
            for (int l = k+1; l < fDimension; l++) {
                TPZFNMatrix<9,STATE> phiPTensor(fDimension, fDimension, 0.);
                phiPTensor(k, l) = datavec[2].phi(i, 0);
                phiPTensor(l, k) = -datavec[2].phi(i, 0);
                TPZManVector<STATE, 9> phiPTensorVoigt(voigtdim, 0.0);
                ToVoigt(phiPTensor, phiPTensorVoigt);
                for (int e = 0; e < voigtdim; e++)
                {
                    PhiPVoight(e, i * nrotations + cont) = phiPTensorVoigt[e];
                }
                cont++;
            }
        }
    }

    factor = weight;
    ek.AddContribution(firstequation_S, firstequation_P, PhiSVoight, true, PhiPVoight, false, factor);
    ek.AddContribution(firstequation_P, firstequation_S, PhiPVoight, true, PhiSVoight, false, factor);

    //Body force contribution
    TPZFNMatrix<3,STATE>  force(fDimension,1, 0.);
    if (this->HasForcingFunction()) {
        fForcingFunction(datavec[0].x, fForce);
    }
    for (int i = 0; i < fDimension; i++) {
        force(i, 0) = fForce[i];
    }

    ef.AddContribution(firstequation_U,0, PhiU, true, force, false, -factor);
}

/** @brief Calculates the element stiffness matrix using 3 spaces - Stress tensor, displacement, and skew-symmetric tensor (for weak symmetry) */
void TPZMixedElasticityDynamics::Contribute_5spaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){

    int nspaces = datavec.size();
    Contribute_3spaces(datavec, weight, ek, ef);

    const int nshapeS = datavec[0].fDeformedDirections.Cols();
    const int nshapeU = datavec[1].phi.Rows();
    const int nshapeP = datavec[2].phi.Rows();
    const int nrotations = fDimension == 3 ? 3 : 1;
    const int firstequation_S = 0;
    const int firstequation_U = firstequation_S + nshapeS * fDimension;
    const int firstequation_P = firstequation_U + nshapeU * fDimension;
    int firstequation_FRB = firstequation_P + nshapeP * nrotations;

    ContributeRigidBodyMode(datavec, weight, ek, ef, firstequation_FRB, 3);
    
    if (nspaces == 5) return;

    //MHM maybe has 7 spaces?
    // if (fDimension == 2) {
    //     firstequation_FRB += 3;
    // }
    // else {
    //     firstequation_FRB += 12;
    // }
    firstequation_FRB += fDimension*nrotations + nrotations;
    ContributeRigidBodyMode(datavec, weight, ek, ef, firstequation_FRB, 5);
}

void TPZMixedElasticityDynamics::ContributeRigidBodyMode(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, int firstRBEq, int firstRBspace){

    const int nshapeS = datavec[0].fDeformedDirections.Cols();
    const int nshapeU = datavec[1].phi.Rows();
    const int nshapeP = datavec[2].phi.Rows();
    const int nshapeFRB = datavec[firstRBspace].phi.Rows();
    const int nshapeURB = datavec[firstRBspace+1].phi.Rows();
    const int nrotations = (fDimension == 3) ? 3 : 1;
    const int firstequation_S = 0;
    const int firstequation_U = firstequation_S + nshapeS*fDimension;
    const int firstequation_P = firstequation_U + nshapeU*fDimension;
    const int firstequation_FRB = firstRBEq;
    const int ncomponents = fDimension+nrotations;
    const int firstequation_URB = firstequation_FRB + nshapeFRB*ncomponents;

    // get distance of integration point to center of element

    TPZManVector<REAL,3> xcenter = datavec[firstRBspace].XCenter;
    TPZManVector<REAL,3> x = datavec[0].x;
    TPZManVector<REAL,3> delx(3,0.);
    for (int i = 0; i < fDimension; i++) {
        delx[i] = x[i] - xcenter[i];
    }

    
}

void TPZMixedElasticityDynamics::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    if (datavec[0].phi.Rows() != 0 && datavec[0].fShapeType != TPZMaterialData::EScalarShape) {
        DebugStop();
    }
    int ndisp = datavec[1].phi.Rows();

    const TPZVec<REAL> v_2 = bc.Val2();
    TPZFNMatrix<9, STATE> v_1 = bc.Val1();

    // Setting forcing function
    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE, 3> res(fDimension);
        TPZFNMatrix<9, STATE> tens(fDimension, fDimension);
        bc.ForcingFunctionBC()(datavec[0].x, res, tens);
        TPZFNMatrix<9, STATE> strain = tens;
        for (int i = 0; i < fDimension; i++) {
            for (int j = 0; j < fDimension; j++) {
                strain(i, j) += tens(j,i);
            }
        }
        strain *= 0.5;
        v_2[0] = res[0];
        v_2[1] = res[1];
        if(fDimension == 3) v_2[2] = res[2];
        if(bc.Type() == 1) { // Neuman condition
            TElasticityAtPoint elast(fE_const, fnu_const);
            if (fElasticity) {
                TPZManVector<STATE, 3> result(2);
                TPZFNMatrix<4, STATE> Dres(0, 0);
                fElasticity(datavec[0].x, result, Dres);
                REAL E = result[0];
                REAL nu = result[1];
                elast = TElasticityAtPoint(E, nu);
            }
            int nterms = fDimension*fDimension;
            TPZManVector<STATE,9> deform(nterms,0.),stressvec(nterms,0.);
            ToVoigt(strain, deform);
            ComputeStressVector(deform, stressvec, elast);
            FromVoigt(stressvec, tens);
            const TPZManVector<STATE, 3>& normal = datavec[0].normal;
            for (int i = 0; i < fDimension; i++)
                v_2[i] = 0.0;
            for (int i = 0; i < fDimension; i++)
                for (int j = 0; j < fDimension; j++)
                    v_2[i] += tens(i,j) * normal[j];
        }
    }

    // Setting the phis
    // E
    TPZFMatrix<REAL> &phiS = datavec[0].phi;

    int nshapeS;
    nshapeS = datavec[0].phi.Rows();


    //    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    //    if(shapetype==data.EVecShape){
    //        ContributeVecShapeBC(data,weight,ek, ef,bc);
    //        return;
    //    }

    REAL R = datavec[0].x[0];
    if (R < 1.e-6) R = 1.e-6;
    
    int nstate = 2;
    if(fDimension == 3) nstate = 3;

    switch (bc.Type()) {
        case 0: // Dirichlet condition
        {
            for (int iq = 0; iq < nshapeS; iq++) {
                for (int idf = 0; idf < nstate; idf++) {
                    ef(nstate * iq + idf, 0) += v_2[idf] * phiS(iq, 0) * weight; // forced v2 displacement
                }
            }
        }
            break;

        case 1: // Neumann condition
        {
            for (int iq = 0; iq < nshapeS; iq++) {
                for (int jq = 0; jq < nshapeS; jq++) {
                    if (fAxisSymmetric) {
                        ek(2 * iq, 2 * jq) += TPZMaterial::fBigNumber * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                        ek(2 * iq + 1, 2 * jq + 1) += TPZMaterial::fBigNumber * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                    } else {
                        for(int idf = 0; idf < nstate; idf++)
                        {
                            ek(nstate * iq + idf, nstate * jq + idf) += TPZMaterial::fBigNumber * phiS(iq, 0) * phiS(jq, 0) * weight;
                        }
                    }
                }
                for(int idf = 0; idf < nstate; idf++)
                {
                    ef(nstate * iq + idf, 0) += TPZMaterial::fBigNumber * v_2[idf] * phiS(iq, 0) * weight; // normal stress in x direction
                }
            }
        }
            break;

        case 2: // Mixed condition
        {
            for (int iq = 0; iq < nshapeS; iq++) {
                for (int jq = 0; jq < nshapeS; jq++) {
                    if (fAxisSymmetric) {
                        ek(2 * iq, 2 * jq) += v_1(0, 0) * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                        ek(2 * iq + 1, 2 * jq + 1) += v_1(1, 1) * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                        ek(2 * iq + 1, 2 * jq) += v_1(1, 0) * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                        ek(2 * iq, 2 * jq + 1) += v_1(0, 1) * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                    } else {
                        for (int idf = 0; idf < nstate; idf++)
                        {
                            for (int jdf = 0; jdf < nstate; jdf++)
                            {
                                ek(nstate * iq + idf, nstate * jq + jdf) += v_1(idf, jdf) * phiS(iq, 0) * phiS(jq, 0) * weight;
                            }
                        }
                    }
                }
                for (int idf = 0; idf < nstate; idf++) {
                    ef(nstate * iq + idf, 0) += v_2[idf] * phiS(iq, 0) * weight; // normal stress
                }
            }
        }
            break;
        case 5:
            /// aplicando multiplicador de lagrange para impor deslocamento
            if (ndisp != 1) {
                DebugStop();
            }
            for (int idf = 0; idf < nstate; idf++)
            {
                if (v_1(idf, idf) > 1.e-1) {
                    ek(nstate * nshapeS + idf, nstate * nshapeS + idf) = 1.;
                }
            }
            break;
        default:
            DebugStop();
            // nulo introduzindo o BIGNUMBER pelos valores da condição
    } // 1 Val1 : a leitura é 00 01 10 11
    
}