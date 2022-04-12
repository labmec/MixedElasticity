/**
 * @file
 * @brief Contains implementations of the TPZMixedElasticityNDAiry methods.
 */

#include "TPZMixedElasticityNDAiry.h"
#include "TPZMaterialDataT.h"
#include "TPZAnalyticSolution.h"

TPZMixedElasticityNDAiry::TPZMixedElasticityNDAiry() : TBase() {
    fE_const = -1.; // Young modulus
    fnu_const = -1.; // Poisson coefficient
    fForce[0] = 0.; // X component of the body force
    fForce[1] = 0.; // Y component of the body force
    fForce[2] = 0.; // Z component of the body force - not used for this class
    flambda_const = 0.;
    fmu_const = 0.;

    fPlaneStress = 0;
    fMatrixA = 0.;
}


TPZMixedElasticityNDAiry::TPZMixedElasticityNDAiry(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress, int dimension) : TBase(id) {
    this->SetElasticity(E, nu);
    fDimension = dimension;
    fForce[0] = fx; // X component of the body force
    fForce[1] = fy; // Y component of the body force
    fForce[2] = 0.;
    fPlaneStress = planestress;
}

TPZMixedElasticityNDAiry::TPZMixedElasticityNDAiry(const TPZMixedElasticityNDAiry &copy) : TBase(copy)
{
    fE_const = copy.fE_const;
    fnu_const = copy.fnu_const;
    flambda_const = copy.flambda_const;
    fmu_const = copy.fmu_const;
    fForce = copy.fForce;
    fPlaneStress = copy.fPlaneStress;
    fDimension = copy.fDimension;
    fMatrixA = copy.fMatrixA;
    fElasticity = copy.fElasticity;
    fAxisSymmetric = copy.fAxisSymmetric;
    fExactSol = copy.fExactSol;
    fExactPOrder = copy.fExactPOrder;
}

TPZMixedElasticityNDAiry::~TPZMixedElasticityNDAiry(){

}

void TPZMixedElasticityNDAiry::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    if (datavec[0].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec[0]);
    }
    int nspaces = datavec.size();
    switch(nspaces){
        case 4:
            Contribute_4spaces(datavec, weight, ek, ef);
            break;
        default:
            DebugStop();
    }    
}

/** @brief Calculates the element stiffness matrix using 3 spaces - Stress tensor, displacement, and skew-symmetric tensor (for weak symmetry) */
void TPZMixedElasticityNDAiry::Contribute_4spaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    REAL R = datavec[0].x[0];
    // Setting the phi's
    // E
    TPZFMatrix<REAL> &phiS = datavec[0].phi;
    TPZFMatrix<REAL> &dphiS = datavec[0].dphix;
    // Airy
    TPZFMatrix<REAL> &phiAiry = datavec[1].phi;
    TPZFMatrix<REAL> &dphiAiry = datavec[1].dphix;
    // U
    TPZFMatrix<REAL> &phiU = datavec[2].phi;
    TPZFMatrix<REAL> &dphiU = datavec[2].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[3].phi;

    TElasticityAtPoint elast(fE_const,fnu_const);
    if(this->fElasticity)
    {
		TPZManVector<STATE, 3> result(2);
		TPZFNMatrix<4,STATE> Dres(0,0);
        this->fElasticity->Execute(datavec[0].x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }
    

    int nshapeS, nshapeA, nshapeU, nshapeP;
    nshapeS = datavec[0].fVecShapeIndex.NElements();
    nshapeA = datavec[1].phi.Rows();
    nshapeU = datavec[2].phi.Rows();
    nshapeP = datavec[3].phi.Rows();
    const int firstequation_S = 0;
    const int firstequation_A = firstequation_S + nshapeS*fDimension;
    const int firstequation_U = firstequation_A + nshapeA*fDimension;
    const int firstequation_P = firstequation_U + nshapeU*fDimension;
    
    // number of asymetric tensors for each shape function
    int nrotations = 1;
    if(fDimension == 3) nrotations = 3;

    TPZManVector<STATE, 3>  force(fDimension, 0.);
    
    TPZFNMatrix<3,STATE>    phiSi(3, 1, 0.), 
                            phiSj(fDimension, 1, 0.);
    
    TPZManVector<STATE, 3>  divSi1x(fDimension, 0.), 
                            divSi1y(fDimension, 0.), 
                            divSi1z(fDimension, 0.);
    
    int voigtdim = fDimension*fDimension;
    TPZManVector<STATE, 9>  phiSi1x(voigtdim, 0.0), phiSi1y(voigtdim, 0.0), phiSi1z(voigtdim, 0.0),
                            phiAi1x(voigtdim, 0.0), phiAi1y(voigtdim, 0.0), phiAi1z(voigtdim, 0.0),
                            phiSj1x(voigtdim, 0.0), phiSj1y(voigtdim, 0.0), phiSj1z(voigtdim, 0.0),
                            phiAj1x(voigtdim, 0.0), phiAj1y(voigtdim, 0.0), phiAj1z(voigtdim, 0.0),
                            phiPj1x(voigtdim, 0.0), phiPj1y(voigtdim, 0.0), phiPj1z(voigtdim, 0.0);
    
    TPZManVector<STATE, 3>  phiUi(fDimension, 0.0), 
                            phiUj(fDimension, 0.0), 
                            phiUj1x(fDimension, 0.0), 
                            phiUj1y(fDimension, 0.0), 
                            phiUj1z(fDimension, 0.0);;
    
    TPZFNMatrix<3,STATE>    phiPi(fDimension, 1, 0.0), 
                            phiPj(fDimension, 1, 0.0);
    
    TPZFNMatrix<3, REAL>    ivecS(3, 1, 0.);

    force = fForce;
    if (this->HasForcingFunction()) {
        fForcingFunction(datavec[0].x, force);
#ifdef LOG4CXX
        if (logdata->isDebugEnabled()) {
            std::stringstream sout;
            sout << " x = " << datavec[0].x << " force = " << force << std::endl;
            LOGPZ_DEBUG(logdata, sout.str())
        }
#endif
    }
//    datavec[0].ComputeFunctionDivergence();
    for (int i = 0; i < nshapeS; i++) {
        int iphi = datavec[0].fVecShapeIndex[i].second;
        int ivec = datavec[0].fVecShapeIndex[i].first;
        TPZFNMatrix<4> GradSi(fDimension, fDimension);

        REAL divSi = 0.;

        for (int e = 0; e < 3; e++) {

            ivecS(e, 0) = datavec[0].fDeformedDirections(e, ivec);
            phiSi(e, 0) = phiS(iphi, 0) * ivecS(e, 0);

        }
//        TPZFNMatrix<3, REAL> axesvec(fDimension, 1, 0.);
//        datavec[0].axes.Multiply(ivecS, axesvec);
//
//        //calculando div(Si)
//        for (int f = 0; f < fDimension; f++) {
//            divSi += axesvec(f, 0) * dphiS(f, iphi);
//        }
        divSi = datavec[0].divphi(i);

        TPZFNMatrix<9, STATE> phiTensx(fDimension, fDimension, 0.), phiTensy(fDimension, fDimension, 0.),
        phiTensz(fDimension,fDimension,0.);

        for(int d = 0; d< fDimension; d++)
        {
            phiTensx(0,d) = phiSi(d,0);
            phiTensy(1,d) = phiSi(d,0);
            if(fDimension == 3)
            {
                phiTensz(2,d) = phiSi(d,0);
            }
        }
        ToVoigt(phiTensx, phiSi1x);
        ToVoigt(phiTensy, phiSi1y);
        if(fDimension == 3)
        {
            ToVoigt(phiTensz, phiSi1z);
        }

        divSi1x[0] = divSi;
        divSi1y[1] = divSi;
        if(fDimension == 3) divSi1z[2] = divSi;

        // matrix K11 - (Matrix A * stress tensor) x test-function stress tensor
        for (int j = 0; j < nshapeS; j++) {
            int jphi = datavec[0].fVecShapeIndex[j].second;
            int jvec = datavec[0].fVecShapeIndex[j].first;

            for (int e = 0; e < fDimension; e++) {
                phiSj(e, 0) = phiS(jphi, 0) * datavec[0].fDeformedDirections(e, jvec);
            }

            
            TPZFNMatrix<9, STATE> phjTensx(fDimension, fDimension, 0.), phjTensy(fDimension, fDimension, 0.)
            , phjTensz(fDimension, fDimension, 0.);

            for(int d = 0; d< fDimension; d++)
            {
                phjTensx(0,d) = phiSj(d,0);
                phjTensy(1,d) = phiSj(d,0);
                if(fDimension == 3)
                {
                    phjTensz(2,d) = phiSj(d,0);
                }
            }
            ToVoigt(phjTensx, phiSj1x);
            ToVoigt(phjTensy, phiSj1y);
            if(fDimension == 3)
            {
                ToVoigt(phjTensz, phiSj1z);
            }

            //Multiply by Lamé parameters
            TPZManVector<STATE, 9> AphiSi1x(voigtdim, 0.0), AphiSi1y(voigtdim, 0.0), AphiSi1z(voigtdim, 0.0);

            ComputeDeformationVector(phiSi1x, AphiSi1x,elast);
            ComputeDeformationVector(phiSi1y, AphiSi1y,elast);
            if(fDimension == 3)
            {
                ComputeDeformationVector(phiSi1z, AphiSi1z,elast);
            }

            STATE valxx = InnerVec(AphiSi1x, phiSj1x);
            STATE valxy = InnerVec(AphiSi1x, phiSj1y);
            STATE valyx = InnerVec(AphiSi1y, phiSj1x);
            STATE valyy = InnerVec(AphiSi1y, phiSj1y);

            //Matrix K11
            if (fAxisSymmetric) {
                valxx /= R;
                valxy /= R;
                valyx /= R;
                valyy /= R;
            }
            ek(fDimension * i    , fDimension * j    ) += weight * valxx;
            ek(fDimension * i    , fDimension * j + 1) += weight * valxy;
            ek(fDimension * i + 1, fDimension * j    ) += weight * valyx;
            ek(fDimension * i + 1, fDimension * j + 1) += weight * valyy;
            if(fDimension == 3)
            {
                STATE valxz = InnerVec(AphiSi1x, phiSj1z);
                ek(fDimension * i    , fDimension * j + 2) += weight * valxz;
                STATE valyz = InnerVec(AphiSi1y, phiSj1z);
                ek(fDimension * i + 1, fDimension * j + 2) += weight * valyz;
                STATE valzx = InnerVec(AphiSi1z, phiSj1x);
                ek(fDimension * i + 2, fDimension * j    ) += weight * valzx;
                STATE valzy = InnerVec(AphiSi1z, phiSj1y);
                ek(fDimension * i + 2, fDimension * j + 1) += weight * valzy;
                STATE valzz = InnerVec(AphiSi1z, phiSj1z);
                ek(fDimension * i + 2, fDimension * j + 2) += weight * valzz;

            }
        }


        //Here we have to put the airy functions - Jeferson
        for (int j = 0; j < nshapeA; j++) {
            phiAj1x[Exx] = phiAiry(2*j  ,0);
            phiAj1x[Exy] = phiAiry(2*j  ,1);    
            phiAj1x[Eyy] = phiAiry(2*j  ,2);
            phiAj1x[Eyy] = phiAiry(2*j  ,3);

            phiAj1y[Exx] = phiAiry(2*j  ,0);
            phiAj1y[Exy] = phiAiry(2*j  ,1);    
            phiAj1y[Eyy] = phiAiry(2*j  ,2);
            phiAj1y[Eyy] = phiAiry(2*j  ,3);

            //Multiply by Lamé parameters
            TPZManVector<STATE, 9> AphiSi1x(voigtdim, 0.0), AphiSi1y(voigtdim, 0.0), AphiSi1z(voigtdim, 0.0);

            ComputeDeformationVector(phiSi1x, AphiSi1x,elast);
            ComputeDeformationVector(phiSi1y, AphiSi1y,elast);

            STATE valxx = InnerVec(AphiSi1x, phiAj1x);
            STATE valxy = InnerVec(AphiSi1x, phiAj1y);
            STATE valyx = InnerVec(AphiSi1y, phiAj1x);
            STATE valyy = InnerVec(AphiSi1y, phiAj1y);

            ek(firstequation_A + fDimension * j    , fDimension * i    ) += weight * valxx;
            ek(firstequation_A + fDimension * j    , fDimension * i + 1) += weight * valxy;
            ek(firstequation_A + fDimension * j + 1, fDimension * i    ) += weight * valyx;
            ek(firstequation_A + fDimension * j + 1, fDimension * i + 1) += weight * valyy;

            ek(fDimension * i    , firstequation_A + fDimension * j    ) += weight * valxx;
            ek(fDimension * i    , firstequation_A + fDimension * j + 1) += weight * valxy;
            ek(fDimension * i + 1, firstequation_A + fDimension * j    ) += weight * valyx;
            ek(fDimension * i + 1, firstequation_A + fDimension * j + 1) += weight * valyy;
        }   

        // matrix K21 and K12 - divergent of test-function stress tensor * displacement vector
        for (int j = 0; j < nshapeU; j++) {
            phiUj1x[0] = phiU(j, 0);
            phiUj1y[1] = phiU(j, 0);
            
            STATE valx = weight * InnerVec(divSi1x, phiUj1x);
            STATE valy = weight * InnerVec(divSi1y, phiUj1y);

            //position K21
            ek(fDimension * j     + firstequation_U, fDimension * i    ) += valx;
            ek(fDimension * j + 1 + firstequation_U, fDimension * i + 1) += valy;

            //position K12
            ek(fDimension * i    , fDimension * j     + firstequation_U) += valx;
            ek(fDimension * i + 1, fDimension * j + 1 + firstequation_U) += valy;
            
            if(fDimension == 3)
            {
                phiUj1z[2] = phiU(j,0);
                REAL valz = weight * InnerVec(divSi1z, phiUj1z);
                ek(fDimension * j + 2 + firstequation_U, fDimension * i + 2) += valz;
                ek(fDimension * i + 2, fDimension * j + 2 + firstequation_U) += valz;
            }

        }

        // matrix K31 and K13 - test-function stress tensor x rotation tensor p
        for (int j = 0; j < nshapeP; j++) {
            TPZFNMatrix<9, STATE> phiPTensor(fDimension, fDimension, 0.);
            phiPTensor(0, 1) = phiP(j, 0);
            phiPTensor(1, 0) = -phiP(j, 0);
            ToVoigt(phiPTensor, phiPj1x);

            STATE valxx = InnerVec(phiSi1x, phiPj1x);
            STATE valxy = InnerVec(phiSi1y, phiPj1x);
            //Matrix K31
            ek(nrotations*j + firstequation_P, fDimension * i    ) += weight * valxx;
            ek(nrotations*j + firstequation_P, fDimension * i + 1) += weight * valxy;

            //Matrix K13
            ek(fDimension * i    , nrotations*j + firstequation_P) += weight * valxx;
            ek(fDimension * i + 1, nrotations*j + firstequation_P) += weight * valxy;
            
            if(fDimension == 3)
            {
                STATE valxz = InnerVec(phiSi1z, phiPj1x);
                ek(nrotations*j + firstequation_P, fDimension * i + 2) += weight * valxz;
                ek(fDimension * i + 2, nrotations*j + firstequation_P) += weight * valxz;
                for(int d=0; d<2; d++)
                {
                    phiPTensor.Zero();
                    phiPTensor(d, 2) = phiP(j, 0);
                    phiPTensor(2, d) = -phiP(j, 0);
                    ToVoigt(phiPTensor, phiPj1x);

                    STATE valxx = InnerVec(phiSi1x, phiPj1x);
                    STATE valxy = InnerVec(phiSi1y, phiPj1x);
                    STATE valxz = InnerVec(phiSi1z, phiPj1x);

                    //Matrix K31
                    ek(nrotations*j + d + 1 + firstequation_P, fDimension * i) += weight * valxx;
                    ek(nrotations*j + d + 1 + firstequation_P, fDimension * i + 1) += weight * valxy;
                    ek(nrotations*j + d + 1 + firstequation_P, fDimension * i + 2) += weight * valxz;

                    //Matrix K13
                    ek(fDimension * i, nrotations*j + d + 1 + firstequation_P) += weight * valxx;
                    ek(fDimension * i + 1, nrotations*j + d + 1 + firstequation_P) += weight * valxy;
                    ek(fDimension * i + 2, nrotations*j + d + 1 + firstequation_P) += weight * valxz;

                }
            }
        }
    }

    //Diagonal contribution of Airy functions
    for (int i = 0; i < nshapeA; i++) {
        phiAi1x[Exx] = phiAiry(2*i  ,0);
        phiAi1x[Exy] = phiAiry(2*i  ,1);    
        phiAi1x[Eyy] = phiAiry(2*i  ,2);
        phiAi1x[Eyy] = phiAiry(2*i  ,3);
        
        phiAi1y[Exx] = phiAiry(2*i  ,0);
        phiAi1y[Exy] = phiAiry(2*i  ,1);    
        phiAi1y[Eyy] = phiAiry(2*i  ,2);
        phiAi1y[Eyy] = phiAiry(2*i  ,3);

        for (int j = 0; j < nshapeA; j++) {
            phiAj1x[Exx] = phiAiry(2*j  ,0);
            phiAj1x[Exy] = phiAiry(2*j  ,1);    
            phiAj1x[Eyy] = phiAiry(2*j  ,2);
            phiAj1x[Eyy] = phiAiry(2*j  ,3);

            phiAj1y[Exx] = phiAiry(2*j  ,0);
            phiAj1y[Exy] = phiAiry(2*j  ,1);    
            phiAj1y[Eyy] = phiAiry(2*j  ,2);
            phiAj1y[Eyy] = phiAiry(2*j  ,3);

            //Multiply by Lamé parameters
            TPZManVector<STATE, 9> AphiAi1x(voigtdim, 0.0), AphiAi1y(voigtdim, 0.0), AphiAi1z(voigtdim, 0.0);

            ComputeDeformationVector(phiAi1x, AphiAi1x,elast);
            ComputeDeformationVector(phiAi1y, AphiAi1y,elast);

            STATE valxx = InnerVec(AphiAi1x, phiAj1x);
            STATE valxy = InnerVec(AphiAi1x, phiAj1y);
            STATE valyx = InnerVec(AphiAi1y, phiAj1x);
            STATE valyy = InnerVec(AphiAi1y, phiAj1y);

            ek(firstequation_A + fDimension * j    , firstequation_A + fDimension * i    ) += weight * valxx;
            ek(firstequation_A + fDimension * j    , firstequation_A + fDimension * i + 1) += weight * valxy;
            ek(firstequation_A + fDimension * j + 1, firstequation_A + fDimension * i    ) += weight * valyx;
            ek(firstequation_A + fDimension * j + 1, firstequation_A + fDimension * i + 1) += weight * valyy;
        }   
    }

    for (int i = 0; i < nshapeU; i++) {
        phiUj1x[0] = phiU(i, 0);
        phiUj1y[1] = phiU(i, 0);

        // Load vector f:
        STATE factfx = -weight * phiUj1x[0] * force[0];
        STATE factfy = -weight * phiUj1y[1] * force[1];
        if (fAxisSymmetric) {
            factfx *= R;
            factfy *= R;
        }
        //if(factfx != 0) DebugStop();
        ef(firstequation_U + fDimension * i, 0) += factfx;
        ef(firstequation_U + fDimension * i + 1, 0) += factfy;
        
        if(fDimension == 3)
        {
            phiUj1z[2] = phiU(i, 0);
            STATE factfz = -weight * phiUj1z[2] * force[2];
            ef(firstequation_U + fDimension * i + 2, 0) += factfz;

        }

        //if(ef(nshapeS*2+2*j) != 0.) DebugStop();
        if (fAxisSymmetric)
        {
            for (int j = 0; j < nshapeU; j++) {
                ek(firstequation_U + 2 * i, firstequation_U + 2 * j) -= weight * phiU(i, 0) * phiU(j, 0) / (elast.fE * R);
            }
        }
    }
}


int TPZMixedElasticityNDAiry::ClassId() const {
    return Hash("TPZMixedElasticityNDAiry") ^ TBase::ClassId() << 1;
}


/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMixedElasticityNDAiry::Solution(const TPZVec<TPZMaterialDataT<STATE>> &data, int var, TPZVec<STATE> &Solout) {
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
    for (int i = 0; i < dim; i++) {
        disp[i] = data[2].sol[0][i];
    }
    TPZFNMatrix<9, STATE> antisym(dim, dim, 0.);
    antisym(0, 1) = data[3].sol[0][0];
    antisym(1, 0) = -antisym(0, 1);
    if(dim == 3)
    {
        antisym(0, 2) = data[3].sol[0][1];
        antisym(2, 0) = -antisym(0, 2);
        antisym(1, 2) = data[3].sol[0][2];
        antisym(2, 1) = -antisym(1, 2);
    }

    TElasticityAtPoint elast(fE_const,fnu_const);
    if(fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
		TPZManVector<STATE, 3> result(dim);
        TPZFNMatrix<9,STATE> Dres(0,0);
        fElasticity->Execute(x, result, Dres);
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
    // SigmaX                
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
    //  TauXY
    if (var == 8) {
        Solout[0] = 0.5 * (sigma(0, 1) + sigma(1, 0));
        return;
    }
    // SigmaR              
    if (var == 30) {
        REAL r = sqrt(x[0]*x[0] + x[1]*x[1]);
        REAL theta = atan2(x[1],x[0]);
        Solout[0] = sigma(0,0)*cos(theta)*cos(theta) + 2.*sigma(0,1)*cos(theta)*sin(theta) + sigma(1,1)*sin(theta)*sin(theta);
        return;
    }
    // SigmaT         
    if (var == 31) {
        REAL r = sqrt(x[0]*x[0] + x[1]*x[1]);
        REAL theta = atan2(x[1],x[0]);
        Solout[0] = sigma(1,1)*cos(theta)*cos(theta) - 2.*sigma(0,1)*cos(theta)*sin(theta) + sigma(0,0)*sin(theta)*sin(theta);
        return;
    }
    // TauRT         
    if (var == 32) {
        REAL r = sqrt(x[0]*x[0] + x[1]*x[1]);
        REAL theta = atan2(x[1],x[0]);
        Solout[0] = sigma(0,1)*cos(theta)*cos(theta) + (sigma(1,1)-sigma(0,0))*cos(theta)*sin(theta) - sigma(0,1)*sin(theta)*sin(theta);
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


void TPZMixedElasticityNDAiry::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {
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
        for (int j = 0; j < dim; j++) {
            divSigma[i] += data[0].dsol[0](i * 3 + j, j);
        }
    }

    TPZManVector<STATE, 3> disp(dim);
    for (int i = 0; i < dim; i++) {
        disp[i] = data[2].sol[0][i];
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
        fElasticity->Execute(x, result, Dres);
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
    rotation[0] = data[3].sol[0][0];
    rotationExact[0] = 0.5*(du_exact(1, 0)-du_exact(0, 1));
    if(dim == 3)
    {
        rotation[1] = data[3].sol[0][1];
        rotationExact[1] = 0.5*(du_exact(2, 0)-du_exact(0, 2));
        rotation[2] = data[3].sol[0][2];
        rotationExact[2] = 0.5*(du_exact(1, 2)-du_exact(2, 1));

    }
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
    //    std::cout << "x " << data[0].x << std::endl;
    //    std::cout << "disp " << u_exact << std::endl;
    //    du_exact.Print("du ",std::cout);
    //    std::cout << "sigma_exact " << sigma_exactV << std::endl;
    //    std::cout << errors << std::endl;
    //or we can compute as this methods
    //TPZFMatrix<STATE> MatrixElast(4,4,0.);
    //ElasticityModulusTensor(MatrixElast);
    //errors[1] = 0;
    //for(int i==0; i<4； i++)
    //{
    //for(int j=0; j<4; j++)
    //{
    //errors[1]+=SIGMA[j]*MatrixElast(j,i);
    //}
    //}

    // SemiH1 norm
    //TPZFNMatrix<4,REAL> antisym(2,2,0.);
    //antisym(0,1) = data[2].sol[0][0];
    //antisym(1,0) = -antisym(0,1);
    //TPZManVector<REAL,4> P(4,0.),DU(4,0.);
    //ToVoigt(antisym,P);
    //ToVoigt(du_exact,DU);
    //errors[2]=0;
    //for(int i=0;i<4;i++)
    //{
    //   errors[2]+=(P[i]+EPSZ[i]-DU[i])*(P[i]+EPSZ[i]-DU[i]);
    //}
    //values[0] = calculo do erro estimado em norma Energia
    //values[0] = fE*(sigx*sigx + sigy*sigy + 2*fnu*sigx*sigy)/(1-fnu*fnu);
    //values[0] = (values[0] + .5*fE*sigxy*sigxy/(1+fnu));

    //values[1] : erro em norma L2 em tens�s
    //values[1] = sigx*sigx + sigy*sigy + sigxy*sigxy;

    //values[1] : erro em norma L2 em deslocamentos
    //values[1] = pow((REAL)fabs(u[0] - u_exact[0]),(REAL)2.0)+pow((REAL)fabs(u[1] - u_exact[1]),(REAL)2.0);

    //values[2] : erro estimado na norma H1
    //REAL SemiH1 =0.;
    //for(int i = 0; i < 2; i++) for(int j = 0; j < 2; j++) SemiH1 += (du(i,j) - du_exact(i,j)) * (du(i,j) - du_exact(i,j));
    //values[2] = values[1] + SemiH1;
}
