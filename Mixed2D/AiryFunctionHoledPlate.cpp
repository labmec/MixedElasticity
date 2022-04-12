#include "AiryFunctionHoledPlate.h"


AiryFunctionHoledPlate::AiryFunctionHoledPlate(TPZFMatrix<REAL> center, TPZVec<REAL> radius){
    auto nHoles = radius.size();
    gCenter.Resize(nHoles,3);
    gCenter.Zero();
    gRadius.resize(nHoles);

    for (int i = 0; i < nHoles; i++)
    {
        gCenter(i,0) = center(i,0);
        gCenter(i,1) = center(i,1);
        gRadius[i] = radius[i];
    }

}

void AiryFunctionHoledPlate::SetElasticConstants(REAL elastic, REAL poisson, REAL force){
    gE = elastic;
    gPoisson = poisson;
    gForce = force;
}

void AiryFunctionHoledPlate::GetDisplacement(TPZVec<REAL> &x, TPZFMatrix<REAL> &disp){
    
    //Lame constants
    REAL mu = gE/(2.*(1.+gPoisson));
    REAL kappa = (3-gPoisson)/(1.+gPoisson);

    for (int iHole = 0; iHole < gRadius.size(); iHole++)
    {
        //The components are computed in polar coordinates
        REAL theta = atan2(x[1]-gCenter(iHole,1), x[0]-gCenter(iHole,0));
        REAL r = sqrt((x[0]-gCenter(iHole,0))*(x[0]-gCenter(iHole,0)) + (x[1]-gCenter(iHole,1))*(x[1]-gCenter(iHole,1)));    
        REAL rad = gRadius[iHole];

        //Displacement
        disp(iHole,0) = gForce*rad/(8*mu)*(r/rad*(kappa+1)*cos(theta)+2.*rad/r*((1+kappa)*cos(theta)+cos(3.*theta))
                      - 2.*rad*rad*rad*cos(3.*theta)/(r*r*r)) + gCenter(iHole,0);
        disp(iHole,1) = gForce*rad/(8*mu)*(r/rad*(kappa-3)*sin(theta)+2.*rad/r*((1-kappa)*sin(theta)+sin(3.*theta))
                      - 2.*rad*rad*rad*sin(3.*theta)/(r*r*r)) + gCenter(iHole,1);
    }
}

void AiryFunctionHoledPlate::GetStress(TPZVec<REAL> &x, TPZFMatrix<REAL> &stress, TPZFMatrix<REAL> &divStress){
    
    for (int iHole = 0; iHole < gRadius.size(); iHole++)
    {
        //The components are computed in polar coordinates
        REAL theta = atan2(x[1]-gCenter(iHole,1), x[0]-gCenter(iHole,0));
        REAL r = sqrt((x[0]-gCenter(iHole,0))*(x[0]-gCenter(iHole,0)) + (x[1]-gCenter(iHole,1))*(x[1]-gCenter(iHole,1)));    
        REAL rad = gRadius[iHole];

        // Sigma xx
        stress(iHole,0) = gForce - gForce * rad*rad/(r*r) * (1.5*cos(2.*theta)+cos(4.*theta)) + gForce * 1.5 * rad*rad*rad*rad/(r*r*r*r)*cos(4.*theta);
        // Sigma yy
        stress(iHole,1) = -gForce * rad*rad/(r*r) * (0.5*cos(2.*theta)-cos(4.*theta)) - gForce * 1.5 * rad*rad*rad*rad/(r*r*r*r)*cos(4.*theta);
        //Tau xy
        stress(iHole,2) = -gForce * rad*rad/(r*r) * (0.5*sin(2.*theta)+sin(4.*theta)) + gForce * 1.5 * rad*rad*rad*rad/(r*r*r*r)*sin(4.*theta);

        
    }
    
    divStress.Zero();

}

int AiryFunctionHoledPlate::ClassId() const{
    return Hash("AiryFunctionHoledPlate") << 1;
}