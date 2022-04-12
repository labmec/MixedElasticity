#include "AiryFunctionLoadedBeam.h"


AiryFunctionLoadedBeam::AiryFunctionLoadedBeam(TPZVec<REAL> axis, REAL height, REAL length){
    
    gAxis.resize(3);

    gAxis[0] = axis[0];
    gAxis[1] = axis[1];
    gAxis[2] = 0.;
    gHeight = height;
    gLength = length;
    
}

void AiryFunctionLoadedBeam::SetElasticConstants(REAL elastic, REAL poisson, REAL inertia, REAL force){
    gE = elastic;
    gPoisson = poisson;
    gForce = force;
    gInertia = inertia;
}

void AiryFunctionLoadedBeam::GetDisplacement(TPZVec<REAL> &x, TPZFMatrix<REAL> &disp){
    
    REAL l = gLength / 2.;//metade do comprimento
    REAL c = gHeight / 2.;//metade da altura

    REAL deflection =  5.*gForce*l*l*l*l/(24.*gE*gInertia)*(1+12.*c*c/(5.*l*l)*(4./5.+gPoisson/2));
    disp(0,0) = gForce/(2.*gE*gInertia)*((l*l*x[0]-x[0]*x[0]*x[0]/3)*x[1] + x[0]*(2.*x[1]*x[1]*x[1]/3.-2.*c*c*x[1]) + gPoisson*x[0]*(x[1]*x[1]*x[1]/3.-c*c*x[1]+2.*c*c*c/3));
    disp(0,1) = -gForce/(2.*gE*gInertia)*(x[1]*x[1]*x[1]*x[1]/12. - c*c*x[1]*x[1]/2. + 2.*c*c*c*x[1]/3. + gPoisson*((l*l-x[0]*x[0])*x[1]*x[1]/2.+x[1]*x[1]*x[1]*x[1]/6.-c*c*x[1]*x[1]/5.)) 
            -gForce/(2.*gE*gInertia)*(l*l*x[0]*x[0]/2. - x[0]*x[0]*x[0]*x[0]/12. - c*c*x[0]*x[0]/5. + (1.+0.5*gPoisson)*c*c*x[0]*x[0]) + deflection;

  
}

void AiryFunctionLoadedBeam::GetStress(TPZVec<REAL> &x, TPZFMatrix<REAL> &stress, TPZFMatrix<REAL> &divStress){
    
    REAL l = gLength / 2.;//metade do comprimento
    REAL c = gHeight / 2.;//metade da altura

    //Sigma xx
    stress(0,0) = gForce / (2.*gInertia) * ((l*l-x[0]*x[0])*x[1] + 2.*x[1]*x[1]*x[1]/3. - 2.*c*c*x[1]/5.); 
    //Sigma yy
    stress(0,1) = -gForce / (2.*gInertia) * (x[1]*x[1]*x[1]/3. - c*c*x[1] + 2.*c*c*c/3.); 
    //Tau xy
    stress(0,2) = -gForce / (2.*gInertia) * (c*c - x[1]*x[1])*x[0]; 

    divStress.Zero();
}

int AiryFunctionLoadedBeam::ClassId() const{
    return Hash("AiryFunctionLoadedBeam") << 1;
}
