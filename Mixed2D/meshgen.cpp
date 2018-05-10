// 
#include "meshgen.h"

#include "pzgmesh.h"
#include "pzgengrid.h"
#include "pzgeoel.h"
#include "TPZRefPatternTools.h"
#include "pzcheckgeom.h"
#include "TPZVTKGeoMesh.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"

#include "TPZMaterial.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZSSpStructMatrix.h"

#ifndef USING_MKL
#include "pzskylstrmatrix.h"
#endif


#include "pzlog.h"

#ifdef _AUTODIFF
#include "fadType.h"

static FADFADREAL FADsin(FADFADREAL x) {
    FADREAL_ sinaval = sin(x.val());
    FADREAL_ cosaval = cos(x.val());
    FADFADREAL sina(x.size(), sinaval);
    for (int i = 0; i < x.size(); i++) {
        sina.fastAccessDx(i) = cosaval * x.dx(i);
    }
    return sina;
}

static FADFADREAL FADcos(FADFADREAL x) {
    FADREAL_ sinaval = sin(x.val());
    FADREAL_ cosaval = cos(x.val());
    FADFADREAL cosa(x.size(), cosaval);
    for (int i = 0; i < x.size(); i++) {
        cosa.fastAccessDx(i) = -sinaval * x.dx(i);
    }
    return cosa;
}

static FADFADREAL FADexp(FADFADREAL x) {
    FADREAL_ expaval = exp(x.val());
    FADFADREAL expa(x.size(), expaval);
    for (int i = 0; i < x.size(); i++) {
        expa.fastAccessDx(i) = expaval * x.dx(i);
    }
    return expa;
}

static FADFADREAL FADsqrt(FADFADREAL x) {
    FADREAL_ fadres = sqrt(x.val());
    FADFADREAL resa(x.size(), fadres);
    for (int i = 0; i < x.size(); i++) {
        resa.fastAccessDx(i) = REAL(0.5) / fadres * x.dx(i);
    }
    return resa;
}

static FADFADREAL FADatan(FADFADREAL x) {
    FADREAL_ fadres = atan(x.val());
    FADFADREAL resa(x.size(), fadres);
    for (int i = 0; i < x.size(); i++) {
        resa.fastAccessDx(i) = 1. / (1 + x.val() * x.val()) * x.dx(i);
    }
    return resa;
}

template<class TVar>
void TElasticityExample1::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) {
    switch (fProblemType) {
        case Etest1:
        {
            disp[0] = TVar(1. / 27.) * x[0] * x[0] * x[1] * x[1] * cos(TVar(6. * M_PI) * x[0]) * sin(TVar(7. * M_PI) * x[1]);
            disp[1] = TVar(0.2) * exp(x[1]) * sin(TVar(4. * M_PI) * x[0]);
            break;
        }
        case EThiago:
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            //disp[0] = ((1-x[0]*x[0])*(1+x[1]*x[1]*x[1]*x[1]));
            //disp[1] = ((1-x[1]*x[1])*(1+x[0]*x[0]*x[0]*x[0]));
            disp[0] = (TVar) cos(M_PI * x[0])*(TVar) sin(2 * M_PI * x[1]);
            disp[1] = (TVar) cos(M_PI * x[1])*(TVar) sin(M_PI * x[0]);
            break;
        }
        case ERot: //rotation
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] = (TVar) - x[1];
            disp[1] = (TVar) x[0];
            break;
        }
        case EShear://pure shear
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += (TVar) x[1];
            disp[1] += (TVar) 0.;
            break;
        }
        case EStretchx://strech x
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += (TVar) x[0] * x[0];
            disp[1] += (TVar) 0.;
            break;
        }
        case EStretchy://strech y
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += (TVar) 0.;
            disp[1] += (TVar) x[1] * x[1];
            break;
        }
        case EDispx:
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += 1.;
            disp[0] += 0.;
            break;
        }
        case EDispy:
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += (TVar) 0.;
            disp[0] += (TVar) 1.;
            break;
        }
        case EPoly:
        {
            disp[0] = 0. * x[0]; //*x[0]*x[1];
            disp[1] = x[1] * x[0]; //(x[1]-1.)*(x[1]-x[0])*(x[0]+4.*x[1]);
            break;
        }
        default:
        {
            DebugStop();
        }
    }
}

template<>
void TElasticityExample1::uxy(const TPZVec<FADFADREAL > &x, TPZVec<FADFADREAL > &disp) {
    switch (fProblemType) {
        case Etest1:
        {
            FADFADREAL tmp = (FADFADREAL) (1. / 27.) * x[0] * x[0] * x[1] * x[1];
            disp[0] = tmp;
            FADFADREAL a = FADcos((FADFADREAL) (6. * M_PI) * x[0]);
            disp[0] = disp[0] * a;
            disp[0] *= FADsin((FADFADREAL) (7. * M_PI) * x[1]);
            disp[1] = (FADFADREAL) (0.2) * FADexp(x[1]) * FADsin((FADFADREAL) (4. * M_PI) * x[0]);
            break;
        }
        case EThiago:
        {
            disp[0] = FADcos(M_PI * x[0]) * FADsin(2 * M_PI * x[1]);
            disp[1] = FADcos(M_PI * x[1]) * FADsin(M_PI * x[0]);
            break;
        }
        case ERot://rotation
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] = (FADFADREAL) - x[1];
            disp[1] = (FADFADREAL) x[0];
            break;

        }
        case EShear://pure shear
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += (FADFADREAL) x[1];
            disp[1] += (FADFADREAL) 0.;
            break;
        }
        case EStretchx://strech x
        {
            disp[0] += (FADFADREAL) x[0] * x[0];
            disp[1] += (FADFADREAL) x[1]*0.;
            break;
        }
        case EStretchy://strech y
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += (FADFADREAL) x[0]*0.;
            disp[1] += (FADFADREAL) x[1] * x[1];
            break;
        }
        case EDispx:
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += 1.;
            disp[0] += 0.;
            break;
        }
        case EDispy:
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += (FADFADREAL) 0.;
            disp[0] += (FADFADREAL) 1.;
            break;
        }
        case EPoly:
        {
            disp[0] = 0. * x[0]; //*x[0]*x[1];
            disp[1] = x[1] * x[0]; //(x[1]-1.)*(x[1]-x[0])*(x[0]+4.*x[1]);
            break;
        }
        default:
        {
            DebugStop();
        }
    }
}

REAL TElasticityExample1::fElast = 1.;

REAL TElasticityExample1::fNu = 0.;

template<class TVar>
void TElasticityExample1::Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu) {
    //    Elast = (TVar(100.) * (TVar(1.) + TVar(0.3) * sin(TVar(10 * M_PI) * (x[0] - TVar(0.5))) * cos(TVar(10. * M_PI) * x[1])));
    Elast = TVar(fElast); //TVar(206.815026);
    nu = TVar(fNu); //TVar(0.30400395);
}

template<>
void TElasticityExample1::Elastic(const TPZVec<double> &x, double &Elast, double &nu) {
    //  Elast = 1000.;
    //    Elast = (100. * (1. + 0.3 * sin(10 * M_PI * (x[0] - 0.5)) * cos(10. * M_PI * x[1])));
    Elast = fElast; //206.815026;
    nu = fNu; //0.30400395;
}

void TElasticityExample1::ElasticDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv) {
    TPZManVector<STATE> xstate(x.size());
    for (int i = 0; i < xstate.size(); i++) {
        xstate[i] = x[i];
    }
    STATE E, nu;
    Elastic(xstate, E, nu);
    result[0] = E;
    result[1] = nu;
}

TPZAutoPointer<TPZFunction<STATE> > TElasticityExample1::ConstitutiveLawFunction() {
    TPZAutoPointer<TPZFunction<STATE> > result;
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(TElasticityExample1::ElasticDummy);
    dummy->SetPolynomialOrder(1);
    result = TPZAutoPointer<TPZFunction<STATE> >(dummy);
    return result;
}

template<class TVar>
void TElasticityExample1::graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &grad) {
    TPZManVector<Fad<TVar>, 3> xfad(x.size());
    for (int i = 0; i < 2; i++) {

        Fad<TVar> temp = Fad<TVar>(2, i, x[i]);
        xfad[i] = temp;

    }
    xfad[2] = x[2];
    TPZManVector<Fad<TVar>, 3> result(2);
    uxy(xfad, result);
    grad.Resize(2, 2);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            grad(i, j) = result[i].d(j);
        }
    }
}

void TElasticityExample1::GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) {
    TPZManVector<Fad<REAL>, 3> xfad(x.size());
    for (int i = 0; i < 2; i++) {
        Fad<REAL> temp = Fad<REAL>(2, i, x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<REAL>, 3> result(2);
    uxy(xfad, result);
    //std::cout<<"result  "<<result[0].d(0)<<endl;
    gradu.Redim(2, 2);
    u[0] = result[0].val();
    u[1] = result[1].val();
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            gradu(i, j) = result[j].d(i);
        }
    }

}

template<>
void TElasticityExample1::graduxy(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &grad) {
    TPZManVector<Fad<Fad<REAL> >, 3> xfad(x.size());
    for (int i = 0; i < 3; i++) {
        Fad<Fad<REAL> > temp = Fad<Fad<REAL> >(3, Fad<REAL>(3, 0.));
        //      Fad<REAL> temp = Fad<REAL>(2,i,x[i]);
        temp.val() = x[i];
        Fad<REAL> temp3(3, 0.);
        for (int j = 0; j < 3; j++) {
            temp.fastAccessDx(j) = temp3;
        }
        Fad<REAL> temp2(3, 1.);
        temp.fastAccessDx(i) = temp2;
        //      Fad<REAL> temp = Fad<REAL>(2,i,x[i]);
        xfad[i] = temp;
        //      xfad[i] = temp;
    }
    TPZManVector<Fad<Fad<REAL> >, 3> result(2);
    uxy(xfad, result);
    grad.Resize(2, 2);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            grad(i, j) = result[i].d(j);
        }
    }
}

template<class TVar>
void TElasticityExample1::Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) {
    TPZFNMatrix<4, TVar> grad;
    TVar E, nu;
    Elastic(x, E, nu);
    TVar Fac = E / ((TVar) 1. + nu) / ((TVar(1.) - TVar(2.) * nu));
    graduxy(x, grad);
    TPZManVector<TVar, 2> result(2, 0.);
    uxy(x, result);
    sigma.Resize(2, 2);
    if (fStressState == EPlaneStress) {
        sigma(0, 0) = Fac * ((TVar(1.) - nu) * grad(0, 0) + nu * grad(1, 1));
        sigma(1, 1) = Fac * ((TVar(1.) - nu) * grad(1, 1) + nu * grad(0, 0));
        sigma(0, 1) = E / (TVar(2.)*(TVar(1.) + nu))*(grad(0, 1) + grad(1, 0));
        sigma(1, 0) = sigma(0, 1);
    } else if (fStressState == EPlaneStrain) {
        sigma(0, 0) = (E / ((1. + nu)*(1. - 2. * nu)))*((1. - nu) * grad(0, 0) + nu * grad(1, 1));
        sigma(1, 1) = (E / ((1. + nu)*(1. - 2. * nu)))*((1. - nu) * grad(1, 1) + nu * grad(0, 0));
        sigma(0, 1) = E / (TVar(2.)*(TVar(1.) + nu))*(grad(0, 1) + grad(1, 0));
        sigma(1, 0) = sigma(0, 1);
    } else if (fStressState == EAxiSymmetric) {
        sigma(0, 0) = (E / ((1. + nu)*(1. - 2. * nu)))*((1. - nu) * grad(0, 0) + nu * grad(1, 1) + nu * result[0] / x[0]);
        sigma(1, 1) = (E / ((1. + nu)*(1. - 2. * nu)))*((1. - nu) * grad(1, 1) + nu * grad(0, 0) + nu * result[0] / x[0]);
        sigma(0, 1) = E / (TVar(2.)*(TVar(1.) + nu))*(grad(0, 1) + grad(1, 0));
        sigma(1, 0) = sigma(0, 1);
    } else {
        DebugStop();
    }

}

template<>
void TElasticityExample1::Sigma(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &sigma) {
    typedef Fad<REAL> TVar;
    TPZFNMatrix<4, Fad<REAL> > grad;
    Fad<REAL> E, nu;
    Elastic(x, E, nu);
    Fad<REAL> Fac = E / (Fad<REAL>(1.) + nu) / ((Fad<REAL>(1.) - Fad<REAL>(2.) * nu));
    graduxy(x, grad);
    sigma.Resize(2, 2);
    TPZVec<Fad<REAL> > result(2, 0.);
    uxy(x, result);
    if (fStressState == EPlaneStress) {
        sigma(0, 0) = Fac * ((Fad<REAL>(1.) - nu) * grad(0, 0) + nu * grad(1, 1));
        sigma(1, 1) = Fac * ((Fad<REAL>(1.) - nu) * grad(1, 1) + nu * grad(0, 0));
        sigma(0, 1) = E / (Fad<REAL>(2.)*(Fad<REAL>(1.) + nu))*(grad(0, 1) + grad(1, 0));
        sigma(1, 0) = sigma(0, 1);
    } else if (fStressState == EPlaneStrain) {
        sigma(0, 0) = (E / ((1. + nu)*(1. - 2. * nu)))*((1. - nu) * grad(0, 0) + nu * grad(1, 1));
        sigma(1, 1) = (E / ((1. + nu)*(1. - 2. * nu)))*((1. - nu) * grad(1, 1) + nu * grad(0, 0));
        sigma(0, 1) = E / (Fad<REAL>(2.)*(TVar(1.) + nu))*(grad(0, 1) + grad(1, 0));
        sigma(1, 0) = sigma(0, 1);
    } else if (fStressState == EAxiSymmetric) {
        sigma(0, 0) = Fac * ((Fad<REAL>(1.) - nu) * grad(0, 0) + nu * grad(1, 1) + nu * (result[0] / x[0]));
        sigma(1, 1) = Fac * ((Fad<REAL>(1.) - nu) * grad(1, 1) + nu * grad(0, 0) + nu * (result[0] / x[0]));
        sigma(0, 1) = E / (Fad<REAL>(2.)*(Fad<REAL>(1.) + nu))*(grad(0, 1) + grad(1, 0));
        sigma(1, 0) = sigma(0, 1);
    } else {
        DebugStop();
    }


}

template
void TElasticityExample1::Sigma<REAL>(const TPZVec<REAL> &x, TPZFMatrix<REAL> &divsigma);

template<class TVar>
void TElasticityExample1::DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma) {
    TPZVec<TVar> result(2, 0.);
    uxy(x, result);
    TVar E, nu;
    Elastic(x, E, nu);
    TPZManVector<Fad<TVar>, 3> xfad(x.size());
    for (int i = 0; i < 3; i++) {
        xfad[i] = Fad<TVar>(3, i, x[i]);
    }
    TPZFNMatrix<4, Fad<TVar> > sigma(2, 2);
    Sigma(xfad, sigma);
    TPZFNMatrix<4, TVar> sigma1(2, 2);
    Sigma(x, sigma1);
    if (fStressState != EAxiSymmetric) {
        divsigma[0] = sigma(0, 0).dx(0) + sigma(0, 1).dx(1);
        divsigma[1] = sigma(1, 0).dx(0) + sigma(1, 1).dx(1);
    } else {
        divsigma[0] = sigma(0, 0).dx(0) + sigma(0, 1).dx(1) + sigma(0, 0).val() / x[0] - result[0] / (E * x[0] * x[0]);
        divsigma[1] = sigma(1, 0).dx(0) + sigma(1, 1).dx(1) + sigma(1, 0).val() / x[0];
    }

}

TPZAutoPointer<TPZFunction<STATE> > TElasticityExample1::ForcingFunction() {
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Force);
    dummy->SetPolynomialOrder(30);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;
}

TPZAutoPointer<TPZFunction<STATE> > TElasticityExample1::ValueFunction() {
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(GradU);
    dummy->SetPolynomialOrder(10);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;

}


template
void TElasticityExample1::DivSigma<REAL>(const TPZVec<REAL> &x, TPZVec<REAL> &divsigma);
template
void TElasticityExample1::Sigma<Fad<REAL> >(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &sigma);


TElasticityExample1::EDefState TElasticityExample1::fProblemType = TElasticityExample1::EDispx;

TElasticityExample1::EStressState TElasticityExample1::fStressState = TElasticityExample1::EPlaneStrain;

template<class TVar>
void TElasticityExample2::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) {
    if (fProblemType == EDispx)//rotation
    {
        disp[0] = (TVar) - x[1];
        disp[1] = (TVar) x[0];
    }

}

template<class TVar>
void TLaplaceExample1::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) {
    disp[0] = sin((TVar) M_PI * x[0]) * sin((TVar) M_PI * x[1]);
    TVar r = sqrt(x[0] * x[0] + x[1] * x[1]);
    TVar atanco = (r - (TVar) 0.5)*100.;
    TVar freq = 10.;
    TVar mult = (TVar(1) + TVar(0.3) * sin(TVar(M_PI) * x[0] * freq))*(TVar(1) + TVar(0.5) * cos(TVar(M_PI) * r * freq));
    disp[0] = atan(atanco) * mult;
}

template<>
void TLaplaceExample1::uxy(const TPZVec<FADFADREAL > &x, TPZVec<FADFADREAL > &disp) {
    //  disp[0] = FADsin((FADFADREAL)(M_PI)*x[0])*FADsin((FADFADREAL)(M_PI)*x[1]);
    FADFADREAL r = FADsqrt(x[0] * x[0] + x[1] * x[1]);
    FADFADREAL atanco = (r - (FADFADREAL) 0.5)*100.;
    FADFADREAL freq = (FADFADREAL) 10.;
    FADFADREAL mult = ((FADFADREAL) 1. + (FADFADREAL) 0.3 * FADsin((FADFADREAL) M_PI * x[0] * freq))*((FADFADREAL) 1. + (FADFADREAL) 0.5 * FADcos((FADFADREAL) M_PI * r * freq));
    disp[0] = FADatan(atanco) * mult;

}

template<class TVar>
void TLaplaceExample1::Permeability(const TPZVec<TVar> &x, TVar &Perm) {
    Perm = (TVar) (1.);
}

void TLaplaceExample1::PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv) {
    TPZManVector<STATE, 3> xloc(x.size());
    for (auto i : xloc) {
        xloc[i] = x[i];
    }
    STATE Perm;
    Permeability(xloc, Perm);
    deriv.Zero();
    deriv(0, 0) = Perm;
    deriv(1, 1) = Perm;
    deriv(2, 0) = 1. / Perm;
    deriv(3, 1) = 1. / Perm;
}

template<class TVar>
void TLaplaceExample1::graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad) {
    TPZManVector<Fad<TVar>, 3> xfad(x.size());
    for (int i = 0; i < 2; i++) {
        Fad<TVar> temp = Fad<TVar>(2, i, x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<TVar>, 3> result(1);
    uxy(xfad, result);
    grad.resize(2);
    for (int i = 0; i < 2; i++) {
        grad[i] = result[0].d(i);
    }
}

void TLaplaceExample1::GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) {
    TPZManVector<Fad<REAL>, 3> xfad(x.size());
    for (int i = 0; i < 2; i++) {
        Fad<REAL> temp = Fad<REAL>(2, i, x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<REAL>, 3> result(2);
    uxy(xfad, result);
    gradu.Redim(2, 1);
    u[0] = result[0].val();
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            gradu(i, j) = result[j].d(i);
        }
    }

}

template<class TVar>
void TLaplaceExample1::Sigma(const TPZVec<TVar> &x, TPZVec<TVar> &sigma) {
    TPZManVector<TVar, 3> grad;
    TVar Perm;
    Permeability(x, Perm);
    graduxy(x, grad);
    sigma.resize(2);
    sigma[0] = -Perm * grad[0];
    sigma[1] = -Perm * grad[1];

}

template<class TVar>
void TLaplaceExample1::DivSigma(const TPZVec<TVar> &x, TVar &divsigma) {
    TPZManVector<Fad<TVar>, 3> xfad(x.size());
    for (int i = 0; i < 2; i++) {
        xfad[i] = Fad<TVar>(2, i, x[i]);
    }
    TPZManVector<Fad<TVar>, 3> sigma(2);
    Sigma(xfad, sigma);
    divsigma = sigma[0].dx(0) + sigma[1].dx(1);

}

TPZAutoPointer<TPZFunction<STATE> > TLaplaceExample1::ForcingFunction() {
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Force);
    dummy->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;
}

TPZAutoPointer<TPZFunction<STATE> > TLaplaceExample1::ValueFunction() {
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(GradU);
    dummy->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;

}

TPZAutoPointer<TPZFunction<STATE> > TLaplaceExample1::ConstitutiveLawFunction() {
    TPZAutoPointer<TPZFunction<STATE> > result;
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(PermeabilityDummy);
    dummy->SetPolynomialOrder(4);
    result = TPZAutoPointer<TPZFunction<STATE> >(dummy);
    return result;

}

//ExactFunc *Exact();

template<class TVar>
void TLaplaceExampleSmooth::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) {
    disp[0] = x[0];
}

template<>
void TLaplaceExampleSmooth::uxy(const TPZVec<FADFADREAL > &x, TPZVec<FADFADREAL > &disp) {
    disp[0] = x[0];

}

template<class TVar>
void TLaplaceExampleSmooth::Permeability(const TPZVec<TVar> &x, TVar &Perm) {
    Perm = (TVar) (1.);
}

void TLaplaceExampleSmooth::PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv) {
    TPZManVector<STATE, 3> xloc(x.size());
    for (auto i : xloc) {
        xloc[i] = x[i];
    }
    STATE Perm;
    Permeability(xloc, Perm);
    deriv.Zero();
    deriv(0, 0) = Perm;
    deriv(1, 1) = Perm;
    deriv(2, 0) = 1. / Perm;
    deriv(3, 1) = 1. / Perm;
}

template<class TVar>
void TLaplaceExampleSmooth::graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad) {
    TPZManVector<Fad<TVar>, 3> xfad(x.size());
    for (int i = 0; i < 2; i++) {
        Fad<TVar> temp = Fad<TVar>(2, i, x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<TVar>, 3> result(1);
    uxy(xfad, result);
    grad.resize(2);
    for (int i = 0; i < 2; i++) {
        grad[i] = result[0].d(i);
    }
}

void TLaplaceExampleSmooth::GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) {
    TPZManVector<Fad<REAL>, 3> xfad(x.size());
    for (int i = 0; i < 2; i++) {
        Fad<REAL> temp = Fad<REAL>(2, i, x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<REAL>, 3> result(2);
    uxy(xfad, result);
    gradu.Redim(2, 1);
    u[0] = result[0].val();
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            gradu(i, j) = result[j].d(i);
        }
    }

}

template<class TVar>
void TLaplaceExampleSmooth::Sigma(const TPZVec<TVar> &x, TPZVec<TVar> &sigma) {
    TPZManVector<TVar, 3> grad;
    TVar Perm;
    Permeability(x, Perm);
    graduxy(x, grad);
    sigma.resize(2);
    sigma[0] = -Perm * grad[0];
    sigma[1] = -Perm * grad[1];

}

template<class TVar>
void TLaplaceExampleSmooth::DivSigma(const TPZVec<TVar> &x, TVar &divsigma) {
    TPZManVector<Fad<TVar>, 3> xfad(x.size());
    for (int i = 0; i < 2; i++) {
        xfad[i] = Fad<TVar>(2, i, x[i]);
    }
    TPZManVector<Fad<TVar>, 3> sigma(2);
    Sigma(xfad, sigma);
    divsigma = sigma[0].dx(0) + sigma[1].dx(1);

}

TPZAutoPointer<TPZFunction<STATE> > TLaplaceExampleSmooth::ForcingFunction() {
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Force);
    dummy->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;
}

TPZAutoPointer<TPZFunction<STATE> > TLaplaceExampleSmooth::ValueFunction() {
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(GradU);
    dummy->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;

}

TPZAutoPointer<TPZFunction<STATE> > TLaplaceExampleSmooth::ConstitutiveLawFunction() {
    TPZAutoPointer<TPZFunction<STATE> > result;
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(PermeabilityDummy);
    dummy->SetPolynomialOrder(4);
    result = TPZAutoPointer<TPZFunction<STATE> >(dummy);
    return result;

}

#endif
