
#ifndef MESHGEN
#define MESHGEN

#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzfunction.h"

#include <string>

class TPZGeoMesh;

TPZGeoMesh *MalhaGeomFredQuadrada(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, TPZVec<int64_t> &coarseindices, int ndiv);

struct TAnalyticSolution;

typedef void (ExactFunc) (const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);

struct TAnalyticSolution {
    virtual TPZAutoPointer<TPZFunction<STATE> > ForcingFunction() = 0;

    virtual TPZAutoPointer<TPZFunction<STATE> > ValueFunction() = 0;

    virtual TPZAutoPointer<TPZFunction<STATE> > ConstitutiveLawFunction() = 0;

    virtual ExactFunc *Exact() = 0;

    virtual ~TAnalyticSolution() {

    }
};

#ifdef _AUTODIFF

struct TElasticityExample1 : public TAnalyticSolution {

    enum EDefState {
        ENone, EDispx, EDispy, ERot, EStretchx, EStretchy, EShear, Etest1, EThiago, EPoly
    };

    static EDefState fProblemType;

    enum EStressState {
        ENone2, EPlaneStress, EPlaneStrain, EAxiSymmetric
    };

    static REAL fElast;

    static REAL fNu;

    static EStressState fStressState;

    static void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) {
        TPZManVector<REAL, 3> locforce(2);
        DivSigma(x, locforce);
        force[0] = -locforce[0];
        force[1] = -locforce[1];
    }

    virtual TPZAutoPointer<TPZFunction<STATE> > ForcingFunction();

    virtual TPZAutoPointer<TPZFunction<STATE> > ValueFunction();

    virtual TPZAutoPointer<TPZFunction<STATE> > ConstitutiveLawFunction();

    static void GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);

    virtual ExactFunc *Exact() {
        return GradU;
    }

    virtual ~TElasticityExample1() {

    }

    template<class TVar>
    static void Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma);

    template<class TVar>
    static void DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma);

    template<class TVar>
    static void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);

    static void Dirichlet(const TPZVec<REAL> &x, TPZVec<STATE> &disp) {
        TPZManVector<REAL, 3> disploc(2, 0.);
        uxy(x, disploc);
        disp[0] = disploc[0];
        disp[1] = disploc[1];
    }

    template<class TVar>
    static void graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &grad);

    template<class TVar>
    static void Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu);

    static void ElasticDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);

};

struct TElasticityExample2 : public TElasticityExample1 {
    template<class TVar>
    static void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);
};

struct TLaplaceExample : public TAnalyticSolution {
    virtual TPZAutoPointer<TPZFunction<STATE> > ForcingFunction();

    virtual TPZAutoPointer<TPZFunction<STATE> > ValueFunction();

    virtual TPZAutoPointer<TPZFunction<STATE> > ConstitutiveLawFunction();

    virtual ~TLaplaceExample() {

    }

    static void GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);

    virtual ExactFunc *Exact() {
        return GradU;
    }

    template<class TVar>
    static void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);

    template<class TVar>
    static void graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad);

    template<class TVar>
    static void Permeability(const TPZVec<TVar> &x, TVar &Elast);

    static void PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);

    template<class TVar>
    static void Sigma(const TPZVec<TVar> &x, TPZVec<TVar> &sigma);

    template<class TVar>
    static void DivSigma(const TPZVec<TVar> &x, TVar &divsigma);

    static void Dirichlet(const TPZVec<REAL> &x, TPZVec<STATE> &disp) {
        TPZManVector<REAL, 3> disploc(2, 0.);
        uxy(x, disploc);
        for (int i = 0; i < 1; i++) disp[i] = disploc[i];
    }

    static void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) {
        REAL locforce;
        DivSigma(x, locforce);
        force[0] = locforce;
    }

};

struct TLaplaceExampleSmooth : public TAnalyticSolution {
    virtual TPZAutoPointer<TPZFunction<STATE> > ForcingFunction();

    virtual TPZAutoPointer<TPZFunction<STATE> > ValueFunction();

    virtual TPZAutoPointer<TPZFunction<STATE> > ConstitutiveLawFunction();

    virtual ~TLaplaceExampleSmooth() {

    }

    static void GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);

    virtual ExactFunc *Exact() {
        return GradU;
    }

    template<class TVar>
    static void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);

    template<class TVar>
    static void graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad);

    template<class TVar>
    static void Permeability(const TPZVec<TVar> &x, TVar &Elast);

    static void PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);

    template<class TVar>
    static void Sigma(const TPZVec<TVar> &x, TPZVec<TVar> &sigma);

    template<class TVar>
    static void DivSigma(const TPZVec<TVar> &x, TVar &divsigma);

    static void Dirichlet(const TPZVec<REAL> &x, TPZVec<STATE> &disp) {
        TPZManVector<REAL, 3> disploc(2, 0.);
        uxy(x, disploc);
        for (int i = 0; i < 1; i++) disp[i] = disploc[i];
    }

    static void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) {
        REAL locforce;
        DivSigma(x, locforce);
        force[0] = locforce;
    }

};

#endif

#endif
