/**
 * @file
 * @brief 
 */

#ifndef AIRYFUNCTION_LOADEDBEAM
#define AIRYFUNCTION_LOADEDBEAM

#include <iostream>
#include "TPZAnalyticSolution.h"
#include "pzvec.h"

class AiryFunctionLoadedBeam{

public:
    //Default creator
    AiryFunctionLoadedBeam() = default;

    //Creates an Airy stress function analytic solution environment for a given axis point and beam height and lenght
    AiryFunctionLoadedBeam(TPZVec<REAL> axis, REAL height, REAL length);

    void SetElasticConstants(REAL elastic, REAL poisson, REAL inertia, REAL force = 1.);

    void GetDisplacement(TPZVec<REAL> &x, TPZFMatrix<REAL> &disp);

    void GetStress(TPZVec<REAL> &x, TPZFMatrix<REAL> &stress, TPZFMatrix<REAL> &divStress);

    int NShapeF(){return 1;};

    int ClassId() const;

private:

    //Elasticity
    REAL gE;

    //Poisson ration
    REAL gPoisson;

    //Aplied force (suposed to be unitary?)
    REAL gForce;

    //Moment of Inertia
    REAL gInertia;

    //Beam axis coordinates
    TPZVec<REAL> gAxis;

    //Beam height
    REAL gHeight;

    //Beam length
    REAL gLength;

};

#endif