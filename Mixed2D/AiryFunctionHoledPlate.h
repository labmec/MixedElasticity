/**
 * @file
 * @brief 
 */

#ifndef AIRYFUNCTION_HOLEDPLATE
#define AIRYFUNCTION_HOLEDPLATE

#include <iostream>
#include "TPZAnalyticSolution.h"
#include "pzvec.h"

class AiryFunctionHoledPlate{

public:
    //Default creator
    AiryFunctionHoledPlate() = default;

    //Creates an Airy stress function analytic solution environment for a given center point and hole radius
    AiryFunctionHoledPlate(TPZFMatrix<REAL> center, TPZVec<REAL> radius);

    void SetElasticConstants(REAL elastic, REAL poisson, REAL force = 1.);

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

    //Hole center coordinates
    TPZFMatrix<REAL> gCenter;

    //Hole radius
    TPZVec<REAL> gRadius;

    int nShape;

};

#endif