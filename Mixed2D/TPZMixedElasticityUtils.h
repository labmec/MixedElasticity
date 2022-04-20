#ifndef MIXED_ELASTICITY_UTILS
#define MIXED_ELASTICITY_UTILS

#include <iostream>
#include <set>
#include "TPZLinearAnalysis.h"
#include <pzgmesh.h> 

class TPZMixedElasticityUtils{

public:

    void SolveProblem(TPZLinearAnalysis &an, int fDomainMatID);

    void PrintResults(TPZLinearAnalysis &an);

};

#endif
