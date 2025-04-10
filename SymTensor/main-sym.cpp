#include <fstream>

#include "pzlog.h"
#include <iostream>
#include <string>

#include <cmath>
#include <set>


#ifdef PZ_LOG
static TPZLogger logger("testmhm");
#endif


int main(int argc, char *argv[]) {
    //    TPZMaterial::gBigNumber = 1.e16;
    
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    return 0;
}
