#include <fstream>

#include "pzlog.h"
#include <iostream>
#include <string>

#include <cmath>
#include <set>

#ifdef PZ_LOG
static TPZLogger logger("pz.elasticity");
#endif

int main(int argc, char *argv[]) {


    TPZLogger::InitializePZLOG();
    std::cout << "FINISHED!" << std::endl;

    return 0;
}
