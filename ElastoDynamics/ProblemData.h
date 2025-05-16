#pragma once

#include <iostream>
#include "json.hpp"
#include <pzfmatrix.h>
#include <pzvec.h>
#include <string>

// declaration of simulation data class.
// all the data herein used are storaged in a .json file. It can be called and storaged using ReadJson

class ProblemData
{
    // struct responsible to summarize all the data from every domain
    struct DomainData {
        std::string name = "none"; // domains name
        int matID = -1; // domain material ID
        REAL E = -1.; // domain young modulus
        REAL nu = -1.; // poisson
    };
    
    // struct responsible to store boundary condition data
    struct BcData {
        std::string name = "none"; // name of the bc
        int type = 0; // bc type (explained below)
        TPZManVector<REAL,3>  value = {0., 0., 0.}; // bc value
        int matID = 0; // bc material ID
        int domainID = 1; // domain material ID whose this BC is associated
    };

public:
    using json = nlohmann::json; // declaration of json class
    
    std::string MeshName;
    
    int Dimension;

    int PlaneStress; // 1 for plane stress, 0 for plane strain

    std::string AnalyticSolution;
    
    int StressOrder; // polynomial approximation order for velocity
    
    int DispOrder; // polynomial approximation order for traction
    
    int VTKRes;

    bool ShouldCondense;

    bool FiveSpaces;
    
    std::vector<DomainData> Domains; // vector containing every domain created

    std::vector<BcData> Boundaries; // vector containg all the velocity bcs info
    
    int InterfaceID;
    
    int LambdaID;
    
    TPZManVector<TPZCompMesh*,7> MeshVector;
    
    ProblemData();
    
    ~ProblemData();
    
    void ReadJson(std::string jsonfile);
};