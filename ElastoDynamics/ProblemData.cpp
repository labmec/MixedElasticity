#include <iostream>
#include <string>
#include <fstream>
#include <pzerror.h>
#include "ProblemData.h"
#include "dirs_config.h"

using namespace std;

// constructor
ProblemData::ProblemData() : AnalyticSolution("none"), MeshName("none"), StressOrder(0), 
                             DispOrder(0), Dimension(0), VTKRes(0), ShouldCondense(false),
                             FiveSpaces(false), PlaneStress(0), InterfaceID(20), LambdaID(10) {
    
    Boundaries.clear();
    Boundaries.reserve(10);
    Domains.clear();
    Domains.reserve(10);
    MeshVector.resize(3);
}

// deconstructor
ProblemData::~ProblemData() {
}

// readjson function. takes a json function as parameter and completes the required simulation data
void ProblemData::ReadJson(std::string file) {

    std::string path = std::string(INPUTDIR) + "/" + file;
    std::ifstream filejson(path);
    json input = json::parse(filejson, nullptr, true, true); // to ignore comments in json file

    if (input.find("MeshName") == input.end()) DebugStop();
    MeshName = input["MeshName"];

    if (input.find("Dimension") == input.end()) DebugStop();
    Dimension = input["Dimension"];

    if (input.find("StressOrder") == input.end()) DebugStop();
    StressOrder = input["StressOrder"];

    if (input.find("DispOrder") == input.end()) DebugStop();
    DispOrder = input["DispOrder"];

    if (input.find("PlaneStress") == input.end()) DebugStop();
    PlaneStress = input["PlaneStress"];

    if (input.find("VTKRes") == input.end()) DebugStop();
    VTKRes = input["VTKRes"];

    if (input.find("ShouldCondense") == input.end()) DebugStop();
    ShouldCondense = input["ShouldCondense"];

    if (input.find("FiveSpaces") == input.end()) DebugStop();
    FiveSpaces = input["FiveSpaces"];

    if (input.find("AnalyticSolution") == input.end()) DebugStop();
    AnalyticSolution = input["AnalyticSolution"];

    DomainData domaindata;
    for (auto &domainjson : input["Domains"]) {
        if (domainjson.find("name") == domainjson.end()) DebugStop();
        domaindata.name = domainjson["name"];
        if (domainjson.find("matID") == domainjson.end()) DebugStop();
        domaindata.matID = domainjson["matID"];
        if (domainjson.find("E") == domainjson.end()) DebugStop();
        domaindata.E = domainjson["E"];
        if (domainjson.find("nu") == domainjson.end()) DebugStop();
        domaindata.nu = domainjson["nu"];

        Domains.push_back(domaindata);
    }

    BcData bcdata;
    for (auto &bcjson : input["Boundaries"]) {
        if (bcjson.find("name") == bcjson.end()) DebugStop();
        bcdata.name = bcjson["name"];
        if (bcjson.find("type") == bcjson.end()) DebugStop();
        bcdata.type = bcjson["type"];
        if (bcjson.find("value") == bcjson.end()) DebugStop();
        for (int i = 0; i < Dimension; i++)
        {
            bcdata.value[i] = bcjson["value"][i];
        }
        if (bcjson.find("matID") == bcjson.end()) DebugStop();
        bcdata.matID = bcjson["matID"];

        bcdata.domainID = 1;
        if (bcjson.find("domainID") != bcjson.end())
            bcdata.domainID = bcjson["domainID"];

        Boundaries.push_back(bcdata);
    }

    if (input.find("InterfaceID") != input.end()) {
        InterfaceID = input["InterfaceID"];
    }
    if (input.find("LambdaID") != input.end()) {
        LambdaID = input["LambdaID"];
    }

    if (FiveSpaces) MeshVector.resize(5);
}