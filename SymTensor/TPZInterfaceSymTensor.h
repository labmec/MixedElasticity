#ifndef TPZINTERFACESYMTENSOR_H
#define TPZINTERFACESYMTENSOR_H

#include "TPZLagrangeMultiplierCS.h"

class TPZInterfaceSymTensor : public TPZLagrangeMultiplierCS<STATE> {
    
    typedef STATE TVar;
public:
    // Constructor
    TPZInterfaceSymTensor(int id, int dimension);

    // Destructor
    ~TPZInterfaceSymTensor() override;

    // Override any necessary methods from TPZMatInterfaceCombinedSpaces
    virtual void ContributeInterface(const TPZMaterialDataT<TVar> &data,
                                     const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                                     const std::map<int, TPZMaterialDataT<TVar>> &dataright, REAL weight,
                                     TPZFMatrix<TVar> &ek,
                                     TPZFMatrix<TVar> &ef) override;


};

#endif // TPZINTERFACESYMTENSOR_H
