

#pragma once

#include "pzcmesh.h"
#include "pzgmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMHMApproxCreator.h"
#include "TPZAnalyticSolution.h"
#include "TPZHDivApproxCreator.h"

class TPZCompMesh;
class TPZGeoMesh;
class TPZMultiphysicsCompMesh;
class TPZCompEl;
class TPZGeoElSide;

class TPZMHMHDivElasApproxCreator : public TPZHDivApproxCreator, public TPZMHMApproxCreator {

private:
    int fRigidBodyLvl; // Lagrange multiplier level for rigid body motion of the coarse elements
    int fDistForceLevel; // Lagrange multiplier level for distributed forces of the coarse elements
    int fElRigidBodyLvl; // Lagrange multiplier level for rigid body motion of the fine elements
    int fElDistForceLevel; // Lagrange multiplier level for distributed forces of the fine elements
    int fSkeletonOrder;
    int fSkeletonMatId;

public:

    TPZMHMHDivElasApproxCreator(TPZGeoMesh *gmesh_coarse);
    
    TPZMHMHDivElasApproxCreator(TPZGeoMesh* gmesh_fine, TPZVec<int64_t>& elPartition);

    ~TPZMHMHDivElasApproxCreator(){};
    
    virtual TPZMultiphysicsCompMesh* CreateApproximationSpace() override;

    void PutinSubstructures(TPZCompMesh &cmesh);

    void CondenseElements(TPZCompMesh &cmesh);

    /// Set skeleton default polynomial order
    void SetSkeletonOrder(const int ord) {
        fSkeletonOrder = ord;
    }

    /// Get skeleton default polynomial order
    int GetPOrderSkeleton() {return fSkeletonOrder;}

    /// Set skeleton matid
    void SetSkeletonMatId(const int matid) {
        fSkeletonOrder = matid;
    }

    /// Get skeleton matid
    int GetSkeletonMatId() {return fSkeletonMatId;}

    /// Set rigid body level
    void SetRigidBodyLevel(const int level) {
        fRigidBodyLvl = level;
    }
    /// Get rigid body level
    int GetRigidBodyLevel() const {
        return fRigidBodyLvl;
    }
    /// Set distributed force level
    void SetDistForceLevel(const int level) {
        fDistForceLevel = level;
    }
    /// Get distributed force level
    int GetDistForceLevel() const {
        return fDistForceLevel;
    }
    /// Set element rigid body level
    void SetElRigidBodyLevel(const int level) {
        fElRigidBodyLvl = level;
    }
    /// Get element rigid body level
    int GetElRigidBodyLevel() const {
        return fElRigidBodyLvl;
    }
    /// Set element distributed force level
    void SetElDistForceLevel(const int level) {
        fElDistForceLevel = level;
    }
    /// Get element distributed force level
    int GetElDistForceLevel() const {
        return fElDistForceLevel;
    }
};

#endif
