
#ifndef MIXED_ELASTICITY_CMESH_CREATOR
#define MIXED_ELASTICITY_CMESH_CREATOR

#include <iostream>
#include <set>
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZAnalyticSolution.h"

class TPZMixedElasticityCMeshCreator{

public:

    TPZMixedElasticityCMeshCreator() = default;

    //Sets the material ids
    void SetDomainMaterialId(int matid){ fDomainMatId = matid; }
    void SetDirichletMaterialId(std::set<int> &matid){ fBCDirichlet = matid; }
    void SetNeumannMaterialId(std::set<int> &matid){ fBCNeumann = matid; }
    void SetMixedMaterialId(std::set<int> &matid){ fBCMixed = matid; }
    void SetDirichletVarMaterialId(std::set<int> &matid){ fBCDirichletVar = matid; }
    void SetPointMaterialId(std::set<int> &matid){ fBCPoint = matid; }

    /**
     * @brief Funcao para criar a malha computacional da tensao
     * @note Responsavel pela criacao dos espacos de aproximacao do problema
     * @param gmesh malha geometrica
     * @param pOrder ordem polinomial de aproximacao
     */
    TPZCompMesh *CMesh_S(TPZGeoMesh *gmesh, int pOrder, HDivFamily hdivfam, TPZAnalyticSolution * gAnalytic = 0);

    /**
     * @brief Funcao para criar a malha computacional do deslocamento
     * @note Responsavel pela criacao dos espacos de aproximacao do problema
     * @param gmesh malha geometrica
     * @param pOrder ordem polinomial de aproximacao
     */
    TPZCompMesh *CMesh_U(TPZGeoMesh *gmesh, int pOrder);

    /**
     * @brief Funcao para criar a malha computacional da rotacao
     * @note Responsavel pela criacao dos espacos de aproximacao do problema
     * @param gmesh malha geometrica
     * @param pOrder ordem polinomial de aproximacao
     * @param elementdimension tamanho de elementos
     */
    TPZCompMesh *CMesh_P(TPZGeoMesh *gmesh, int pOrder, REAL elementdimension);

    /**
     * @brief Funcao para criar a malha computacional multi-fisica
     * @note Responsavel pela criacao dos espacos de aproximacao do problema
     * @param gmesh malha geometrica
     * @param pOrder ordem polinomial de aproximacao
     */
    TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder, TPZAnalyticSolution * gAnalytic = 0);


    /// change the order of the internal connect to the given order
    void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder);

    void CreateGeometricWrapBCElements(TPZGeoMesh* gmesh);

    TPZCompMesh *CMesh_BCHybrid(TPZGeoMesh *gmesh, int pOrder, TPZAnalyticSolution * gAnalytic = 0);

private:

    int fDomainMatId;
    const int matLagrange = -10; // Material of the Lagrange multipliers
    const int fEWrap = -15;
    const int fEInterface = -16;

    std::set<int> fBCDirichlet = {};
    std::set<int> fBCNeumann = {};
    std::set<int> fBCMixed = {};
    std::set<int> fBCDirichletVar = {};
    std::set<int> fBCPoint = {};

};

#endif
