#include "TPZMixedElasticityCMeshCreator.h"

#include "TPZNullMaterial.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "TPZCompelDiscScaled.h"
#include "TPZMixedElasticityNDAiry.h"
#include "pzintel.h"
#include "pzshapequad.h"
using namespace pzshape;
#include "TPZCompElKernelHDiv.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"
#include "pzgeoel.h"
#include "pzgeopoint.h"


TPZCompMesh * TPZMixedElasticityCMeshCreator::CMesh_S(TPZGeoMesh *gmesh, int pOrder, HDivFamily hdivfam, TPZAnalyticSolution * gAnalytic) {
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    int dim = gmesh->Dimension();
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    cmesh->SetDimModel(dim); //Dimesion of the model

    //Definition of the approximation space:
    cmesh->ApproxSpace().SetHDivFamily(hdivfam);
    cmesh->SetAllCreateFunctionsHDiv(); //Creating H(div) functions:

    //Criando material cujo nSTATE = 2:

    TPZNullMaterial<> * material = new TPZNullMaterial<>(fDomainMatId);
    material->SetNStateVariables(dim);
    material->SetDimension(dim);
    cmesh->InsertMaterialObject(material); //Insere material na malha

    TPZNullMaterial<> * matWrap = new TPZNullMaterial<>(fEWrap);
    matWrap->SetNStateVariables(dim);
    matWrap->SetDimension(dim);
    cmesh->InsertMaterialObject(matWrap); //Insere material na malha

    //Boundary conditions:
    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE> val2s(2, 0.), val2(2, 0.);
    // val2s[0] = 10.0; // vx -> 0
    // val2s[1] = 0.0; // vy -> 0


    //Dirichlet Boundary Conditions
    for (auto matId : fBCDirichlet)
    {
        TPZBndCondT<STATE> * BCond = material->CreateBC(material, matId, 0, val1, val2);
        BCond->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
        cmesh->InsertMaterialObject(BCond);
    }

    //Neumann Boundary Conditions
    for (auto matId : fBCNeumann)
    {
        TPZBndCondT<STATE> * BCond = material->CreateBC(material, matId, 1, val1, val2);
        BCond->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
        cmesh->InsertMaterialObject(BCond);
    }

    //Mixed Boundary Conditions
    for (auto matId : fBCMixed)
    {
        DebugStop();
    }
    //Dirichlet Var Boundary Conditions
    for (auto matId : fBCDirichletVar)
    {
        DebugStop();
    }
    //Point Boundary Conditions
    for (auto matId : fBCPoint)
    {
        DebugStop();
    }

    val2s[0] = 0;
    val2s[1] = 0;
    cmesh->InsertMaterialObject(material->CreateBC(material, matLagrange, 0, val1, val2s)); //Insere material na malha

    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    int ncel = cmesh->NElements();
    for (int i = 0; i < ncel; i++) {
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if (!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *> (compEl);
        if (facel)DebugStop();
    }

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;

}




TPZCompMesh * TPZMixedElasticityCMeshCreator::CMesh_U(TPZGeoMesh *gmesh, int pOrder) {
    
    int dim = gmesh->Dimension();

    if (pOrder == 0) {
        // Constant-per-element
        TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
        cmesh->SetDimModel(dim); //Insere dimensão do modelo

        cmesh->SetAllCreateFunctionsDiscontinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);

        TPZNullMaterial<> * material = new TPZNullMaterial<>(fDomainMatId);
        material->SetDimension(dim);
        material->SetNStateVariables(dim);

        cmesh->InsertMaterialObject(material); //Insere material na malha

        std::set<int> materialids;
        materialids.insert(fDomainMatId);
        {
            gmesh->ResetReference();
            int64_t nel = gmesh->NElements();
            for (int64_t el = 0; el < nel; el++) {
                TPZGeoEl *gel = gmesh->Element(el);
                if (!gel)continue;
                int matid = gel->MaterialId();
                if (materialids.find(matid) == materialids.end()) {
                    continue;
                }
                new TPZCompElDiscScaled(*cmesh, gel);
                gel->ResetReference();
            }
        }

        int ncon = cmesh->NConnects();
        for (int i = 0; i < ncon; i++) {
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }

        int64_t nelem = cmesh->NElements();
        for (int64_t el = 0; el < nelem; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZCompElDiscScaled *disc = dynamic_cast<TPZCompElDiscScaled *> (cel);
            if (!disc) {
                continue;
            }
            disc->SetTotalOrderShape();
            disc->SetTrueUseQsiEta();
        }

        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
        return cmesh;
    } else {
        //Criando malha computacional:
        TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
        cmesh->SetDimModel(dim); //Insere dimensão do modelo

        cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
        cmesh->ApproxSpace().CreateDisconnectedElements(true);


        TPZNullMaterial<> * material = new TPZNullMaterial<>(fDomainMatId);
        material->SetDimension(dim);
        material->SetNStateVariables(dim);

        cmesh->InsertMaterialObject(material); //Insere material na malha

        //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha

        int ncel = cmesh->NElements();
        for (int i = 0; i < ncel; i++) {
            TPZCompEl * compEl = cmesh->ElementVec()[i];
            if (!compEl) continue;
            TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *> (compEl);
            if (facel)DebugStop();

        }
        std::set<int> materialids;
        materialids.insert(fDomainMatId);
        cmesh->AutoBuild(materialids);
        cmesh->LoadReferences();
        cmesh->ApproxSpace().CreateDisconnectedElements(false);
        cmesh->AutoBuild();

        int ncon = cmesh->NConnects();
        for (int i = 0; i < ncon; i++) {
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }

        return cmesh;
    }
}



TPZCompMesh * TPZMixedElasticityCMeshCreator::CMesh_P(TPZGeoMesh *gmesh, int pOrder, REAL elementdim) {
    
    int dim = gmesh->Dimension();
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo

    cmesh->SetAllCreateFunctionsDiscontinuous();

    //    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    //Criando material cujo nSTATE = 1:
    TPZNullMaterial<> *material = new TPZNullMaterial<>(fDomainMatId); //criando material que implementa a formulacao fraca do problema modelo
    material->SetDimension(dim);
    if(dim == 3)
    {
        material->SetNStateVariables(3);
    }

    cmesh->InsertMaterialObject(material); //Insere material na malha

    std::set<int> materialids;
    materialids.insert(fDomainMatId);
    //materialids.insert(3);
    {
        gmesh->ResetReference();
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel)continue;
            int matid = gel->MaterialId();
            if (materialids.find(matid) == materialids.end()) {
                continue;
            }
            new TPZCompElDiscScaled(*cmesh, gel);
            gel->ResetReference();
        }
    }

    //cmesh->LoadReferences();
    //    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    //    cmesh->AutoBuild();


    int ncon = cmesh->NConnects();
    for (int i = 0; i < ncon; i++) {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }

    int64_t nelem = cmesh->NElements();
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZCompElDiscScaled *disc = dynamic_cast<TPZCompElDiscScaled *> (cel);
        if (!disc) {
            continue;
        }
        disc->SetTotalOrderShape();
        disc->SetFalseUseQsiEta();
        disc->SetConstC(elementdim);
        disc->SetScale(1./elementdim);
    }

    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    return cmesh;

}



TPZCompMesh * TPZMixedElasticityCMeshCreator::CMesh_m(TPZGeoMesh *gmesh, int pOrder, TPZAnalyticSolution * gAnalytic) {

    int dim = gmesh->Dimension();

    //Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();

    // Criando material:

    // example is initialized in the calling method
    //    example.fProblemType = TElasticityExample1::EThiago;
    //    example.fStressState   = TElasticityExample1::EPlaneStrain;

    REAL E = 0.; //* @param E elasticity modulus
    REAL nu = 0.; //* @param nu poisson coefficient

    TElasticity2DAnalytic * analytic2D = 0;
    TElasticity3DAnalytic * analytic3D = 0;
    if(dim == 2)
    {
        analytic2D = dynamic_cast<TElasticity2DAnalytic *>(gAnalytic);
        if(!analytic2D) DebugStop();
        TPZVec<REAL> x;
        x.Resize(3);
        x.Fill(0.);
        analytic2D->Elastic<REAL>(x, E, nu);
    }
    else if(dim == 3)
    {
        analytic3D = dynamic_cast<TElasticity3DAnalytic *>(gAnalytic);
        if(!analytic3D) DebugStop();
        TPZManVector<REAL,3> x(3,0.);
        analytic3D->Elastic(x, E, nu);
    }

    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = 0.; //* @param fx forcing function \f$ -x = fx \f$
    int plainStress = 0; //* @param plainstress = 1 \f$ indicates use of plainstress
    if(dim == 2)
    {
        if (analytic2D->fPlaneStress == 0) {
            plainStress = 0;
        } else {
            plainStress = 1;
        }
    }
    // TPZMixedElasticityND * material = new TPZMixedElasticityND(matID, E, nu, fx, fy, plainStress, dim);
    TPZMixedElasticityNDAiry * material = new TPZMixedElasticityNDAiry(fDomainMatId, E, nu, fx, fy, plainStress, dim);

    material->SetForcingFunction(gAnalytic->ForceFunc(),4);
    material->SetExactSol(gAnalytic->ExactSolution(),4);     

    cmesh->InsertMaterialObject(material);


    //Condições de contorno:

    TPZFMatrix<REAL> val1(dim, dim, 0.);
    TPZManVector<REAL> val2(dim, 0.);
    

    //Dirichlet Boundary Conditions
    for (auto matId : fBCDirichlet)
    {
        TPZBndCondT<STATE> * BCond = material->CreateBC(material, matId, 0, val1, val2);
        BCond->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
        cmesh->InsertMaterialObject(BCond);
    }

    //Neumann Boundary Conditions
    for (auto matId : fBCNeumann)
    {
        TPZBndCondT<STATE> * BCond = material->CreateBC(material, matId, 1, val1, val2);
        BCond->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
        cmesh->InsertMaterialObject(BCond);
    }

    //Mixed Boundary Conditions
    for (auto matId : fBCMixed)
    {
        DebugStop();
    }
    //Dirichlet Var Boundary Conditions
    for (auto matId : fBCDirichletVar)
    {
        DebugStop();
    }
    //Point Boundary Conditions
    for (auto matId : fBCPoint)
    {
        DebugStop();
    }


    auto * BCond7 = material->CreateBC(material, matLagrange, 0, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond7->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond7); //Insere material na malha


    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;

}




/// change the order of the internal connect to the given order

void TPZMixedElasticityCMeshCreator::ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder) {
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) continue;

        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != cmesh->Dimension()) {
            continue;
        }
        int nc = cel->NConnects();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel) DebugStop();

        intel->ForceSideOrder(gel->NSides() - 1, pOrder);

        if (cmesh->ApproxSpace().HDivFam() == HDivFamily::EHDivKernel){
            TPZCompElKernelHDiv<TPZShapeQuad> *celKernel = dynamic_cast<TPZCompElKernelHDiv<TPZShapeQuad> *> (cel);
            for (int icon = 0; icon < celKernel->NConnects(); icon++)
            {
                TPZConnect &c = celKernel->Connect(icon);
                int nShapeF = celKernel->NConnectShapeF(icon,c.Order());
                c.SetNShape(nShapeF);
                int64_t seqnum = c.SequenceNumber();
                int nvar = 1;
                TPZMaterial * mat = celKernel->Material();
                if (mat) nvar = mat->NStateVariables();
                celKernel->Mesh()->Block().Set(seqnum, nvar * nShapeF);
                celKernel->AdjustIntegrationRule();
            }    
        }
        
    }
    cmesh->ExpandSolution();
}




void TPZMixedElasticityCMeshCreator::CreateGeometricWrapBCElements(TPZGeoMesh* gmesh){

    std::set<int> matIdBC = fBCNeumann;

    int dim = gmesh->Dimension();

    //Disconnect elements
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) DebugStop();
        auto type = gel -> Type();
        auto matid = gel->MaterialId();

        using namespace pzgeom;
        using namespace pzshape;

        if (gel->Dimension() != dim) continue;
        
        int nsides = gel->NSides();

        for (int side = 0; side < nsides; side++) {

            if(gel->SideDimension(side) != dim-1) continue;// onlu edges for 2D and faces for 3D
            TPZGeoElSide geoside(gel,side);
            TPZGeoElSide neighbour = geoside.Neighbour();

        
            //Creates interface and wrap geometric elements for hybridized BC
            if (matIdBC.find(neighbour.Element()->MaterialId()) != matIdBC.end()){
                TPZGeoElBC gelbcWrap(geoside, fEWrap);
                TPZGeoElSide gelWrapSide(gelbcWrap.CreatedElement(),gelbcWrap.CreatedElement()->NSides()-1);
                TPZGeoElBC gelbc(gelWrapSide, fEInterface);
                gelbcWrap.CreatedElement()->ResetReference();
                gelbc.CreatedElement()->ResetReference();
                gel->ResetReference();
            }
            
        }

    }

}


TPZCompMesh *CMesh_BCHybrid(TPZGeoMesh *gmesh, int pOrder, TPZAnalyticSolution * gAnalytic){

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    DebugStop();


    return cmesh;
}