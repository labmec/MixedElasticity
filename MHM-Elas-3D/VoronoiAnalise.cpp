//
//  VoronoiAnalise.cpp
//  mhm-elas-3d
//
//  Created by Philippe Devloo on 21/08/23.
//

#include "VoronoiAnalise.h"
#include "pzcmesh.h"
#include "TPZShapeDisc.h"
#include "pzvec_extras.h"

#include "pzcompel.h"
#include "pzmultiphysicselement.h"
#include "pzconnect.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "tpztriangle.h"
#include "pzquad.h"

#include "TPZMaterialDataT.h"


void TPZFaceDefinition::Shape(const TPZVec<REAL> &x, TPZFMatrix<REAL> &phi)
{

    
//    static void Shape(int dimension, REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, MShapeType type);
    int dim = 2;
    TPZManVector<REAL,3> delx = x-fCenter;
    TPZFMatrix<REAL> delxM(3,1,&delx[0],3);
    TPZFNMatrix<6,REAL> coord(2,1);
    fAxes.Multiply(delxM, coord);
    TPZManVector<REAL,2> xplane(2);
    for (int i=0; i<2; i++) {
        xplane[i] = coord(i);
    }
    int nshape =  pzshape::TPZShapeDisc::NShapeF(fPOrder, 2, pzshape::TPZShapeDisc::EOrdemTotal);
    TPZFNMatrix<24,REAL> dphi(2,nshape);
    phi.Resize(nshape, 1);
    pzshape::TPZShapeDisc::Shape(2, fScale, fCenter, xplane, fPOrder, phi, dphi, pzshape::TPZShapeDisc::EOrdemTotal);
}

int TPZFaceDefinition::NShape() {
    int nshape =  pzshape::TPZShapeDisc::NShapeF(fPOrder, 2, pzshape::TPZShapeDisc::EOrdemTotal);
    return nshape;
}
               // project and restrain the shape functions
void TPZFaceDefinition::Project(TPZCompEl *cel) {
    TPZGeoEl *gel = cel->Reference();
    if(gel->Type() != ETriangle) DebugStop();
    if(cel->NConnects() != 1) DebugStop();
    TPZManVector<REAL,3> X(3);
    TPZFNMatrix<9> gradx(3,2),jac(2,2),jacinv(2,2),axes(2,3);
    REAL detjac;
    int nshape = NShape();
    TPZFNMatrix<10> phi(nshape,1);
    TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
    // get the hdiv element
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(mphys->Element(0));
    if(!intel) DebugStop();
    TPZMaterialDataT<REAL> data;
    intel->InitMaterialData(data);
    int64_t nshapeL2 = data.phi.Rows();
    TPZFNMatrix<9,REAL> matL2(nshapeL2,nshapeL2,0.);
    TPZFNMatrix<9,REAL> L2proj(nshapeL2,nshape,0.);
    
    pztopology::TPZTriangle::IntruleType therule(2*2);
    int np = therule.NPoints();
    for (int ip = 0; ip<np; ip++) {
        TPZManVector<REAL> qsi(2);
        REAL weight;
        therule.Point(ip, qsi, weight);
        intel->ComputeRequiredData(data,qsi);
        gel->X(qsi,X);
        gel->GradX(qsi, gradx);
        gel->Jacobian(gradx, jac, axes, detjac, jacinv);
        Shape(X, phi);
        for(int i=0; i<nshapeL2; i++) {
            for (int j=0; j<nshapeL2; j++) {
                matL2(i,j) += data.phi(i,0)*data.phi(j,0)*detjac*weight;
            }
            for (int j=0; j<nshape; j++) {
                L2proj(i,j) += data.phi(i,0)*phi(j,0)*detjac*weight;
            }
        }
    }
//    matL2.Print("MatL2");
    matL2.SolveDirect(L2proj, ELU);
    
    TPZConnect &c = cel->Connect(0);
    
    c.AddDependency(cel->ConnectIndex(0), fConnectIndex, L2proj, 0, 0, (int)nshapeL2, nshape);
}

static void ComputeNormal(TPZGeoEl *gel, TPZVec<REAL> &normal) {
    TPZManVector<REAL,2> center(2);
    TPZFNMatrix<9> gradx(3,2),jac(2,2),axes(2,3),jacinv(2,2);
    REAL detjac;
    gel->CenterPoint(gel->NSides()-1, center);
    gel->GradX(center, gradx);
    gel->Jacobian(gradx, jac, axes, detjac, jacinv);
    for(int i=0; i<3; i++) {
        normal[i] = axes(0,(i+1)%3)*axes(1,(i+2)%3)-axes(1,(i+1)%3)*axes(0,(i+2)%3);
    }
}

static void ComputeAxes(TPZGeoEl *gel, TPZFMatrix<REAL> &axes) {
    TPZManVector<REAL,2> center(2);
    TPZFNMatrix<9> gradx(3,2),jac(2,2),jacinv(2,2);
    REAL detjac;
    gel->CenterPoint(gel->NSides()-1, center);
    gel->GradX(center, gradx);
    gel->Jacobian(gradx, jac, axes, detjac, jacinv);
}
// Initialize geometric datastructure
void TPZFaceDefinition::InitializeDataStructure(const int pord)
{
    if(fElementList.size() == 0) DebugStop();
    fPOrder = pord;
    TPZCompEl *cel = *fElementList.begin();
    TPZGeoEl *gel = cel->Reference();
    TPZGeoElSide gelside(gel);
    TPZManVector<REAL,3> center(2),normal(3);
    gelside.CenterPoint(center);
    fNormal.resize(3);
    ComputeNormal(gel, fNormal);
    fAxes.Redim(2, 3);
    ComputeAxes(gel, fAxes);
    ComputeCenter();
    TPZCompMesh *cmesh = cel->Mesh();
    int64_t newcon = cmesh->AllocateNewConnect(NShape(), 3, fPOrder);
    fConnectIndex = newcon;
}

/// Compute the center of of the polygon and its area
void TPZFaceDefinition::ComputeCenter() {
    TPZManVector<REAL,3> Xacc(3,0.);
    REAL totalarea = 0.;
    for (auto it : fElementList) {
        TPZCompEl *cel = it;
        TPZGeoEl *gel = cel->Reference();
        REAL area = gel->Volume();
        TPZGeoElSide gelside(gel);
        TPZManVector<REAL,3> centerx(3,0.);
        gelside.CenterX(centerx);
        for(int i = 0; i<3; i++) {
            Xacc[i] += centerx[i]*area;
        }
        totalarea += area;
    }
    fCenter.Resize(3, 0.);
    for (int i=0; i<3; i++) {
        fCenter[i] = Xacc[i]/totalarea;
    }
    fScale = sqrt(totalarea);
    fTotalArea = totalarea;
}
