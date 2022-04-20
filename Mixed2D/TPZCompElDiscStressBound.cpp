#include "TPZCompElDiscStressBound.h"
#include "pzelchdiv.h"
#include "TPZMaterial.h"
#include "pzaxestools.h"

/** @brief Constructor of the discontinuous element associated with geometric element */
template<class StressDef>
TPZCompElDiscStressBound<StressDef>::TPZCompElDiscStressBound(TPZCompMesh &mesh,TPZGeoEl *ref) : TPZCompElDiscStress<StressDef>(mesh,ref),TPZRegisterClassId(&TPZCompElDiscStressBound::ClassId) {
    // this->fConnectIndex = this->CreateMidSideConnect();
    this->fConnectIndex = 0; //Bug here!!! temos que encontrar o numero do conect criado pelo elemento descontinuo.
}

template<class StressDef>
TPZCompElDiscStressBound<StressDef>::TPZCompElDiscStressBound(TPZCompMesh &mesh, const TPZCompElDiscStressBound &copy) :
TPZRegisterClassId(&TPZCompElDiscStressBound::ClassId), TPZCompElDiscStress<StressDef>(mesh,copy){

    this->fConnectIndex = copy.fIndex;
    this->fCenterPoint = copy.fCenterPoint;
    this->fShapefunctionType = copy.fShapefunctionType;
    //criando nova malha computacional
    this->Reference()->SetReference(this);
    //TPZMaterial * mat = copy.Material();
    this->fConstC = copy.fConstC;
    this->fConnectIndex = copy.fConnectIndex;
    this->SetDegree( copy.Degree() );
    //as interfaces foram clonadas
    if (copy.fIntRule){
        this->fIntRule = copy.GetIntegrationRule().Clone();
    }
    else{
        this->fIntRule = NULL;
    }
    this->SetExternalShapeFunction(copy.fExternalShape);
    this->fUseQsiEta = copy.fUseQsiEta;
    this->fStress = copy.fStress;
}


/**
 * @brief Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */
template<class StressDef>
void TPZCompElDiscStressBound<StressDef>::InitMaterialData(TPZMaterialData &data)
{   
    TPZInterpolationSpace::InitMaterialData(data);
    data.XCenter = this->fCenterPoint;

    int nShape = this->fStress->NShapeF();    

    data.fShapeType = TPZMaterialData::EVecShape;

    data.phi.Resize(nShape,2);
    data.divphi.Resize(nShape,1);
    // // não precisa disso. Apenas com base na funcao externa, pega o numero de funções e redimensiona os vetores phi e divphi
    // // Depois no Compute Shape ele vai chamar o método da função externa que calcula as funções e o divergente para o ponto de integração; 
    
}

template<class StressDef>
void TPZCompElDiscStressBound<StressDef>::ComputeShape(TPZVec<REAL> &intpoint,TPZMaterialData &data){

    int nShape = fStress->NShapeF();

    TPZFMatrix<REAL> AiryStress(nShape,3,0.);
    TPZFMatrix<REAL> AiryDivStress(nShape,2,0.);

    TPZFMatrix<REAL> normal(1,3,0.), dnormal(3,3,1.);
    normal(0,0) = data.axes(0,1);
    normal(0,1) = -data.axes(0,0);
    // TPZAxesTools<REAL>::Axes2XYZ(normal,dnormal,data.axes);
    // std::cout << "Mat Id = " << this->Reference()->MaterialId() << std::endl;
    // std::cout << "Axes = " << data.axes << std::endl;
    // std::cout << "Normal = " << normal << std::endl;
    // std::cout << "DNormal = " << dnormal << std::endl;

    fStress->GetStress(data.x,AiryStress,AiryDivStress);
    // std::cout << "X = " << intpoint << std::endl;
    // std::cout << "AiryStress = " << AiryStress << std::endl;
    // std::cout << "AiryDiv = " << AiryDivStress << std::endl;

    

    for (int iShape = 0; iShape < nShape; iShape++)
    {
        data.phi(iShape,0) = AiryStress(iShape,0) * normal[0] + AiryStress(iShape,2) * normal[1];
        data.phi(iShape,1) = AiryStress(iShape,2) * normal[0] + AiryStress(iShape,1) * normal[1];
        // data.phi(iShape,2) = AiryStress(iShape,2);

        data.divphi(iShape,0) = AiryDivStress(iShape,0);
        // data.divphi(iShape,1) = AiryDivStress(iShape,1);
    }
    // std::cout << "phi " << data.phi << std::endl;
}

//External shape retornando a função e seu divergente para este caso.

/**
 * returns the unique identifier for reading/writing objects to streams
 */
template<class StressDef>
int TPZCompElDiscStressBound<StressDef>::ClassId() const{
    return Hash("TPZCompElDiscStressBound") ^ TPZInterpolationSpace::ClassId() << 1;
}


template<class StressDef>
int  TPZCompElDiscStressBound<StressDef>::NShapeF()const{

    return fStress->NShapeF();

}
template<class StressDef>
int TPZCompElDiscStressBound<StressDef>::MaxOrder(){
	int result = TPZInterpolationSpace::MaxOrder();
	result += 3;
	return result;
}

// todos os elementos computacionais da malha Airy terão o mesmo connect.


#include "AiryFunctionLoadedBeam.h"
#include "AiryFunctionHoledPlate.h"

template class TPZCompElDiscStressBound<AiryFunctionLoadedBeam>;
template class TPZCompElDiscStressBound<AiryFunctionHoledPlate>;