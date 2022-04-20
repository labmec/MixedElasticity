#include "TPZCompElDiscStress.h"
#include "pzelchdiv.h"
#include "TPZMaterial.h"

/** @brief Constructor of the discontinuous element associated with geometric element */
template<class StressDef>
TPZCompElDiscStress<StressDef>::TPZCompElDiscStress(TPZCompMesh &mesh,TPZGeoEl *ref) : TPZCompElDisc(mesh,ref),TPZRegisterClassId(&TPZCompElDiscStress::ClassId) {
    
    this->fConnectIndex = 0;
    this->SetConnectIndex(0,0);

}

template<class StressDef>
TPZCompElDiscStress<StressDef>::TPZCompElDiscStress(TPZCompMesh &mesh, const TPZCompElDiscStress &copy) :
TPZRegisterClassId(&TPZCompElDiscStress::ClassId), TPZCompElDisc(mesh,copy){

    this->fConnectIndex = copy.fIndex;
    this->fCenterPoint = copy.fCenterPoint;
    fShapefunctionType = copy.fShapefunctionType;
    //criando nova malha computacional
    Reference()->SetReference(this);
    //TPZMaterial * mat = copy.Material();
    fConstC = copy.fConstC;
    fConnectIndex = copy.fConnectIndex;
    this->fStress = copy.fStress;
    this->SetDegree( copy.Degree() );
    //as interfaces foram clonadas
    if (copy.fIntRule){
        this->fIntRule = copy.GetIntegrationRule().Clone();
    }
    else{
        this->fIntRule = NULL;
    }
    this->SetExternalShapeFunction(copy.fExternalShape);
    fUseQsiEta = copy.fUseQsiEta;
    
}


/**
 * @brief Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */
template<class StressDef>
void TPZCompElDiscStress<StressDef>::InitMaterialData(TPZMaterialData &data)
{   
    TPZInterpolationSpace::InitMaterialData(data);
    data.XCenter = fCenterPoint;

    int nShape = this->fStress->NShapeF();    

    data.fShapeType = TPZMaterialData::ETensorShape;

    data.phi.Resize(nShape,3);
    data.divphi.Resize(nShape,2);
    // // não precisa disso. Apenas com base na funcao externa, pega o numero de funções e redimensiona os vetores phi e divphi
    // // Depois no Compute Shape ele vai chamar o método da função externa que calcula as funções e o divergente para o ponto de integração; 
    
}

template<class StressDef>
void TPZCompElDiscStress<StressDef>::ComputeShape(TPZVec<REAL> &intpoint,TPZMaterialData &data){

    int nShape = fStress->NShapeF();

    TPZFMatrix<REAL> AiryStress(nShape,3,0.);
    TPZFMatrix<REAL> AiryDivStress(nShape,2,0.);
    this->Reference()->X(intpoint,data.x);

    fStress->GetStress(data.x,AiryStress,AiryDivStress);
    // std::cout << "X = " << intpoint << std::endl;
    // std::cout << "AiryStress = " << AiryStress << std::endl;
    // std::cout << "AiryDiv = " << AiryDivStress << std::endl;
    for (int iShape = 0; iShape < nShape; iShape++)
    {
        data.phi(iShape,0) = AiryStress(iShape,0);
        data.phi(iShape,1) = AiryStress(iShape,1);
        data.phi(iShape,2) = AiryStress(iShape,2);

        data.divphi(iShape,0) = AiryDivStress(iShape,0);
        data.divphi(iShape,1) = AiryDivStress(iShape,1);
    }
    // std::cout << "COmpute shape x = " << data.x << std::endl;
    // std::cout << "COmpute shape Phi = " << data.phi << std::endl;
}

//External shape retornando a função e seu divergente para este caso.

/**
 * returns the unique identifier for reading/writing objects to streams
 */
template<class StressDef>
int TPZCompElDiscStress<StressDef>::ClassId() const{
    return Hash("TPZCompElDiscStress") ^ TPZInterpolationSpace::ClassId() << 1;
}


template<class StressDef>
int  TPZCompElDiscStress<StressDef>::NShapeF()const{

    return fStress->NShapeF();

}

template<class StressDef>
void TPZCompElDiscStress<StressDef>::ReallyComputeSolution(TPZMaterialDataT<STATE>& data){
    const TPZFMatrix<REAL> &phi = data.phi;
    const TPZFMatrix<REAL> &dphix = data.dphix;
    const TPZFMatrix<REAL> &axes = data.axes;
    TPZSolVec<STATE> &sol = data.sol;
    TPZGradSolVec<STATE> &dsol = data.dsol;
	const int nstate = this->Material()->NStateVariables();
	const int ncon = this->NConnects();
	TPZBlock &block = Mesh()->Block();
	TPZFMatrix<STATE> &MeshSol = Mesh()->Solution();
    const int64_t numbersol = MeshSol.Cols();
	// std::cout << "Meshsol" << MeshSol << std::endl;
    const int nTens = data.phi.Cols();
	const int solVecSize = ncon? nTens : 0;
    
    if (nstate != 1){
        DebugStop();//Stress variables and stress solution, 
    }

    sol.resize(numbersol);
    dsol.resize(numbersol);
    for (int is = 0; is<numbersol; is++) {
        sol[is].Resize(solVecSize);
        sol[is].Fill(0.);
        dsol[is].Redim(dphix.Rows(), solVecSize);
        dsol[is].Zero();
    }	
	int64_t iv = 0;


	for(int in=0; in<ncon; in++) {
		TPZConnect *df = &Connect(in);
		const int64_t dfseq = df->SequenceNumber();
		const int dfvar = block.Size(dfseq);
		const int64_t pos = block.Position(dfseq);
		for(int jn=0; jn<nTens; jn++) {
            for (int64_t is=0; is<numbersol; is++) {
                sol[is][iv/nstate] +=
                    (STATE)phi.Get(0,iv/nstate)*MeshSol(pos,is);
                // for(auto d=0; d<dphix.Rows(); d++){
                //     dsol[is](d,iv%nstate) +=
                //         (STATE)dphix.Get(d,iv/nstate)*MeshSol(pos,is);
                // }
            }
			iv++;
		}
	}
	
}//method


template<class StressDef>
int TPZCompElDiscStress<StressDef>::NConnectShapeF(int inod, int order) const {

	return this->NShapeF();
}

template<class StressDef>
int TPZCompElDiscStress<StressDef>::MaxOrder(){
	int result = TPZInterpolationSpace::MaxOrder();
	result += 3;
	return result;
}

#include "AiryFunctionLoadedBeam.h"
#include "AiryFunctionHoledPlate.h"

template class TPZCompElDiscStress<AiryFunctionLoadedBeam>;
template class TPZCompElDiscStress<AiryFunctionHoledPlate>;