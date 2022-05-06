#include <iostream>

#include "pzlog.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZGenGrid3D.h"
#include "Elasticity/TPZElasticity3D.h"
#include "TPZAnalyticSolution.h"
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>
#include <pzstepsolver.h>


using namespace std;

enum EMatId {ENone, EVolume, EFaceBC};

int main() {
	
#ifdef PZ_LOG
	TPZLogger::InitializePZLOG();
#endif
	
	std::vector<int> refvec = {2,4,8,16};
	const int nref = refvec.size();
	
	std::vector<REAL> L2err(nref), H1err(nref);
	
	for (int iref = 0; iref < nref; iref++) {
		const int ref = refvec[iref];
		cout << "\n\t----------- Starting Convergence with mesh size " << ref << " -----------\n" << endl;
		
		// -------------------- Create GMesh --------------------
		TPZGeoMesh* gmesh = nullptr;
		const TPZVec<REAL> minX = {-1.,-1.,-1.}, maxX = {1.,1.,1.};
		const TPZVec<int> nelDiv = {ref,ref,ref};
		const MMeshType elType = MMeshType::ETetrahedral;
		TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
		gmesh = gen3d.BuildVolumetricElements(EVolume);
		gmesh = gen3d.BuildBoundaryElements(EFaceBC,EFaceBC,EFaceBC,EFaceBC,EFaceBC,EFaceBC);
		gmesh->BuildConnectivity();
		const int dim = gmesh->Dimension();
		
		// -------------------- Create Analytical Sol --------------------
		TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
		elas->fE = 1.;
		elas->fPoisson = 0.3;
		elas->fProblemType = TElasticity3DAnalytic::EYotov;
		
		// -------------------- Create CompMesh --------------------
		TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
		const REAL Ela = 1., poisson = 0.3;
		TPZVec<STATE> force(3,0.);
		TPZElasticity3D* material = new TPZElasticity3D(EVolume, Ela, poisson, force);
		material->SetForcingFunction(elas->ForceFunc(),3);
		cmesh->InsertMaterialObject(material);
		const int dirichlet = 0;
		TPZFMatrix<REAL> val1(dim, dim, 0.);
		TPZManVector<REAL> val2(dim, 0.);
		auto bnd = material->CreateBC(material, EFaceBC, dirichlet, val1, val2);
		bnd->SetForcingFunctionBC(elas->ExactSolution(),3);
		cmesh->InsertMaterialObject(bnd);
		cmesh->SetDefaultOrder(1);
		cmesh->AutoBuild();
		
		
		// -------------------- Create Analysis --------------------
		cout << "----------- Opt Band -----------" << endl;
		cout << "\nNumber of elements: " << cmesh->NElements() << endl;
		cout << "Number of equations: " << cmesh->NEquations() << endl;
		auto start_time_anal = std::chrono::steady_clock::now();
		TPZLinearAnalysis an(cmesh,true);
		constexpr int nThreads{8};
		TPZSSpStructMatrix<STATE> matskl(cmesh); // fast - works great with mkl
		matskl.SetNumThreads(nThreads);
		an.SetStructuralMatrix(matskl);
		TPZStepSolver<STATE> step;
		step.SetDirect(ELDLt);
		an.SetSolver(step);
		auto total_time_anal = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_anal).count()/1000.;
		cout << "Total time Opt Band = " << total_time_anal << " seconds" << endl;
		
		
		// -------------------- Solve Problem --------------------
		cout << "\n----------- Solving -----------" << endl;
		auto start_time_ass = std::chrono::steady_clock::now();
		cout << "\nDoing assemble..." << endl;
		an.Assemble();
		auto total_time_ass = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_ass).count()/1000.;
		cout << "Total time assemble = " << total_time_ass << " seconds" << endl;
		
		auto start_time_solve = std::chrono::steady_clock::now();
		cout << "\nDoing solve..." << endl;
		an.Solve();
		auto total_time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_solve).count()/1000.;
		cout << "Total time solve = " << total_time_solve << " seconds" << endl;
		
		// -------------------- Post Process --------------------
		TPZManVector<std::string,1> scalnames(0), vecnames(1);
		vecnames[0] = "Displacement";
		int div = 0;
		string plotfile = "PostProcess.vtk";
		an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
		an.PostProcess(div,dim);
		
		// -------------------- Compute Error --------------------
		cout << "\n----------- Computing Error -----------" << endl;
		TPZMatErrorSingleSpace<STATE> *materr = dynamic_cast<TPZMatErrorSingleSpace<STATE> *>(cmesh->FindMaterial(EVolume));
		an.SetExact(elas->ExactSolution());
		an.SetThreadsForError(nThreads);
		TPZManVector<REAL, 10> Errors(materr->NEvalErrors());
		const bool store_errors = false;
		an.PostProcessError(Errors, store_errors);
		std::cout << "Errors = " << Errors << std::endl;
		
		L2err[iref] = Errors[1];
		H1err[iref] = Errors[2];
		
		delete elas;
		delete cmesh;
		delete gmesh;
	}
	
	cout << "\n----------- Finished! -----------" << endl;
	cout << "Errors:\n" << endl;
	cout << "ErrorsL2 = {" << L2err[0];
	
	for (int i = 1; i < nref; i++) {
		cout << "," << L2err[i];
	}
	cout << "};" << endl;
	
	cout << "ErrorsH1 = {" << H1err[0];
	for (int i = 1; i < nref; i++) {
		cout << "," << H1err[i];
	}
	cout << "};\n" << endl;

	return 0;
}
