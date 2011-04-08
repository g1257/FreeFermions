

// SAmple of how to use FreeFermions core engine to calculate
// < n_i n_j >

#include "Engine.h"
#include "ObservableLibrary.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
//#include "EtoTheIhTime.h"
//#include "DiagonalOperator.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef std::vector<bool> LevelsType;
typedef PsimagLite::ConcurrencySerial<ComplexType> ConcurrencyType;
typedef PsimagLite::Matrix<ComplexType> MatrixType;
typedef FreeFermions::Engine<RealType,ComplexType,LevelsType,ConcurrencyType> EngineType;
typedef EngineType::HilbertVectorType HilbertVectorType;
typedef typename HilbertVectorType::HilbertTermType HilbertTermType;
typedef typename HilbertVectorType::FlavoredStateType FlavoredStateType;
typedef EngineType::FreeOperatorType FreeOperatorType;
typedef FreeFermions::GeometryLibrary<MatrixType> GeometryLibraryType;
typedef FreeFermions::ObservableLibrary<EngineType> ObservableLibraryType;
//typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
//typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;

int main(int argc,char* argv[])
{
	if (argc!=7) throw std::runtime_error("Needs 6 argument(s)\n");
	size_t n = atoi(argv[1]); // 16 sites
	size_t dof = 2; // spin up and down
	MatrixType t(n,n);
	//std::vector<MatrixType> t;
	GeometryLibraryType geometry(n,GeometryLibraryType::LADDER);
	//GeometryLibraryType geometry(n,GeometryLibraryType::CHAIN);
	//geometry.setGeometry(t,GeometryLibraryType::OPTION_PERIODIC);
	
	//GeometryLibraryType geometry(n,GeometryLibraryType::FEAS);
	//size_t leg = atoi(argv[3]);
	geometry.setGeometry(t,2); //,argv[4],leg); //,GeometryLibraryType::OPTION_PERIODIC);
	
	//std::cerr<<t;
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(t,concurrency,dof,true);
	std::vector<size_t> ne(dof,atoi(argv[2])); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	RealType gsEnergy = 0;
	for (size_t i=0;i<ne[0];i++) gsEnergy += engine.eigenvalue(i);
	gsEnergy *= dof;
	std::cerr<<"Energy="<<gsEnergy<<"\n";
	
	size_t site1 = atoi(argv[3]);
	size_t site2 = atoi(argv[4]);

	enum {SPIN_UP,SPIN_DOWN};

	FreeOperatorType ci = engine.newSimpleOperator("destruction",site1,SPIN_UP);
	HilbertVectorType phiI = engine.newState();
	ci.apply(phiI,gs,FreeOperatorType::SIMPLIFY);

	FreeOperatorType cj = engine.newSimpleOperator("destruction",site2,SPIN_UP);
	HilbertVectorType phiJ = engine.newState();
	cj.apply(phiJ,gs,FreeOperatorType::SIMPLIFY);

	RealType delta = 0.1;
	for (int x=0;x<atoi(argv[5]);x++) {
		RealType omega = x*atof(argv[6]);
		ComplexType z(omega,delta);
		ComplexType sum(0,0);
		for (size_t i=0;i<engine.size();i++) {
			FlavoredStateType nLevel(dof,n);
			nLevel.fill(SPIN_UP,i);
			for (size_t j=0;j<ne[1];j++) nLevel.fill(SPIN_DOWN,j);
			RealType val = 1.0;
			HilbertTermType h(nLevel,val);
			HilbertVectorType nLevelV(n,dof);
			nLevelV.add(h);
			ComplexType tmp1 = scalarProduct(phiI,nLevelV);
			ComplexType tmp2 = scalarProduct(phiJ,nLevelV);
			//std::cout<<"tmp1="<<tmp1<<" tmp2="<<tmp2<<"\n";
			sum += tmp1*tmp2/
				(z + gsEnergy/dof - engine.eigenvalue(i));
		}
		std::cout<<omega<<" "<<std::real(sum)<<" "<<std::imag(sum)<<"\n";
	}
}

