

// SAmple of how to use FreeFermions core engine to calculate
// < n_i n_j >

#include "Engine.h"
#include "ObservableLibrary.h"

typedef double RealType;
typedef double FieldType;
typedef size_t UnsignedIntegerType;
typedef FreeFermions::Engine<RealType,FieldType,UnsignedIntegerType> EngineType;
typedef EngineType::HilbertVectorType HilbertVectorType;
typedef EngineType::FreeOperatorType FreeOperatorType;

//typedef FreeFermions::ObservableLibrary<EngineType> ObservableLibraryType;

int main()
{
	size_t n = 16; // 16 sites
	size_t dof = 2; // spin up and down
	bool isPeriodic = false;
	psimag::Matrix<FieldType> t(n,n);
	for (size_t i=0;i<n;i++) {
		for (size_t j=0;j<n;j++) {
			if (i-j==1 || j-i==1) t(i,j) = 1.0;
		}
	}
	if (isPeriodic) t(0,n-1) = t(n-1,0) = 1.0;
	
	EngineType engine(t,dof,false);
	std::vector<size_t> ne(dof,8); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	
	//ObservableLibraryType library(engine);
	size_t sigma = 0;
	for (size_t site = 0; site<n ; site++) {
		HilbertVectorType phi = engine.newState();
		FreeOperatorType myOp = engine.newSimpleOperator("destruction",site,sigma);
		myOp.apply(phi,gs);
		for (size_t site2=0; site2<n; site2++) {
			HilbertVectorType phi2 = engine.newState();
			FreeOperatorType myOp2 = engine.newSimpleOperator("creation",site2,sigma);
			myOp2.apply(phi2,phi);
			std::cout<<scalarProduct(phi2,gs)<<" ";
		}
		std::cout<<"\n";
	}
}
