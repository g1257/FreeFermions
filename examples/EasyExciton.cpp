

// <c^dagger_i exp(iHt) n_j exp(-iHt) c_i>

#include "Engine.h"
#include "ObservableLibrary.h"

typedef double RealType;
typedef double FieldType;
typedef size_t UnsignedIntegerType;
typedef FreeFermions::Engine<RealType,FieldType,UnsignedIntegerType> EngineType;
typedef EngineType::HilbertVectorType HilbertVectorType;
typedef EngineType::FreeOperatorType FreeOperatorType;

typedef FreeFermions::ObservableLibrary<EngineType> ObservableLibraryType;



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
	std::cout<<ne;
	size_t flavor =0;
	size_t site = 6;
	FreeOperatorType myOp = engine.newSimpleOperator("destruction",site,flavor);
	HilbertVectorType phi = engine.newState();
	
	myOp.apply(phi,gs);
	FieldType density = scalarProduct(phi,phi);
	std::cerr<<"density="<<density<<"\n";	
	for (size_t site2=0;site2<16;site2++) {
		
		FreeOperatorType myOp2 = engine.newSimpleOperator("destruction",site2,flavor);
		HilbertVectorType phi2 = engine.newState();
		myOp2.apply(phi2,phi);
		//phi2.simplify();
		std::cout<<site2<<" "<<scalarProduct(phi2,phi2)/density<<"\n";
	}
}
