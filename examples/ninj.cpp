

// SAmple of how to use FreeFermions core engine to calculate
// < n_i n_j >

#include "Engine.h"
#include "ObservableLibrary.h"

typedef double FieldType;
typedef FreeFermions::Engine<FieldType> EngineType;
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
	
	ObservableLibraryType library(engine);
	
	for (size_t site = 0; site<n ; site++) {
		for (size_t site2=0; site2<n; site2++) {
			HilbertVectorType phi = engine.newState();
			library.applyNiAllFlavors(phi,gs,site);
			HilbertVectorType phi2 = engine.newState();
			library.applyNiAllFlavors(phi2,phi,site2);
			std::cout<<scalarProduct(phi2,gs)<<" ";
		}
		std::cout<<"\n";
	}
}
