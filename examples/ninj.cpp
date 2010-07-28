

// SAmple of how to use FreeFermions core engine to calculate
// < n_i n_j >

#include "Includes.h"

int main(int argc,char* argv[])
{
	if (argc!=3) throw std::runtime_error("Needs 2 argument(s)\n");
        size_t n = atoi(argv[1]); 

	size_t dof = 2; // spin up and down
	MatrixType t(n,n);
	GeometryLibraryType geometry(n,GeometryLibraryType::LADDER);
	//geometry.setGeometry(t,GeometryLibraryType::CHAIN);
	geometry.setGeometry(t,2);
	
	EngineType engine(t,dof,false);
	std::vector<size_t> ne(dof,atoi(argv[2])); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	ObservableLibraryType library(engine);
	
	for (size_t site = 0; site<n ; site++) {
		HilbertVectorType phi = engine.newState();
		library.applyNiAllFlavors(phi,gs,site);
		for (size_t site2=0; site2<n; site2++) {
			HilbertVectorType phi2 = engine.newState();
			library.applyNiAllFlavors(phi2,phi,site2);
			std::cout<<scalarProduct(phi2,gs)<<" ";
		}
		std::cout<<"\n";
	}
}
