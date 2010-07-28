

// SAmple of how to use FreeFermions core engine to calculate
// < n_i n_j >

#include "Includes.h"

int main(int argc,char* argv[])
{
	if (argc!=3) throw std::runtime_error("Needs 2 argument(s)\n");
        size_t n = atoi(argv[1]); 

	size_t dof = 2; // spin up and down
	std::vector<MatrixType> t;
	GeometryLibraryType geometry(n);
	//geometry.setGeometry(t,GeometryLibraryType::CHAIN);
	geometry.setGeometry(t,GeometryLibraryType::FEAS,"feasHoppings.inp",GeometryLibraryType::OPTION_PERIODIC);
	
	EngineType engine(t,dof,true);
	std::vector<size_t> ne(dof,atoi(argv[2])); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	
	ObservableLibraryType library(engine);
	FieldType sum = 0.0, sum2 = 0.0;
	for (size_t site = 0; site<n ; site++) {
		HilbertVectorType phi = engine.newState();
		library.applyNiAllFlavors(phi,gs,site);
		FieldType y = scalarProduct(phi,gs);
		//size_t site2 = site;
		//for (size_t site2=0; site2<n; site2++) {
			//HilbertVectorType phi2 = engine.newState();
			//library.applyNiAllFlavors(phi2,phi,site2);
			FieldType x = scalarProduct(phi,phi);
			std::cout<<x<<" "<<y;
			sum += x;
			sum2 += y*y;
		//}
		std::cout<<"\n";
	}
	std::cout<<"Sum="<<sum<<" "<<sum2<<"\n";
}
