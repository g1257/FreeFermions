

// SAmple of how to use FreeFermions core engine to calculate
// < n_i n_j >

#include "Includes.h"

int main(int argc,char* argv[])
{
	if (argc!=3) throw std::runtime_error("Needs 2 argument(s)\n");
        size_t n = atoi(argv[1]); 

	size_t dof = 2; // spin up and down
	std::vector<MatrixType> t;
	GeometryLibraryType geometry(n,GeometryLibraryType::FEAS);
	//geometry.setGeometry(t,GeometryLibraryType::CHAIN);
	geometry.setGeometry(t,"feasHoppings.inp",GeometryLibraryType::OPTION_PERIODIC);
	
	EngineType engine(t,dof,true);
	std::vector<size_t> ne(dof,atoi(argv[2])); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	
	ObservableLibraryType library(engine);
	MatrixType m(n,n);
	FieldType sum2 = 0;
	for (size_t site = 0; site<n ; site++) {
		FieldType y = 0;
		for (size_t orb1 = 0;orb1<2; orb1++) {
			HilbertVectorType phi = engine.newState();
			library.applyNiAllFlavors(phi,gs,site+orb1*n);
			y += scalarProduct(phi,gs);
			for (size_t site2=0; site2<n; site2++) {
				for (size_t orb2=0;orb2<2;orb2++) {
					HilbertVectorType phi2 = engine.newState();
					library.applyNiAllFlavors(phi2,phi,site2+orb2*n);
					FieldType x = scalarProduct(phi2,gs);
					std::cout<<x<<" ";
					m(site,site2) += x;
					sum2 += y*y;
				}
			}
			std::cout<<"\n";
		}
	}
	std::vector<std::complex<FieldType> > fm(n);
	geometry.fourierTransform(fm,m);
	std::cout<<"Fourier transform:\n";
	std::cout<<fm;
}
