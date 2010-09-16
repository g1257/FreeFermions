

// SAmple of how to use FreeFermions core engine to calculate
// < n_i n_j >

#include "Includes.h"

int main(int argc,char* argv[])
{
	if (argc!=4) throw std::runtime_error("Needs 3 argument(s)\n");
        size_t n = atoi(argv[1]); 

	size_t dof = 2; // spin up and down
	std::vector<MatrixType> t;
	GeometryLibraryType geometry(n,GeometryLibraryType::FEAS);
	//geometry.setGeometry(t,GeometryLibraryType::CHAIN);
	size_t leg = atoi(argv[3]);
	geometry.setGeometry(t,"feasHoppings.inp",leg,GeometryLibraryType::OPTION_PERIODIC);
	
	EngineType engine(t,dof,true);
	std::vector<size_t> ne(dof,atoi(argv[2])); // argv[2] up and argv[2] down
	HilbertVectorType gs = engine.newGroundState(ne);
	
	ObservableLibraryType library(engine);
	MatrixType m(n,n);
	FieldType sum2 = 0;
	size_t effectiveN  = n;
	for (size_t site = 0; site<effectiveN ; site++) {
		FieldType y = 0;
		for (size_t orb1 = 0;orb1<2; orb1++) {
			HilbertVectorType phi = engine.newState();
			library.applyNiAllFlavors(phi,gs,site+orb1*n); // phi = operator |gs>
			y += scalarProduct(phi,gs);
			for (size_t site2=0; site2<effectiveN; site2++) {
				for (size_t orb2=0;orb2<2;orb2++) {
					HilbertVectorType phi2 = engine.newState();
					library.applyNiAllFlavors(phi2,gs,site2+orb2*n); // phi2 = operator2 |phi>
					FieldType x = scalarProduct(phi2,phi);
					std::cout<<x<<" ";
					std::cout.flush();
					m(site,site2) += x;
					sum2 += y*y;
				}
			}
			std::cout<<"\n";
		}
	}
	std::vector<std::complex<FieldType> > fm(n);
	geometry.fourierTransform(fm,m,leg);
	std::cout<<"Fourier transform:\n";
	std::cout<<fm;
}
