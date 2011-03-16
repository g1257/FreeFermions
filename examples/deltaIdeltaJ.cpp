

// SAmple of how to use FreeFermions core engine to calculate
// < \delta_{i,\gamma} \delta^\dagger_{j,\gamma'} >

#include "Engine.h"
#include "ObservableLibrary.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
//#include "EtoTheIhTime.h"
//#include "DiagonalOperator.h"

typedef double RealType;
typedef std::complex<double> FieldType;
typedef std::vector<bool> LevelsType;
typedef Dmrg::ConcurrencySerial<FieldType> ConcurrencyType;
typedef PsimagLite::Matrix<FieldType> MatrixType;
typedef FreeFermions::Engine<RealType,FieldType,LevelsType,ConcurrencyType> EngineType;
typedef EngineType::HilbertVectorType HilbertVectorType;
typedef EngineType::FreeOperatorType FreeOperatorType;
typedef FreeFermions::GeometryLibrary<MatrixType> GeometryLibraryType;
typedef FreeFermions::ObservableLibrary<EngineType> ObservableLibraryType;
//typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
//typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;

int main(int argc,char* argv[])
{
	if (argc!=4) throw std::runtime_error("Needs 3 argument(s)\n");
        size_t n = atoi(argv[1]); 

	size_t dof = 2; // spin up and down
	std::vector<MatrixType> t;
	GeometryLibraryType geometry(n,GeometryLibraryType::FEAS);
	//geometry.setGeometry(t,GeometryLibraryType::CHAIN);
	size_t leg = atoi(argv[3]);
	geometry.setGeometry(t,"feasHoppings.inp",leg); //,GeometryLibraryType::OPTION_PERIODIC);
	
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(t,concurrency,dof,false);
	std::vector<size_t> ne(dof,atoi(argv[2])); // argv[2] up and argv[2] down
	HilbertVectorType gs = engine.newGroundState(ne);
	
	ObservableLibraryType library(engine);
	MatrixType m(n,n);
	FieldType sum2 = 0;
	size_t effectiveN  = n;
	for (size_t site = 0; site<effectiveN ; site++) {
		FieldType y = 0;
		size_t orb1 = 0;
		//for (size_t orb1 = 0;orb1<2; orb1++) {
			HilbertVectorType phi = engine.newState();
			library.applyDelta(phi,gs,site+orb1*n);  // phi = operator |gs>
			y += scalarProduct(phi,gs);
			size_t orb2 = orb1;
			for (size_t site2=0; site2<effectiveN; site2++) {
				//for (size_t orb2=0;orb2<2;orb2++) {
					HilbertVectorType phi2 = engine.newState();
					library.applyDelta(phi2,gs,site2+orb2*n); // phi2 = operator2 |phi>
					FieldType x = scalarProduct(phi2,phi);
					std::cout<<x<<" ";
					std::cout.flush();
					m(site,site2) += x;
					sum2 += y*y;
				//}
			}
			std::cout<<"\n";
		//}
	}
}
