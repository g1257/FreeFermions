// SAmple of how to use FreeFermions core engine to calculate
// < n_i n_j >

#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "LibraryOperator.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef RealType FieldType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::Engine<RealType,FieldType,ConcurrencyType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::HilbertState<OperatorType> HilbertStateType;
typedef FreeFermions::GeometryLibrary<MatrixType> GeometryLibraryType;
typedef FreeFermions::LibraryOperator<OperatorType> LibraryOperatorType;

int main(int argc,char* argv[])
{
	if (argc!=3) throw std::runtime_error("Needs 2 argument(s)\n");
        size_t n = atoi(argv[1]); 

	size_t dof = 2; // spin up and down
	MatrixType t(n,n);
	//GeometryLibraryType geometry(n,GeometryLibraryType::LADDER);
	GeometryLibraryType geometry(n,GeometryLibraryType::CHAIN);
	geometry.setGeometry(t); //,2);
	//t(3,5) = t(5,3) = 0;
	//t(2,4) = t(4,2) = 2.0;

	ConcurrencyType concurrency(argc,argv);
	EngineType engine(t,concurrency,dof,true);
	std::vector<size_t> ne(dof,atoi(argv[2])); // 8 up and 8 down
	HilbertStateType gs(engine.size(),ne[0]);

	RealType sum = 0;
	for (size_t i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";	
	size_t sigma = 0;
	for (size_t site = 0; site<n ; site++) {
		LibraryOperatorType myOp(engine,LibraryOperatorType::N,site,sigma);
		for (size_t site2=0; site2<n; site2++) {
			HilbertStateType phi = gs;
			myOp.applyTo(phi);
			LibraryOperatorType myOp2(engine,LibraryOperatorType::N,site2,sigma);
			myOp2.applyTo(phi);
			std::cout<<scalarProduct(gs,phi)<<" ";
		}
		std::cout<<"\n";
	}
}
