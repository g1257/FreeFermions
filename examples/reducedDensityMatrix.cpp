
#include "ReducedDensityMatrix.h"
#ifdef USE_MPI
#include "ConcurrencyMpi.h"
#else
#include "ConcurrencySerial.h"
#endif

// TBW FIXME
typedef double RealType;
typedef std::complex<double> ComplexType;
typedef RealType FieldType;
#ifndef USE_MPI
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
#else
typedef PsimagLite::ConcurrencyMpi<RealType> ConcurrencyType;
#endif

typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::Engine<RealType,FieldType,ConcurrencyType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
//typedef FreeFermions::HilbertState<OperatorType> HilbertStateType;
typedef FreeFermions::RealSpaceState<OperatorType> HilbertStateType;
typedef FreeFermions::GeometryLibrary<MatrixType> GeometryLibraryType;
typedef OperatorType::FactoryType OpNormalFactoryType;


int main(int argc,char* argv[])
{
	if (argc<3) throw std::runtime_error("Needs 2 or more argument(s)\n");
	size_t n = atoi(argv[1]);
	size_t nup =  atoi(argv[2]); 
	typedef FreeFermions::ReducedDensityMatrix<EngineType> ReducedDensityMatrixType;
	ConcurrencyType concurrency(argc,argv);
	MatrixType t(2*n,2*n);
	size_t geometryType = GeometryLibraryType::CHAIN;
	size_t leg = GeometryLibraryType::OPTION_NONE;
	if (argc>=4) {
		leg = atoi(argv[3]);
		if (leg>0) geometryType = GeometryLibraryType::LADDER;
	}
	GeometryLibraryType geometry(2*n,geometryType);
	geometry.setGeometry(t,leg);
	if (concurrency.root()) std::cerr<<t;
	
	size_t dof = 1;
	EngineType engine(t,concurrency,dof,false);
	
	ReducedDensityMatrixType reducedDensityMatrix(engine,n,nup);
	
	std::vector<double> e(reducedDensityMatrix.rank());
	reducedDensityMatrix.diagonalize(e);
	if (concurrency.root()) {
		std::cout<<"DensityMatrixEigenvalues:\n";
		std::cout<<e;
	}
}
