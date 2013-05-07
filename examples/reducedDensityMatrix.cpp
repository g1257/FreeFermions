
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
typedef FreeFermions::GeometryParameters<RealType> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType,ConcurrencyType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
//typedef FreeFermions::HilbertState<OperatorType> HilbertStateType;
typedef FreeFermions::RealSpaceState<OperatorType> HilbertStateType;
typedef OperatorType::FactoryType OpNormalFactoryType;


int main(int argc,char* argv[])
{
	if (argc<3) throw std::runtime_error("Needs 2 or more argument(s)\n");
	size_t n = atoi(argv[1]);
	size_t nup =  atoi(argv[2]); 
	typedef FreeFermions::ReducedDensityMatrix<EngineType> ReducedDensityMatrixType;
	ConcurrencyType concurrency(argc,argv);
	
	GeometryParamsType geometryParams;
	geometryParams.type = GeometryLibraryType::CHAIN;
	geometryParams.sites = 2*n;
	if (argc>=4) {
		geometryParams.leg = atoi(argv[3]);
		if (geometryParams.leg>0) 
			  geometryParams.type = GeometryLibraryType::LADDER;
	}

	GeometryLibraryType geometry(geometryParams);
	if (concurrency.root()) std::cerr<<geometry;
	
	size_t dof = 1;
	EngineType engine(geometry,concurrency,dof,false);
	
	ReducedDensityMatrixType reducedDensityMatrix(engine,n,nup);
	
	typename PsimagLite::Vector<double>::Type e(reducedDensityMatrix.rank());
	reducedDensityMatrix.diagonalize(e);
	if (concurrency.root()) {
		std::cout<<"DensityMatrixEigenvalues:\n";
		std::cout<<e;
	}
}
