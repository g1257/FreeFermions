
#include "ReducedDensityMatrix.h"
#include "Concurrency.h"

// TBW FIXME
typedef double RealType;
typedef std::complex<double> ComplexType;
typedef RealType FieldType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::GeometryParameters<RealType> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
//typedef FreeFermions::HilbertState<OperatorType> HilbertStateType;
typedef FreeFermions::RealSpaceState<OperatorType> HilbertStateType;
typedef OperatorType::FactoryType OpNormalFactoryType;
typedef FreeFermions::ReducedDensityMatrix<EngineType> ReducedDensityMatrixType;


int main(int argc,char* argv[])
{
	int opt = 0;
	PsimagLite::String file("");

	while ((opt = getopt(argc, argv, "f:")) != -1) {
		switch (opt) {
		case 'f':
			file=optarg;
			break;
		default: /* '?' */
			throw std::runtime_error("Wrong usage\n");
		}
	}

	if (file=="") {
		throw std::runtime_error("Wrong usage\n");
	}

	GeometryParamsType geometryParams(file);
	size_t electronsUp = GeometryParamsType::readElectrons(file,geometryParams.sites);

	size_t dof = 1; // spinless

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;
	SizeType npthreads = 1;
	PsimagLite::Concurrency concurrency(&argc,&argv,npthreads);
	EngineType engine(geometry.matrix(),dof,true);
	PsimagLite::Vector<size_t>::Type ne(dof,electronsUp); // n. of up (= n. of  down electrons)
	HilbertStateType gs(engine,ne);
	RealType sum = 0;
	for (size_t i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";
	
	ReducedDensityMatrixType reducedDensityMatrix(engine,geometryParams.sites,electronsUp);
	
	PsimagLite::Vector<double>::Type e(reducedDensityMatrix.rank());
	reducedDensityMatrix.diagonalize(e);
	if (concurrency.root()) {
		std::cout<<"DensityMatrixEigenvalues:\n";
		std::cout<<e;
	}
}
