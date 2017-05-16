#define USE_PTHREADS_OR_NOT_NG
#include "ReducedDensityMatrix.h"
#include "Concurrency.h"

// TBW FIXME
typedef double RealType;
typedef std::complex<double> ComplexType;
typedef RealType FieldType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
typedef FreeFermions::GeometryParameters<RealType,InputNgType::Readable> GeometryParamsType;
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

	FreeFermions::InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(file,inputCheck);
	InputNgType::Readable io(ioWriteable);

	GeometryParamsType geometryParams(io);
	SizeType electronsUp = GeometryParamsType::readElectrons(io,geometryParams.sites);

	SizeType dof = 1; // spinless

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;

	SizeType npthreads = 1;
	try {
		io.readline(npthreads,"Threads=");
	} catch (std::exception& e) {}

	PsimagLite::Concurrency concurrency(&argc,&argv,npthreads);
	EngineType engine(geometry.matrix(),dof,true);
	PsimagLite::Vector<SizeType>::Type ne(dof,electronsUp); // n. of up (= n. of  down electrons)
	HilbertStateType gs(engine,ne,0,false);
	RealType sum = 0;
	for (SizeType i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";

	SizeType halfSites = static_cast<SizeType>(0.5*engine.size());	
	ReducedDensityMatrixType reducedDensityMatrix(engine,halfSites,electronsUp);
	
	PsimagLite::Vector<double>::Type e(reducedDensityMatrix.rank());
	reducedDensityMatrix.diagonalize(e);
	if (concurrency.root()) {
		std::cout<<"DensityMatrixEigenvalues:\n";
		std::cout<<e;
	}
}

