// SAmple of how to use FreeFermions core engine to calculate
// <c^\dagger_i c_j >
#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "EtoTheBetaH.h"
#include "DiagonalOperator.h"
#include "LibraryOperator.h"
#include "Combinations.h"
#include "GeometryParameters.h"
#include "Concurrency.h"
#include "InputNg.h"
#include "InputCheck.h"

typedef double RealType;
typedef RealType FieldType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
typedef FreeFermions::GeometryParameters<RealType,InputNgType::Readable> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::EToTheBetaH<EngineType> EtoTheBetaHType;
typedef FreeFermions::DiagonalOperator<EtoTheBetaHType> DiagonalOperatorType;
typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
typedef DiagonalOperatorType::FactoryType OpDiagonalFactoryType;
typedef OperatorType::FactoryType OpNormalFactoryType;
typedef FreeFermions::LibraryOperator<OperatorType> LibraryOperatorType;
typedef LibraryOperatorType::FactoryType OpLibFactoryType;

void doOneBeta(const EngineType& engine,
               SizeType site,
               RealType beta)
{
	RealType density = 0;
	for (SizeType i = 0; i<engine.size(); ++i) {
		RealType fermiFactor = 1.0/(1.0 + exp(beta*engine.eigenvalue(i)));
		density +=  fermiFactor;
	}

	RealType energy = 0;
	for (SizeType i = 0; i<engine.size(); ++i) {
		RealType fermiFactor = engine.eigenvalue(i)/(1.0 + exp(beta*engine.eigenvalue(i)));
		energy +=  fermiFactor;
	}

	RealType sum = 0;
	for (SizeType i = 0; i<engine.size(); ++i) {
		RealType fermiFactor = 1.0/(1.0 + exp(beta*engine.eigenvalue(i)));
		sum += PsimagLite::conj(engine.eigenvector(site,i)) * engine.eigenvector(site,i)*
		        fermiFactor;
	}

	std::cout<<beta<<" "<<sum<<" "<<(density/engine.size())<<" "<<energy<<"\n";
}

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	RealType step = 0;
	RealType offset=0;
	SizeType total=0;
	SizeType site = 0;
	RealType mu = 0;

	while ((opt = getopt(argc, argv, "f:s:t:o:i:m:")) != -1) {
		switch (opt) {
		case 'f':
			file=optarg;
			break;
		case 's':
			site = atoi(optarg);
			break;
		case 't':
			total = atoi(optarg);
			break;
		case 'i':
			step = atof(optarg);
			break;
		case 'o':
			offset = atof(optarg);
			break;
		case 'm':
			mu = atof(optarg);
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

	SizeType dof = 1; // spinless
	GeometryLibraryType geometry(geometryParams);

	geometry.addPotential(mu);

	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc,&argv,npthreads);
	EngineType engine(geometry.matrix(),
	                  geometryParams.outputFile,
	                  dof,
	                  EngineType::VerboseEnum::YES);
	std::cout<<geometry;
	std::cout<<"#site="<<site<<"\n";

	for (SizeType i=0;i<total;++i) {
		RealType beta = i*step + offset;
		doOneBeta(engine,site,beta);
	}
}
