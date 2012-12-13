

// SAmple of how to use FreeFermions core engine to calculate
// <c^\dagger_i c_j >
#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "EtoTheBetaH.h"
#include "DiagonalOperator.h"
#include "Tokenizer.h"
#include "LibraryOperator.h"
#include "Combinations.h"
#include "GeometryParameters.h"
#include "DriverHelper.h"

typedef double RealType;
typedef RealType FieldType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::GeometryParameters<RealType> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType,ConcurrencyType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::EToTheBetaH<EngineType> EtoTheBetaHType;
typedef FreeFermions::DiagonalOperator<EtoTheBetaHType> DiagonalOperatorType;
typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
typedef DiagonalOperatorType::FactoryType OpDiagonalFactoryType;
typedef OperatorType::FactoryType OpNormalFactoryType;
typedef FreeFermions::LibraryOperator<OperatorType> LibraryOperatorType;
typedef LibraryOperatorType::FactoryType OpLibFactoryType;
typedef FreeFermions::DriverHelper<GeometryLibraryType> DriverHelperType;

void doOneBeta(const EngineType& engine,
               size_t site,
               RealType beta)
{
	RealType density = 0;
	for (size_t i = 0; i<engine.size(); ++i) {
		RealType fermiFactor = 1.0/(1.0 + exp(beta*engine.eigenvalue(i)));
		density +=  fermiFactor;
	}

	RealType energy = 0;
	for (size_t i = 0; i<engine.size(); ++i) {
		RealType fermiFactor = engine.eigenvalue(i)/(1.0 + exp(beta*engine.eigenvalue(i)));
		energy +=  fermiFactor;
	}

	RealType sum = 0;
	for (size_t i = 0; i<engine.size(); ++i) {
		RealType fermiFactor = 1.0/(1.0 + exp(beta*engine.eigenvalue(i)));
		sum += std::conj(engine.eigenvector(site,i)) * engine.eigenvector(site,i)* fermiFactor;
	}

	std::cout<<beta<<" "<<sum<<" "<<density<<" "<<energy<<"\n";
}

int main(int argc,char *argv[])
{
	int opt;
	size_t n =0;
	RealType step = 0;
	RealType offset=0;
	size_t total=0;
	size_t site = 0;
	std::vector<RealType> v;
	std::vector<std::string> str;
	GeometryParamsType geometryParams;

	while ((opt = getopt(argc, argv, "n:s:p:t:o:i:g")) != -1) {
		switch (opt) {
			case 'n':
				n = atoi(optarg);
				geometryParams.sites = n;
				break;
			case 's':
				site = atoi(optarg);
				break;
			case 'p':
				DriverHelperType::readPotential(v,optarg);
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
			case 'g':
				str.clear();
				PsimagLite::tokenizer(optarg,str,",");
				DriverHelperType::setMyGeometry(geometryParams,str);
				break;
			default: /* '?' */
				DriverHelperType::usage(argv[0],"-n sites -g geometry,[leg,filename]");
				throw std::runtime_error("Wrong usage\n");
		}
	}

	if (n==0 || total==0) {
		DriverHelperType::usage(argv[0],"-n sites -g geometry,[leg,filename]");
		throw std::runtime_error("Wrong usage\n");
	}

	size_t dof = 1; // spinless
	GeometryLibraryType geometry(geometryParams);

	geometry.addPotential(v);
	
	std::cerr<<v;

	ConcurrencyType concurrency(argc,argv);
	EngineType engine(geometry,concurrency,dof,true);

	std::cout<<geometry;
	std::cout<<"#site="<<site<<"\n";

	for (size_t i=0;i<total;++i) {
		RealType beta = i*step + offset;
		doOneBeta(engine,site,beta);
	}
}
