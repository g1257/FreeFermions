

// SAmple of how to use FreeFermions core engine to calculate
// <c^\dagger_i c_j >
#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"
#include "Tokenizer.h"
#include "GeometryParameters.h"
#include "Concurrency.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef ComplexType FieldType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::GeometryParameters<RealType> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;
typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
typedef DiagonalOperatorType::FactoryType OpDiagonalFactoryType;
typedef OperatorType::FactoryType OpNormalFactoryType;

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	size_t total=0;
	RealType offset = 0;
	RealType step = 0;
	size_t site3 = 0;

	while ((opt = getopt(argc, argv, "f:t:o:i:p:")) != -1) {
		switch (opt) {
		case 'f':
			file=optarg;
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
		case 'p':
			site3 = atoi(optarg);
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
	PsimagLite::Vector<size_t>::Type sites;
	GeometryParamsType::readVector(sites,file,"TSPSites");
	sites.resize(3);
	sites[2] = site3;

	size_t dof = 1; // spinless

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;
	
	SizeType npthreads = 1;
	ConcurrencyType concurrency(argc,argv,npthreads);
	EngineType engine(geometry,dof,true);

	PsimagLite::Vector<size_t>::Type ne(dof,electronsUp); // 8 up and 8 down
	bool debug = false;
	HilbertStateType gs(engine,ne,debug);

	size_t sigma =0;
	OpNormalFactoryType opNormalFactory(engine);
	OperatorType& myOp = opNormalFactory(OperatorType::DESTRUCTION,sites[0],sigma);
	HilbertStateType phi2 = gs;
	myOp.applyTo(phi2);
	
//	FieldType density = scalarProduct(phi2,phi2);
//	std::cerr<<"density="<<density<<"\n";
	
	std::cout<<"#site="<<sites[0]<<"\n";
	std::cout<<"#site2="<<sites[1]<<"\n";
	for (size_t it = 0; it<total; it++) {
		OpDiagonalFactoryType opDiagonalFactory(engine);
		RealType time = it * step + offset;
		EtoTheIhTimeType eih(time,engine,0);
		DiagonalOperatorType& eihOp = opDiagonalFactory(eih);
		HilbertStateType phi = gs;
		myOp.applyTo(phi);
		eihOp.applyTo(phi);
		OperatorType& myOp2 = opNormalFactory(OperatorType::DESTRUCTION,sites[1],sigma);
		myOp2.applyTo(phi);
		std::cout<<time<<" "<<scalarProduct(phi,phi)<<"\n";
	}
}
