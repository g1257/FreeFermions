

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

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef ComplexType FieldType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::Engine<RealType,FieldType,ConcurrencyType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::GeometryLibrary<MatrixType> GeometryLibraryType;
typedef FreeFermions::EToTheBetaH<EngineType> EtoTheBetaHType;
typedef FreeFermions::DiagonalOperator<EtoTheBetaHType> DiagonalOperatorType;
typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
typedef DiagonalOperatorType::FactoryType OpDiagonalFactoryType;
typedef OperatorType::FactoryType OpNormalFactoryType;
typedef FreeFermions::LibraryOperator<OperatorType> LibraryOperatorType;
typedef LibraryOperatorType::FactoryType OpLibFactoryType;

int main(int argc,char *argv[])
{
	int opt;
	size_t n =0,electronsUp=0;
	RealType beta = 0;
	size_t site = 0;
	bool ladder = false;
	RealType step = 0;
	while ((opt = getopt(argc, argv, "n:e:b:s:l")) != -1) {
		switch (opt) {
		case 'n':
			n = atoi(optarg);
			break;
		case 'e':
			electronsUp = atoi(optarg);
			break;
		case 's':
			site = atoi(optarg);
			break;
		case 'b':
			beta = atoi(optarg);
			break;
		case 'l':
			ladder = true;
			break;
		default: /* '?' */
			throw std::runtime_error("Wrong usage\n");
		}
	}
	if (n==0) throw std::runtime_error("Wrong usage\n");

	size_t dof = 1; // spinless
	MatrixType t(n,n);

	GeometryLibraryType* geometry;
	if (!ladder) {
		geometry = new GeometryLibraryType(n,GeometryLibraryType::CHAIN);
		geometry->setGeometry(t);
	} else {
		geometry = new GeometryLibraryType(n,GeometryLibraryType::LADDER);
		geometry->setGeometry(t,2);
		for (size_t ii=0;ii<n;ii+=2)
			t(ii,ii+1) = t(ii+1,ii) = 0.5;
	}
	delete geometry;
	std::cerr<<t;
	
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(t,concurrency,dof,true);

	std::vector<size_t> ne(dof,electronsUp); // 8 up and 8 down
	bool debug = false;
	HilbertStateType gs(engine,ne,debug);

	size_t sigma =0;
	OpNormalFactoryType opNormalFactory(engine);
	OperatorType& myOp = opNormalFactory(OperatorType::DESTRUCTION,sites[0],sigma);
	HilbertStateType phi2 = gs;
	myOp.applyTo(phi2);
	
//	FieldType density = scalarProduct(phi2,phi2);
//	std::cerr<<"density="<<density<<"\n";
	
	std::cout<<"#site="<<site<<"\n";
	std::cout<<"#beta="<<beta<<"\n";
	
	OpLibFactoryType opLibFactory(engine);
	
	RealType sum = 0;
	RealType denominator = 0;
	for (size_t i = 0; i<basis.size(); it++) {
		OpDiagonalFactoryType opDiagonalFactory(engine);
		EtoTheBetaHType ebh(beta,engine,0);
		DiagonalOperatorType& eibOp = opDiagonalFactory(ebh);
		HilbertStateType phi = basis.state(i);
		eibOp.applyTo(phi);
		denominator += scalarProduct(basis.state(i),phi);
		LibraryOperatorType& myOp2 = opLibFactory(LibraryOperatorType::N,site,sigma);
		myOp2.applyTo(phi);
		sum += scalarProduct(basis.state(i),phi);
	}
	std::cout<<sum<<" "<<denominator<<" "<<sum/denominator<<"\n";
}
