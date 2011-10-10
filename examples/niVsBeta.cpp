

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

typedef double RealType;
typedef RealType FieldType;
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

void doOneBeta(const EngineType& engine,
               std::vector<size_t>& ne,
               OpLibFactoryType& opLibFactory,
               size_t site,
               size_t sigma,
               RealType beta)
{
	RealType sum = 0;
	RealType denominator = 0;
	FreeFermions::Combinations combinations(engine.size(),ne[0]);
	size_t factorToMakeItSmaller = 1e6;
	for (size_t i = 0; i<combinations.size(); ++i) {
		OpDiagonalFactoryType opDiagonalFactory(engine);
		EtoTheBetaHType ebh(beta,engine,0);
		DiagonalOperatorType& eibOp = opDiagonalFactory(ebh);
		
		std::vector<size_t> vTmp(engine.size(),0);
		for (size_t j=0;j<combinations(i).size();++j) vTmp[combinations(i)[j]]=1;
		std::vector<std::vector<size_t> > occupations(1,vTmp);
		HilbertStateType thisState(engine,occupations);
		HilbertStateType phi = thisState;
		eibOp.applyTo(phi);
		denominator += scalarProduct(thisState,phi)*factorToMakeItSmaller;
		LibraryOperatorType& myOp2 = opLibFactory(LibraryOperatorType::N,site,sigma);
		myOp2.applyTo(phi);
		sum += scalarProduct(thisState,phi)*factorToMakeItSmaller;
	}
	std::cout<<beta<<" "<<sum<<" "<<denominator<<" "<<sum/denominator<<"\n";	
}

int main(int argc,char *argv[])
{
	int opt;
	size_t n =0;
	size_t electronsUp=0;
	RealType step = 0;
	RealType offset=0;
	size_t total=0;
	size_t site = 0;
	bool ladder = false;
	std::vector<RealType> v;
	std::vector<std::string> str;

	while ((opt = getopt(argc, argv, "n:e:b:s:p:t:o:i:l")) != -1) {
		switch (opt) {
			case 'n':
				n = atoi(optarg);
				v.resize(n,0);
				break;
			case 'e':
				electronsUp = atoi(optarg);
				break;
			case 's':
				site = atoi(optarg);
				break;
			case 'p':
				if (v.size()==0) {
					std::string s(argv[0]);
					s += "-p must come after -n\n";
					throw std::runtime_error(s.c_str());
				}
				PsimagLite::tokenizer(optarg,str,",");
				if (str.size() & 1) {
					std::string s(argv[0]);
					s += " Expecting pairs for -p\n";
					throw std::runtime_error(s.c_str());
				}
				for (size_t i=0;i<str.size();i+=2) {
					v[atoi(str[i].c_str())] = atof(str[i+1].c_str());
				}
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
			case 'l':
				ladder = true;
				break;
			default: /* '?' */
				throw std::runtime_error("Wrong usage\n");
		}
	}
	if (n==0 || total==0) throw std::runtime_error("Wrong usage\n");

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
	geometry->addPotential(t,v);
	delete geometry;
	std::cerr<<t;

	ConcurrencyType concurrency(argc,argv);
	EngineType engine(t,concurrency,dof,true);

	std::vector<size_t> ne(dof,electronsUp); // 8 up and 8 down
	bool debug = false;
	HilbertStateType gs(engine,ne,debug);

	size_t sigma =0;

	std::cout<<"#site="<<site<<"\n";

	OpLibFactoryType opLibFactory(engine);
	for (size_t i=0;i<total;++i) {
		RealType beta = i*step + offset;
		doOneBeta(engine,ne,opLibFactory,site,sigma,beta);
	}
}
