

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
	bool ladder = false;
	std::vector<RealType> v;
	std::vector<std::string> str;

	while ((opt = getopt(argc, argv, "n:s:p:t:o:i:l")) != -1) {
		switch (opt) {
			case 'n':
				n = atoi(optarg);
				v.resize(n,0);
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
	GeometryParamsType geometryParams;
	geometryParams.sites = n;
	GeometryLibraryType* geometry;
	if (!ladder) {
		geometryParams.type = GeometryLibraryType::CHAIN;
		geometry = new GeometryLibraryType(geometryParams);
	} else {
		geometryParams.hopping.resize(2);
		geometryParams.hopping[0] = 1.0;
		geometryParams.hopping[1] = 0.5;
		geometryParams.type = GeometryLibraryType::LADDER;
		geometry = new GeometryLibraryType(geometryParams);
	}
	geometry->addPotential(v);
	
	std::cerr<<v;

	ConcurrencyType concurrency(argc,argv);
	EngineType engine(*geometry,concurrency,dof,true);

	std::cout<<(*geometry);
	std::cout<<"#site="<<site<<"\n";

	for (size_t i=0;i<total;++i) {
		RealType beta = i*step + offset;
		doOneBeta(engine,site,beta);
	}
	delete geometry;
}
