

// SAmple of how to use FreeFermions core engine to calculate
// <c^\dagger_i c_j >
#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"
#include "Tokenizer.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef ComplexType FieldType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::Engine<RealType,FieldType,ConcurrencyType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::GeometryLibrary<MatrixType> GeometryLibraryType;
typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;
typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
typedef typename DiagonalOperatorType::FactoryType OpDiagonalFactoryType;

int main(int argc,char *argv[])
{
	int opt;
	size_t n =0,electronsUp=0,total=0;
	RealType offset = 0;
	std::vector<size_t> sites;
	std::vector<std::string> str;
	bool ladder = false;
	RealType step = 0;
	while ((opt = getopt(argc, argv, "n:e:b:s:t:o:i:l")) != -1) {
		switch (opt) {
		case 'n':
			n = atoi(optarg);
			break;
		case 'e':
			electronsUp = atoi(optarg);
			break;
		case 's':
			PsimagLite::tokenizer(optarg,str,",");
			for (size_t i=0;i<str.size();i++)
				sites.push_back(atoi(str[i].c_str()));
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
	

	/*GeometryLibraryType geometry

	*/
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
	HilbertStateType gs(engine.size(),ne,debug);

	size_t sigma =0;
	OperatorType myOp(engine,OperatorType::DESTRUCTION,sites[0],sigma);
	HilbertStateType phi2 = gs;
	myOp.applyTo(phi2);
	
//	FieldType density = scalarProduct(phi2,phi2);
//	std::cerr<<"density="<<density<<"\n";
	
	std::cout<<"#site="<<sites[0]<<"\n";
	std::cout<<"#site2="<<sites[1]<<"\n";
	for (size_t it = 0; it<total; it++) {
		OpDiagonalFactoryType opDiagonalFactory;
		RealType time = it * step + offset;
		EtoTheIhTimeType eih(time,engine,0);
		DiagonalOperatorType* eihOp = opDiagonalFactory(eih);
		HilbertStateType phi = gs;
		myOp.applyTo(phi);
		eihOp->applyTo(phi);
		OperatorType myOp2(engine,OperatorType::DESTRUCTION,sites[1],sigma);
		myOp2.applyTo(phi);
		std::cout<<time<<" "<<scalarProduct(phi,phi)<<"\n";
	}
}
