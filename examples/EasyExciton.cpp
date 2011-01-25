
// Calculates
// <c^dagger_i exp(iHt) n_j exp(-iHt) c_i>

#include "Engine.h"
#include "ObservableLibrary.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"

typedef double RealType;
typedef std::complex<double> FieldType;
typedef std::vector<bool> LevelsType;
typedef PsimagLite::ConcurrencySerial<FieldType> ConcurrencyType;
typedef psimag::Matrix<FieldType> MatrixType;
typedef FreeFermions::Engine<RealType,FieldType,LevelsType,ConcurrencyType> EngineType;
typedef EngineType::HilbertVectorType HilbertVectorType;
typedef EngineType::FreeOperatorType FreeOperatorType;
typedef FreeFermions::GeometryLibrary<MatrixType> GeometryLibraryType;
typedef FreeFermions::ObservableLibrary<EngineType> ObservableLibraryType;
typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;

int main(int argc,char *argv[])
{
	if (argc!=5) throw std::runtime_error("Needs 5 arguments\n");
	size_t n = 16; 
	size_t electronsUp = 8;
	size_t dof = 2; // spin up and down
	MatrixType t(n,n);
	
	GeometryLibraryType geometry(n,GeometryLibraryType::CHAIN);
	geometry.setGeometry(t);
	//geometry.setGeometry(t,GeometryLibraryType::LADDER,2);
	std::cerr<<t;
	
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(t,concurrency,dof,true);

	std::vector<size_t> ne(dof,electronsUp); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	std::cout<<ne;
	size_t flavor =0;
	size_t site = atoi(argv[1]);
	FreeOperatorType myOp = engine.newSimpleOperator("destruction",site,flavor);
	HilbertVectorType phi = engine.newState();
	myOp.apply(phi,gs,FreeOperatorType::DO_NOT_SIMPLIFY);
	
	FieldType density = scalarProduct(phi,phi);
	std::cerr<<"density="<<density<<"\n";	
	
	size_t site2=atoi(argv[2]);
	std::cout<<"#site="<<site<<"\n";
	std::cout<<"#site2="<<site2<<"\n";
	for (size_t it = 0; it<size_t(atoi(argv[3])); it++) {
		RealType time = it * atof(argv[4]);
		EtoTheIhTimeType eih(time,engine);
		DiagonalOperatorType eihOp(eih);
		HilbertVectorType phi2 = engine.newState();
		eihOp.apply(phi2,phi);
		FreeOperatorType myOp2 = engine.newSimpleOperator("destruction",site2,flavor);
		HilbertVectorType phi3 = engine.newState();
		myOp2.apply(phi3,phi2,FreeOperatorType::DO_NOT_SIMPLIFY);
		phi2.clear();
		std::cout<<time<<" "<<real(scalarProduct(phi3,phi3))<<"\n";
	}
}
