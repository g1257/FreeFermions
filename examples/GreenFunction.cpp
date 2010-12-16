
// Calculates
// <c^dagger_i exp(iHt) c_i>

#include "Engine.h"
#include "ObservableLibrary.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"

typedef double RealType;
typedef std::complex<double> FieldType;
typedef std::vector<bool> LevelsType;
typedef Dmrg::ConcurrencySerial<FieldType> ConcurrencyType;
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
	size_t nargs = 7;
	std::string s = "Needs " + utils::ttos(nargs) + " arguments\n";
	if (argc!=7) throw std::runtime_error(s.c_str());
	size_t n = atoi(argv[1]); 
	size_t electronsUp = atoi(argv[2]);
	size_t dof = 2; // spin up and down
	MatrixType t(n,n);
	
	GeometryLibraryType geometry(n,GeometryLibraryType::LADDER);
	geometry.setGeometry(t,2);
	size_t nb = 2; // basis sites per site
	RealType tb = 0.5; // hopping for basis sites
	geometry.bathify(t,nb,tb);
	//geometry.setGeometry(t,GeometryLibraryType::LADDER,2);
	std::cerr<<t;
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(t,concurrency,dof,true);
	std::vector<size_t> ne(dof,electronsUp); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	std::cout<<ne;

	size_t flavor =0;

	size_t site = atoi(argv[3]);
	FreeOperatorType myOp = engine.newSimpleOperator("destruction",site,flavor);
	HilbertVectorType phi = engine.newState();
	myOp.apply(phi,gs,FreeOperatorType::DO_NOT_SIMPLIFY);
	
	size_t site2=atoi(argv[4]);
	FieldType density = scalarProduct(phi,phi);
	std::cerr<<"density="<<density<<"\n";
	FreeOperatorType myOp2 = engine.newSimpleOperator("destruction",site2,flavor);
	HilbertVectorType phi2 = engine.newState();
	myOp.apply(phi2,gs,FreeOperatorType::DO_NOT_SIMPLIFY);
	
	std::cout<<"#site="<<site<<"\n";
	std::cout<<"#site2="<<site2<<"\n";
	for (size_t it = 0; it<size_t(atoi(argv[5])); it++) {
		RealType time = it * atof(argv[6]);
		EtoTheIhTimeType eih(time,engine);
		DiagonalOperatorType eihOp(eih);
		HilbertVectorType phi3 = engine.newState();
		eihOp.apply(phi3,phi);
		std::cout<<time<<" "<<real(scalarProduct(phi2,phi3))<<"\n";
	}
}
