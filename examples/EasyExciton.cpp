

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
	if (argc!=6) throw std::runtime_error("Needs 6 arguments\n");
	size_t n = 50;
	size_t electronsUp = 25;
	size_t dof = 1; // spinless
	MatrixType t(n,n);
	
	/*GeometryLibraryType geometry(n,GeometryLibraryType::LADDER);
	geometry.setGeometry(t,2);
	for (size_t ii=0;ii<n;ii+=2) 
		t(ii,ii+1) = t(ii+1,ii) = 0.5;
	*/
	GeometryLibraryType geometry(n,GeometryLibraryType::CHAIN);
	geometry.setGeometry(t);
	std::cerr<<t;
	
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(t,concurrency,dof,true);

	std::vector<size_t> ne(dof,electronsUp); // 8 up and 8 down
	bool debug = false;
	HilbertStateType gs(engine.size(),ne,debug);

	size_t sigma =0;
	size_t site = atoi(argv[1]);
	OperatorType myOp(engine,OperatorType::DESTRUCTION,site,sigma);
	HilbertStateType phi2 = gs;
	myOp.applyTo(phi2);
	
//	FieldType density = scalarProduct(phi2,phi2);
//	std::cerr<<"density="<<density<<"\n";
	
	size_t site2=atoi(argv[2]);
	std::cout<<"#site="<<site<<"\n";
	std::cout<<"#site2="<<site2<<"\n";
	for (size_t it = 0; it<size_t(atoi(argv[3])); it++) {
		OpDiagonalFactoryType opDiagonalFactory;
		RealType time = it * atof(argv[4]) + atof(argv[5]);
		EtoTheIhTimeType eih(time,engine,0);
		DiagonalOperatorType* eihOp = opDiagonalFactory(eih);
		HilbertStateType phi = gs;
		myOp.applyTo(phi);
		eihOp->applyTo(phi);
		OperatorType myOp2(engine,OperatorType::DESTRUCTION,site2,sigma);
		myOp2.applyTo(phi);
		std::cout<<time<<" "<<scalarProduct(phi,phi)<<"\n";
	}
}
