
// Calculates
// <c^dagger_i exp(iHt) n_j exp(-iHt) c_i>

#include "Engine.h"
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"

using namespace FreeFermions;

typedef double RealType;
typedef std::complex<double> FieldType;
typedef size_t UnsignedIntegerType;
typedef Engine<RealType,FieldType,UnsignedIntegerType> EngineType;
typedef EngineType::HilbertVectorType HilbertVectorType;
typedef EngineType::FreeOperatorType FreeOperatorType;
typedef EToTheIhTime<EngineType> EtoTheIhTimeType;
typedef DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;


int main(int argc,char *argv[])
{
	if (argc!=5) throw std::runtime_error("Needs 5 arguments\n");
	size_t n = 16; // 16 sites
	size_t dof = 2; // spin up and down
	bool isPeriodic = false;
	psimag::Matrix<FieldType> t(n,n);
	for (size_t i=0;i<n;i++) {
		for (size_t j=0;j<n;j++) {
			if (i-j==1 || j-i==1) t(i,j) = 1.0;
		}
	}
	if (isPeriodic) t(0,n-1) = t(n-1,0) = 1.0;
	
	EngineType engine(t,dof,false);
	std::vector<size_t> ne(dof,8); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	std::cout<<ne;
	size_t flavor =0;
	size_t site = atoi(argv[1]);
	FreeOperatorType myOp = engine.newSimpleOperator("destruction",site,flavor);
	HilbertVectorType phi = engine.newState();
	myOp.apply(phi,gs);
	
	FieldType density = scalarProduct(phi,phi);
	std::cerr<<"density="<<density<<"\n";	
	
	size_t site2=atoi(argv[2]);
		
	for (size_t it = 0; it<size_t(atoi(argv[3])); it++) {
		RealType time = it * atof(argv[4]);
		EtoTheIhTimeType eih(time,engine);
		DiagonalOperatorType eihOp(eih);
		HilbertVectorType phi2 = engine.newState();
		eihOp.apply(phi2,phi);
	
		FreeOperatorType myOp2 = engine.newSimpleOperator("destruction",site2,flavor);
		HilbertVectorType phi3 = engine.newState();
		myOp2.apply(phi3,phi2);
		
		std::cout<<site2<<" "<<scalarProduct(phi3,phi3)/density<<"\n";
	}
}
