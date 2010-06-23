

// Calculates <phi | phi>
// where 
// |phi> = c_{p up} exp(iHt) c_{i\sigma} nbar_{i \bar{sigma}} c^dagger_{j sigma'} n_{j \bar{sigma'}} |gs>

#include "Engine.h"
#include "ObservableLibrary.h"
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
typedef ObservableLibrary<EngineType> ObservableLibraryType;

int main(int argc,char *argv[])
{
	if (argc!=6) throw std::runtime_error("Needs 5 arguments\n");
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
	ObservableLibraryType library(engine);
	
	std::vector<size_t> ne(dof,8); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	
	size_t site = atoi(argv[1]);
	size_t site2 = atoi(argv[2]);
	size_t site3 = atoi(argv[3]);
	size_t sigma3 = 0;
	for (size_t it=0;it<size_t(atoi(argv[4]));it++) {
		RealType time = it * atof(argv[5]);
		HilbertVectorType timeVector = engine.newState();
		for (size_t sigma = 0;sigma<2;sigma++) {
			HilbertVectorType phi = engine.newState();
			library.applyNiOneFlavor(phi,gs,site,1-sigma);
	
			FreeOperatorType myOp2 = engine.newSimpleOperator("creation",site,sigma);
			HilbertVectorType phi2 = engine.newState();
			myOp2.apply(phi2,phi);
	
		
			for (size_t sigma2 = 0;sigma2 < 2;sigma2++) {
				HilbertVectorType phi3 = engine.newState();
				library.applyNiBarOneFlavor(phi,gs,site2,1-sigma2);
	
				FreeOperatorType myOp4 = engine.newSimpleOperator("destruction",site2,sigma2);
				HilbertVectorType phi4 = engine.newState();
				myOp4.apply(phi4,phi3);
	
				EtoTheIhTimeType eih(time,engine);
				DiagonalOperatorType eihOp(eih);
				HilbertVectorType phi5 = engine.newState();
				eihOp.apply(phi5,phi4);
	
		
		
				FreeOperatorType myOp6 = engine.newSimpleOperator("destruction",site3,sigma3);
				HilbertVectorType phi6 = engine.newState();
				myOp6.apply(phi6,phi5);
				timeVector.add(phi6);
			}
		}
		std::cout<<site3<<" "<<scalarProduct(timeVector,timeVector)<<"\n";
	}
}
