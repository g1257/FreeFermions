

// Calculates <phi | phi>
// where 
// |phi> = c_{p up} exp(iHt) c_{i\sigma} nbar_{i \bar{sigma}} c^dagger_{j sigma'} n_{j \bar{sigma'}} |gs>

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
typedef PsimagLite::Matrix<FieldType> MatrixType;
typedef FreeFermions::Engine<RealType,FieldType,LevelsType,ConcurrencyType> EngineType;
typedef EngineType::HilbertVectorType HilbertVectorType;
typedef EngineType::FreeOperatorType FreeOperatorType;
typedef FreeFermions::GeometryLibrary<MatrixType> GeometryLibraryType;
typedef FreeFermions::ObservableLibrary<EngineType> ObservableLibraryType;
typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;

FieldType calcSuperDensity(size_t site, size_t site2,HilbertVectorType gs,const EngineType& engine,const ObservableLibraryType& library)
{
	size_t DO_NOT_SIMPLIFY = FreeOperatorType::DO_NOT_SIMPLIFY;
	FieldType savedValue = 0;
	
	HilbertVectorType tmpV = engine.newState();
	bool verbose = false;
	std::vector<HilbertVectorType> savedVector(4,engine.newState(verbose));
	
	for (size_t sigma = 0;sigma<2;sigma++) {
		HilbertVectorType phi = engine.newState();
		library.applyNiOneFlavor(phi,gs,site,1-sigma);
		
		FreeOperatorType myOp2 = engine.newSimpleOperator("creation",site,sigma);
		HilbertVectorType phi2 = engine.newState();
		myOp2.apply(phi2,phi,DO_NOT_SIMPLIFY);
		
		tmpV.add(phi2);
		phi.clear();

		for (size_t sigma2 = 0;sigma2 < 2;sigma2++) {
			HilbertVectorType phi3 = engine.newState();
			library.applyNiBarOneFlavor(phi3,phi2,site2,1-sigma2);
			
			FreeOperatorType myOp4 = engine.newSimpleOperator("destruction",site2,sigma2);
			myOp4.apply(savedVector[sigma+sigma2*2],phi3,DO_NOT_SIMPLIFY);
			phi3.clear();
			
		}
	}
	FieldType sum = 0;
	for (size_t sigma=0;sigma<4;sigma++)
		for (size_t sigma2=0;sigma2<4;sigma2++)
			sum += scalarProduct(savedVector[sigma],savedVector[sigma2]);
	std::cerr<<"#sum2="<<scalarProduct(tmpV,tmpV)<<"\n";
	return sum;
}


int main(int argc,char *argv[])
{
	if (argc!=7) throw std::runtime_error("Needs 7 arguments\n");
	size_t n = 12; 
	size_t electronsUp = 6;
	size_t dof = 2; // spin up and down
	
	MatrixType t(n,n);
	GeometryLibraryType geometry(n,GeometryLibraryType::CHAIN);
	geometry.setGeometry(t);
	//geometry.setGeometry(t,GeometryLibraryType::LADDER,2);
	std::cerr<<t;
	
	bool verbose = true;
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(t,concurrency,dof,false);
	ObservableLibraryType library(engine);
	
	std::vector<size_t> ne(dof,electronsUp); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	
	size_t site = atoi(argv[1]);
	size_t site2 = atoi(argv[2]);
	size_t site3 = atoi(argv[3]);
	//size_t sigma3 = 0;
	size_t DO_NOT_SIMPLIFY = FreeOperatorType::DO_NOT_SIMPLIFY;
	
	FieldType superdensity = calcSuperDensity(site,site2,gs,engine,library);
	std::cout<<"#superdensity="<<superdensity<<"\n";
	std::cout<<"#site="<<site<<" site2="<<site2<<"\n";
	
	for (size_t it=0;it<size_t(atoi(argv[4]));it++) {
		RealType time = it * atof(argv[5]) + atof(argv[6]);
		EtoTheIhTimeType eih(time,engine);
		DiagonalOperatorType eihOp(eih);

		std::vector<HilbertVectorType> savedVector(4,engine.newState(verbose));
		FieldType savedValue = 0;
		
		for (size_t sigma = 0;sigma<2;sigma++) {
			HilbertVectorType phi = engine.newState();
			library.applyNiOneFlavor(phi,gs,site,1-sigma);
	
			FreeOperatorType myOp2 = engine.newSimpleOperator("creation",site,sigma);
			HilbertVectorType phi2 = engine.newState();
			myOp2.apply(phi2,phi,DO_NOT_SIMPLIFY);
			phi.clear();
			
			for (size_t sigma2 = 0;sigma2 < 2;sigma2++) {
				HilbertVectorType phi3 = engine.newState();
				library.applyNiBarOneFlavor(phi3,phi2,site2,1-sigma2);
				
				FreeOperatorType myOp4 = engine.newSimpleOperator("destruction",site2,sigma2);
				HilbertVectorType phi4 = engine.newState(verbose);
				myOp4.apply(phi4,phi3,DO_NOT_SIMPLIFY);
				phi3.clear();
				
				if (verbose) std::cerr<<"Applying exp(iHt)\n";
				HilbertVectorType phi5 = engine.newState(verbose);
				eihOp.apply(phi5,phi4);
				phi4.clear();
				
				if (verbose) std::cerr<<"Applying c_{p down}\n";
				FreeOperatorType myOp6 = engine.newSimpleOperator("destruction",site3,1);
				HilbertVectorType phi6 = engine.newState(verbose);
				myOp6.apply(phi6,phi5,DO_NOT_SIMPLIFY);
				phi5.clear();
				
				if (verbose) std::cerr<<"Applying c_{p up}\n";
				FreeOperatorType myOp7 = engine.newSimpleOperator("destruction",site3,0);
				myOp7.apply(savedVector[sigma+sigma2*2],phi6,DO_NOT_SIMPLIFY);
				phi6.clear();
			}
		}
		FieldType sum = 0;
		for (size_t sigma=0;sigma<4;sigma++)
			for (size_t sigma2=0;sigma2<4;sigma2++)
				sum += scalarProduct(savedVector[sigma],savedVector[sigma2]);
			
		std::cout<<time<<" "<<sum<<"\n";
	}
}
