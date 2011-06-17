

// Calculates <phi | phi>
// where 
// |phi> = c_{p up} exp(iHt) c_{i\sigma} nbar_{i \bar{sigma}} c^dagger_{j sigma'} n_{j \bar{sigma'}} |gs>

#include "Engine.h"
#include "ObservableLibrary.h"
#include "GeometryLibrary.h"
#ifdef USE_MPI
#include "ConcurrencyMpi.h"
#else
#include "ConcurrencySerial.h"
#endif
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"
#include "MemoryUsage.h" // in PsimagLite

typedef PsimagLite::MemoryUsage MemoryUsageType;

typedef float RealType;
typedef std::complex<RealType> FieldType;
typedef std::vector<bool> LevelsType;
#ifdef USE_MPI
typedef PsimagLite::ConcurrencyMpi<FieldType> ConcurrencyType;
#else
typedef PsimagLite::ConcurrencySerial<FieldType> ConcurrencyType;
#endif
typedef PsimagLite::Matrix<FieldType> MatrixType;
typedef FreeFermions::Engine<RealType,FieldType,LevelsType,ConcurrencyType> EngineType;
typedef EngineType::HilbertVectorType HilbertVectorType;
typedef EngineType::FreeOperatorType FreeOperatorType;
typedef FreeFermions::GeometryLibrary<MatrixType> GeometryLibraryType;
typedef FreeFermions::ObservableLibrary<EngineType> ObservableLibraryType;
typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;

FieldType calcSuperDensity(size_t site, size_t site2,const HilbertVectorType& gs,const EngineType& engine,const ObservableLibraryType& library)
{
	size_t DO_NOT_SIMPLIFY = FreeOperatorType::DO_NOT_SIMPLIFY;
	HilbertVectorType savedVector = engine.newState();
	FieldType savedValue = 0;
	FieldType sum = 0;
	HilbertVectorType tmpV = engine.newState();
	MemoryUsageType myUsage;

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
			std::cerr<<" sigma= "<<sigma<<" sigma2="<<sigma2<<" phi3.terms="<<phi3.terms()<<"\n";
			std::cerr<<"memory="<<myUsage.vmSize()<<"\n";
			std::cerr.flush();

			FreeOperatorType myOp4 = engine.newSimpleOperator("destruction",site2,sigma2);
			HilbertVectorType phi4 = engine.newState();
			myOp4.apply(phi4,phi3,DO_NOT_SIMPLIFY);
			phi3.clear();
			std::cerr<<" sigma= "<<sigma<<" sigma2="<<sigma2<<" phi4.terms="<<phi4.terms()<<"\n";
			std::cerr<<"memory="<<myUsage.vmSize()<<"\n";

			sum += scalarProduct(phi4,phi4);
			if (sigma ==0 && sigma2 ==0) savedVector = phi4;
			if (sigma ==1 && sigma2 ==1) {
				savedValue = scalarProduct(phi4,savedVector);
				savedVector.clear();
			}
			std::cerr<<" sigma= "<<sigma<<" sigma2="<<sigma2<<" saved.terms="<<savedVector.terms()<<"\n";
			std::cerr<<"memory="<<myUsage.vmSize()<<"\n";
		}
	}
	sum += 2*real(savedValue);
	std::cerr<<"#sum2="<<scalarProduct(tmpV,tmpV)<<"\n";
	return sum;
}


int main(int argc,char *argv[])
{
	if (argc!=7) throw std::runtime_error("Needs 7 arguments\n");
	size_t n = 32; 
	size_t electronsUp = 16;
	size_t dof = 2; // spin up and down
	
	MatrixType t(n,n);
	GeometryLibraryType geometry(n,GeometryLibraryType::LADDER);
	geometry.setGeometry(t,2);
	//geometry.setGeometry(t,GeometryLibraryType::LADDER,2);
	//
	for (size_t ii=0;ii<n;ii+=2)
		t(ii,ii+1) = t(ii+1,ii) = 0.5;

	std::cerr<<t;
	
	bool verbose = false;
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(t,concurrency,dof,false);
	ObservableLibraryType library(engine);
	
	std::vector<size_t> ne(dof,electronsUp); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	
	size_t site = atoi(argv[1]);
	size_t site2 = atoi(argv[2]);
	size_t site3 = atoi(argv[3]);
	size_t sigma3 = 0;
	//size_t DO_NOT_SIMPLIFY = FreeOperatorType::DO_NOT_SIMPLIFY;
	
	FieldType superdensity = calcSuperDensity(site,site2,gs,engine,library);
	std::cout<<"#superdensity="<<superdensity<<"\n";
	std::cout<<"#site="<<site<<" site2="<<site2<<"\n";	
	concurrency.loopCreate(size_t(atoi(argv[4])));
	size_t it = 0;
	while(concurrency.loop(it)) {
	//for (size_t it=0;it<size_t(atoi(argv[4]));it++) {
		RealType time = it * atof(argv[5]) + atof(argv[6]);
		EtoTheIhTimeType eih(time,engine,0);
		DiagonalOperatorType eihOp(eih);
				
		HilbertVectorType savedVector = engine.newState(verbose);
		FieldType savedValue = 0;
		FieldType sum = 0;
		MemoryUsageType myUsage;
		for (size_t sigma = 0;sigma<2;sigma++) {
			HilbertVectorType phi = engine.newState();
			library.applyNiOneFlavor(phi,gs,site,1-sigma);

			FreeOperatorType myOp2 = engine.newSimpleOperator("creation",site,sigma);
			HilbertVectorType phi2 = engine.newState();
			myOp2.apply(phi2,phi,FreeOperatorType::SIMPLIFY);
			phi.clear();
			std::cerr<<"time="<<time<<" sigma= "<<sigma<<" phi2.terms="<<phi2.terms()<<"\n";
			std::cerr<<"memory="<<myUsage.vmSize()<<"\n";
			std::cerr.flush();
			size_t simplifyOrNot = FreeOperatorType::SIMPLIFY; 
			for (size_t sigma2 = 0;sigma2 < 2;sigma2++) {
				HilbertVectorType phi3 = engine.newState();
				library.applyNiBarOneFlavor(phi3,phi2,site2,1-sigma2);
				
				FreeOperatorType myOp4 = engine.newSimpleOperator("destruction",site2,sigma2);
				HilbertVectorType phi4 = engine.newState(verbose);
				myOp4.apply(phi4,phi3,simplifyOrNot);
				phi3.clear();
				std::cerr<<"time="<<time<<" sigma= "<<sigma<<" sigma2="<<sigma2<<" phi4.terms="<<phi4.terms()<<"\n";
				std::cerr<<"memory="<<myUsage.vmSize()<<"\n";
				std::cerr.flush();

				if (verbose) std::cerr<<"Applying exp(iHt)\n";
				HilbertVectorType phi5 = engine.newState(verbose);
				eihOp.apply(phi5,phi4);
				phi4.clear();
				std::cerr<<"time="<<time<<" sigma= "<<sigma<<" sigma2="<<sigma2<<" phi5.terms="<<phi5.terms()<<"\n";
				std::cerr<<"memory="<<myUsage.vmSize()<<"\n";
				std::cerr.flush();

				if (verbose) std::cerr<<"Applying c_p\n";
				FreeOperatorType myOp6 = engine.newSimpleOperator("destruction",site3,sigma3);
				HilbertVectorType phi6 = engine.newState(verbose);
				myOp6.apply(phi6,phi5,simplifyOrNot);
				phi5.clear();
				std::cerr<<"time="<<time<<" sigma= "<<sigma<<" sigma2="<<sigma2<<" phi6.terms="<<phi6.terms()<<"\n";
				std::cerr<<"memory="<<myUsage.vmSize()<<"\n";
				std::cerr.flush();

				if (verbose) std::cerr<<"Adding "<<sigma<<" "<<sigma2<<" "<<it<<"\n";
				sum += scalarProduct(phi6,phi6);
				if (verbose) std::cerr<<"Done with scalar product\n";
				if (sigma ==0 && sigma2 ==0) savedVector = phi6;
				if (sigma ==1 && sigma2 ==1) {
					savedValue = scalarProduct(phi6,savedVector);
					savedVector.clear();
				}
				std::cerr<<"time="<<time<<" sigma= "<<sigma<<" sigma2="<<sigma2<<" saved.terms="<<savedVector.terms()<<"\n";
				std::cerr<<"memory="<<myUsage.vmSize()<<"\n";
				std::cerr.flush();
				phi6.clear();
				//timeVector.add(phi6);
			}
		}
		sum += 2*real(savedValue);
		std::cout<<time<<" "<<real(sum)<<"\n";
	}
}
