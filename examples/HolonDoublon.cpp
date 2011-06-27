

// Calculates <phi | phi>
// where 
// |phi> = c_{p up} exp(iHt) c_{i\sigma} nbar_{i \bar{sigma}} c^dagger_{j sigma'} n_{j \bar{sigma'}} |gs>

#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"
#include "LibraryOperator.h"

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
typedef FreeFermions::LibraryOperator<OperatorType> LibraryOperatorType;
typedef typename OperatorType::FactoryType OpNormalFactoryType;
typedef typename LibraryOperatorType::FactoryType OpLibFactoryType;

FieldType calcSuperDensity(size_t site,
                             size_t site2,
                             const HilbertStateType& gs,
                             const EngineType& engine)
{
	HilbertStateType savedVector = gs;
	FieldType savedValue = 0;
	FieldType sum = 0;
	OpNormalFactoryType opNormalFactory;
	OpLibFactoryType opLibFactory;

	for (size_t sigma = 0;sigma<2;sigma++) {
		HilbertStateType phi = gs;

		LibraryOperatorType* myOp = opLibFactory(engine,
		                         LibraryOperatorType::N,site,1-sigma);

		myOp->applyTo(phi);
		OperatorType* myOp2 = opNormalFactory(engine,OperatorType::CREATION,
		                         site,sigma);

		myOp2->applyTo(phi);
		
		for (size_t sigma2 = 0;sigma2 < 2;sigma2++) {
			HilbertStateType phi3 = phi;
			LibraryOperatorType* myOp3 = opLibFactory(engine,
			                         LibraryOperatorType::NBAR,site2,1-sigma2);
			myOp3->applyTo(phi3);

			OperatorType* myOp4 = opNormalFactory(engine,
			                         OperatorType::DESTRUCTION,site2,sigma2);

			myOp4->applyTo(phi3);

			if (sigma ==0 && sigma2 ==0) savedVector = phi3;
			sum += scalarProduct(phi3,phi3);

			if (sigma ==1 && sigma2 ==1) {
				savedValue = scalarProduct(phi3,savedVector);
			}
		}
	}
	sum += 2*real(savedValue);
	//std::cerr<<"#sum2="<<scalarProduct(tmpV,tmpV)<<"\n";
	return sum;
}


int main(int argc,char *argv[])
{
	if (argc!=7) throw std::runtime_error("Needs 7 arguments\n");
	size_t n = 32;
	size_t electronsUp = 16;
	size_t dof = 2; // spin up and down
	
	MatrixType t(n,n);
	GeometryLibraryType geometry(n,GeometryLibraryType::CHAIN);
	geometry.setGeometry(t);

	/* GeometryLibraryType geometry(n,GeometryLibraryType::LADDER);
	geometry.setGeometry(t,2);

	for (size_t ii=0;ii<n;ii+=2)
		t(ii,ii+1) = t(ii+1,ii) = 0.5;
	*/
	std::cerr<<t;
	
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(t,concurrency,dof,false);
	
	std::vector<size_t> ne(dof,electronsUp); // 8 up and 8 down
	bool debug = false;
	bool verbose = false;
	HilbertStateType gs(engine.size(),ne,debug);
	
	size_t site = atoi(argv[1]);
	size_t site2 = atoi(argv[2]);
	size_t site3 = atoi(argv[3]);
	size_t sigma3 = 0;
	
	FieldType superdensity = calcSuperDensity(site,site2,gs,engine);
	std::cout<<"#superdensity="<<superdensity<<"\n";
	std::cout<<"#site="<<site<<" site2="<<site2<<"\n";

	concurrency.loopCreate(size_t(atoi(argv[4])));
	size_t it = 0;

	while(concurrency.loop(it)) {
		OpNormalFactoryType opNormalFactory;
		OpLibFactoryType opLibFactory;

		RealType time = it * atof(argv[5]) + atof(argv[6]);
		EtoTheIhTimeType eih(time,engine,0);
		DiagonalOperatorType eihOp(eih);

		HilbertStateType savedVector = gs;
		FieldType savedValue = 0;
		FieldType sum = 0;

		for (size_t sigma = 0;sigma<2;sigma++) {
			HilbertStateType phi = gs;
			LibraryOperatorType* myOp = opLibFactory(engine,
					                      LibraryOperatorType::N,site,1-sigma);
			myOp->applyTo(phi);

			OperatorType* myOp2 = opNormalFactory(engine,OperatorType::CREATION,
					                         site,sigma);
			myOp2->applyTo(phi);

			for (size_t sigma2 = 0;sigma2 < 2;sigma2++) {
				HilbertStateType phi3 = phi;


				LibraryOperatorType* myOp3 = opLibFactory(engine,
				                 LibraryOperatorType::NBAR,site2,1-sigma2);
				myOp3->applyTo(phi3);

				OperatorType* myOp4 = opNormalFactory(engine,
				                 OperatorType::DESTRUCTION,site2,sigma2);
				myOp4->applyTo(phi3);

				if (verbose) std::cerr<<"Applying exp(iHt)\n";
				eihOp.applyTo(phi3);

				if (verbose) std::cerr<<"Applying c_p\n";
				OperatorType* myOp6 = opNormalFactory(engine,
								  OperatorType::DESTRUCTION,site3,sigma3);
				myOp6->applyTo(phi3);

				if (verbose) std::cerr<<"Adding "<<sigma<<" "<<sigma2<<" "<<it<<"\n";

				if (sigma ==0 && sigma2 ==0) savedVector = phi3;
				if (sigma ==1 && sigma2 ==1) {
					savedValue = scalarProduct(phi3,savedVector);
				}
				sum += scalarProduct(phi3,phi3);
				if (verbose) std::cerr<<"Done with scalar product\n";
			}
		}
		sum += 2*real(savedValue);
		std::cout<<time<<" "<<real(sum)<<"\n";
	}
}
