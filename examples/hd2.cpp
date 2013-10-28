

// Calculates <phi | phi>
// where 
// |phi> = c_{p up} exp(iHt) c_{i\sigma} nbar_{i \bar{sigma}} c^dagger_{j sigma'} n_{j \bar{sigma'}} |gs>

typedef double RealType;
#include <cstdlib>
#include "unistd.h"
#include "Engine.h"
#include "GeometryLibrary.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"
#include "LibraryOperator.h"
#include "Tokenizer.h" // in PsimagLite
#include "GeometryParameters.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "InputNg.h"
#include "InputCheck.h"

typedef std::complex<double> ComplexType;
typedef ComplexType FieldType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
typedef FreeFermions::GeometryParameters<RealType,InputNgType::Readable> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;
typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
typedef FreeFermions::LibraryOperator<OperatorType> LibraryOperatorType;
typedef OperatorType::FactoryType OpNormalFactoryType;
typedef LibraryOperatorType::FactoryType OpLibFactoryType;
typedef DiagonalOperatorType::FactoryType OpDiagonalFactoryType;

FieldType calcSuperDensity(size_t site,
						   size_t site2,
						   const HilbertStateType& gs,
						   const EngineType& engine)
{
	HilbertStateType savedVector = gs;
	FieldType savedValue = 0;
	FieldType sum = 0;
	OpNormalFactoryType opNormalFactory(engine);
	OpLibFactoryType opLibFactory(engine);
	
	for (size_t sigma = 0;sigma<2;sigma++) {
		HilbertStateType phi = gs;
		
		LibraryOperatorType& myOp = opLibFactory(
			LibraryOperatorType::N,site,1-sigma);
		
		myOp.applyTo(phi);
		OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,
											  site,sigma);
		
		myOp2.applyTo(phi);
		
		for (size_t sigma2 = 0;sigma2 < 2;sigma2++) {
			HilbertStateType phi3 = phi;
			LibraryOperatorType& myOp3 = opLibFactory(
				LibraryOperatorType::NBAR,site2,1-sigma2);
			myOp3.applyTo(phi3);
			
			OperatorType& myOp4 = opNormalFactory(
				OperatorType::DESTRUCTION,site2,sigma2);
			
			myOp4.applyTo(phi3);
			
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

class MyLoop {

	typedef PsimagLite::Concurrency ConcurrencyType;

	enum {SPIN_UP,SPIN_DOWN};

public:

	MyLoop(EngineType& engine,
	       SizeType step,
	       SizeType offset,
	       const HilbertStateType& gs,
	       PsimagLite::Vector<SizeType>::Type& sites,
	       bool verbose)
	    : engine_(engine),
	      step_(step),
	      offset_(offset),
	      gs_(gs),
	      sites_(sites),
	      verbose_(verbose)
	{}

	void thread_function_(SizeType threadNum,
	                      SizeType blockSize,
	                      SizeType total,
	                      ConcurrencyType::MutexType* myMutex)
	{
		SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
		SizeType npthreads = ConcurrencyType::npthreads;
		for (SizeType p=0;p<blockSize;p++) {
			SizeType it = (threadNum+npthreads*mpiRank)*blockSize + p;
			if (it>=total) break;

			OpNormalFactoryType opNormalFactory(engine_);
			OpLibFactoryType opLibFactory(engine_);
			OpDiagonalFactoryType opDiagonalFactory(engine_);

			RealType time = it * step_ + offset_;
			EtoTheIhTimeType eih(time,engine_,0);
			DiagonalOperatorType& eihOp = opDiagonalFactory(eih);

			PsimagLite::Vector<HilbertStateType*>::Type savedVector(4);

			for (size_t sigma = 0;sigma<2;sigma++) {
				HilbertStateType phi = gs_;
				LibraryOperatorType& myOp = opLibFactory(
				            LibraryOperatorType::N,sites_[0],1-sigma);
				myOp.applyTo(phi);

				OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,
				                                      sites_[0],sigma);
				myOp2.applyTo(phi);

				for (size_t sigma2 = 0;sigma2 < 2;sigma2++) {
					savedVector[sigma+sigma2*2] = new HilbertStateType(phi);


					LibraryOperatorType& myOp3 = opLibFactory(
					            LibraryOperatorType::NBAR,sites_[1],1-sigma2);
					myOp3.applyTo(*savedVector[sigma+sigma2*2]);

					OperatorType& myOp4 = opNormalFactory(
					            OperatorType::DESTRUCTION,sites_[1],sigma2);
					myOp4.applyTo(*savedVector[sigma+sigma2*2]);

					if (verbose_) std::cerr<<"Applying exp(iHt)\n";
					eihOp.applyTo(*savedVector[sigma+sigma2*2]);

					if (verbose_) std::cerr<<"Applying c_{p down}\n";
					OperatorType& myOp6 = opNormalFactory(
					            OperatorType::DESTRUCTION,sites_[2],SPIN_DOWN);
					myOp6.applyTo(*savedVector[sigma+sigma2*2]);

					if (verbose_) std::cerr<<"Applying c_{p up}\n";
					OperatorType& myOp7 = opNormalFactory(
					            OperatorType::DESTRUCTION,sites_[2],SPIN_UP);
					myOp7.applyTo(*savedVector[sigma+sigma2*2]);
				}
			}
			FieldType sum = 0;
			size_t total = savedVector.size()*savedVector.size()/2;
			for (size_t x=0;x<total;x++) {
				size_t sigma = (x & 3);
				size_t sigma2 = (x & 12);
				sigma2 >>= 2;
				sum += scalarProduct(*savedVector[sigma],*savedVector[sigma2]);
			}
			for (size_t x=0;x<savedVector.size();x++) delete savedVector[x];

			std::cout<<time<<" "<<(2.0*sum)<<"\n";
		}
	}

private:

	EngineType& engine_;
	SizeType step_;
	SizeType offset_;
    const HilbertStateType& gs_;
    PsimagLite::Vector<SizeType>::Type& sites_;
    bool verbose_;
};

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	size_t total=0;
	RealType offset = 0;
	RealType step = 0;
	size_t site3=0;

	while ((opt = getopt(argc, argv, "f:p:t:o:i:")) != -1) {
		switch (opt) {
		case 'f':
			file=optarg;
			break;
		case 'p':
			site3 = atoi(optarg);
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
		default: /* '?' */
			throw std::runtime_error("Wrong usage\n");
		}
	}
	if (file=="") throw std::runtime_error("Wrong usage\n");
	
	FreeFermions::InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(file,inputCheck);
	InputNgType::Readable io(ioWriteable);

	GeometryParamsType geometryParams(io);
	size_t electronsUp = GeometryParamsType::readElectrons(io,geometryParams.sites);

	PsimagLite::Vector<size_t>::Type sites;
	GeometryParamsType::readVector(sites,file,"TSPSites");
	SizeType nthreads = 1;
	GeometryParamsType::readLabel(nthreads,file,"Threads=");
	sites.resize(3);
	sites[2]=site3;

	size_t dof = 2; // spin up and down
	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;

	typedef PsimagLite::Concurrency ConcurrencyType;
	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc,&argv,npthreads);

	EngineType engine(geometry.matrix(),dof,false);
	
	PsimagLite::Vector<size_t>::Type ne(dof,electronsUp); // 8 up and 8 down
	bool debug = false;
	bool verbose = false;
	HilbertStateType gs(engine,ne,debug);

	RealType sum = 0;
	for (size_t i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";

	FieldType superdensity = calcSuperDensity(sites[0],sites[1],gs,engine);
	std::cout<<"#superdensity="<<superdensity<<"\n";

	std::cout<<"#sites= ";
	for (size_t i=0;i<sites.size();i++) std::cout<<sites[i]<<" ";
	std::cout<<"\n";
	
	typedef MyLoop MyLoopType;
	typedef PsimagLite::Parallelizer<MyLoopType> ParallelizerType;
	ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
	                              PsimagLite::MPI::COMM_WORLD);

	MyLoopType myLoop(engine,step,offset,gs,sites,verbose);

	std::cout<<"Using "<<threadObject.name();
	std::cout<<" with "<<threadObject.threads()<<" threads.\n";
	threadObject.loopCreate(total,myLoop);

}
