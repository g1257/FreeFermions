

// Calculates <phi | phi>
// where
// |phi> = c_{p up} exp(iHt) c_{i\sigma} nbar_{i \bar{sigma}} c^dagger_{j sigma'} n_{j \bar{sigma'}} |gs>

typedef double RealType;
#include <cstdlib>
#include <unistd.h>
#define USE_PTHREADS_OR_NOT_NG
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

FieldType calcSuperDensity(SizeType site,
                           SizeType site2,
                           const HilbertStateType& gs,
                           const EngineType& engine)
{
	HilbertStateType savedVector = gs;
	FieldType savedValue = 0;
	FieldType sum = 0;
	OpNormalFactoryType opNormalFactory(engine);
	OpLibFactoryType opLibFactory(engine);

	for (SizeType sigma = 0;sigma<2;sigma++) {
		HilbertStateType phi = gs;

		LibraryOperatorType& myOp = opLibFactory(
		            LibraryOperatorType::N,site,1-sigma);

		myOp.applyTo(phi);
		OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,
		                                      site,sigma);

		myOp2.applyTo(phi);

		for (SizeType sigma2 = 0;sigma2 < 2;sigma2++) {
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
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorComplexType;

	enum {SPIN_UP,SPIN_DOWN};

public:

	MyLoop(EngineType& engine,
	       RealType step,
	       RealType offset,
	       const HilbertStateType& gs,
	       PsimagLite::Vector<SizeType>::Type& sites,
	       SizeType total,
	       bool verbose)
	    : engine_(engine),
	      step_(step),
	      offset_(offset),
	      gs_(gs),
	      sites_(sites),
	      data_(total),
	      verbose_(verbose)
	{}

	SizeType tasks() const { return data_.size(); }

	void doTask(SizeType taskNumber, SizeType)
	{
		OpNormalFactoryType opNormalFactory(engine_);
		OpLibFactoryType opLibFactory(engine_);
		OpDiagonalFactoryType opDiagonalFactory(engine_);

		RealType time = taskNumber * step_ + offset_;
		EtoTheIhTimeType eih(time,engine_,0);
		DiagonalOperatorType& eihOp = opDiagonalFactory(eih);

		PsimagLite::Vector<HilbertStateType*>::Type savedVector(4);

		for (SizeType sigma = 0;sigma<2;sigma++) {
			HilbertStateType phi = gs_;
			LibraryOperatorType& myOp = opLibFactory(
			            LibraryOperatorType::N,sites_[0],1-sigma);
			myOp.applyTo(phi);

			OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,
			                                      sites_[0],sigma);
			myOp2.applyTo(phi);

			for (SizeType sigma2 = 0;sigma2 < 2;sigma2++) {
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
		SizeType total = savedVector.size()*savedVector.size()/2;
		for (SizeType x=0;x<total;x++) {
			SizeType sigma = (x & 3);
			SizeType sigma2 = (x & 12);
			sigma2 >>= 2;
			sum += scalarProduct(*savedVector[sigma],*savedVector[sigma2]);
		}

		for (SizeType x=0;x<savedVector.size();x++) delete savedVector[x];

		data_[taskNumber] = 2.0*sum;
	}

	void printTasks(std::ostream& os) const
	{
		SizeType total = data_.size();
		for (SizeType it = 0; it < total; ++it) {
			RealType time = it * step_ + offset_;
			os<<time<<" "<<data_[it]<<"\n";
		}
	}

private:

	EngineType& engine_;
	RealType step_;
	RealType offset_;
	const HilbertStateType& gs_;
	PsimagLite::Vector<SizeType>::Type& sites_;
	VectorComplexType data_;
	bool verbose_;
};

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	SizeType total=0;
	RealType offset = 0;
	RealType step = 0;
	SizeType site3=0;

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
	SizeType electronsUp = GeometryParamsType::readElectrons(io,geometryParams.sites);

	PsimagLite::Vector<SizeType>::Type sites;
	GeometryParamsType::readVector(sites,file,"TSPSites");
	SizeType nthreads = 1;
	GeometryParamsType::readLabel(nthreads,file,"Threads=");
	sites.resize(3);
	sites[2]=site3;

	SizeType dof = 2; // spin up and down
	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;

	typedef PsimagLite::Concurrency ConcurrencyType;
	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc,&argv,npthreads);

	EngineType engine(geometry.matrix(),dof,false);

	PsimagLite::Vector<SizeType>::Type ne(dof,electronsUp); // 8 up and 8 down
	bool debug = false;
	bool verbose = false;
	HilbertStateType gs(engine,ne,debug);

	RealType sum = 0;
	for (SizeType i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";

	FieldType superdensity = calcSuperDensity(sites[0],sites[1],gs,engine);
	std::cout<<"#superdensity="<<superdensity<<"\n";

	std::cout<<"#sites= ";
	for (SizeType i=0;i<sites.size();i++) std::cout<<sites[i]<<" ";
	std::cout<<"\n";

	typedef MyLoop MyLoopType;
	typedef PsimagLite::Parallelizer<MyLoopType> ParallelizerType;
	ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
	                              PsimagLite::MPI::COMM_WORLD);

	MyLoopType myLoop(engine,step,offset,gs,sites,total,verbose);

	std::cout<<"Using "<<threadObject.name();
	std::cout<<" with "<<threadObject.threads()<<" threads.\n";
	threadObject.loopCreate(myLoop);
	myLoop.printTasks(std::cout);
}
