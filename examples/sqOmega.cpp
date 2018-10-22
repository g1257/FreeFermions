// SAmple of how to use FreeFermions core engine to calculate
// <A_i 1/(z-H) A_j>
#include <cstdlib>
#define USE_PTHREADS_OR_NOT_NG
#include "Engine.h"
#include "GeometryLibrary.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "OneOverZminusH.h"
#include "DiagonalOperator.h"
#include "GeometryParameters.h"
#include "Concurrency.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "LibraryOperator.h"
#include "Parallelizer.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef ComplexType FieldType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
typedef FreeFermions::GeometryParameters<RealType,InputNgType::Readable>
GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType>
GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::OneOverZminusH<EngineType> OneOverZminusHType;
typedef FreeFermions::DiagonalOperator<OneOverZminusHType> DiagonalOperatorType;
typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
typedef DiagonalOperatorType::FactoryType OpDiagonalFactoryType;
typedef OperatorType::FactoryType OpNormalFactoryType;
typedef FreeFermions::LibraryOperator<OperatorType> LibraryOperatorType;
typedef LibraryOperatorType::FactoryType OpLibFactoryType;
typedef PsimagLite::Matrix<ComplexType> MatrixComplexType;

enum ObservableEnum {OBS_SZ, OBS_C};

void usage(const PsimagLite::String& thisFile)
{
	std::cout<<thisFile<<": USAGE IS "<<thisFile<<" ";
	std::cout<<" -n sites -e electronsUp -g geometry,[leg,filename]\n";
}

struct SqOmegaParams {
	SqOmegaParams(const EngineType& engine_,
	              const HilbertStateType& gs_,
	              RealType Eg_,
	              SizeType sites_,
	              SizeType centralSite_,
	              ObservableEnum what_)
	    : engine(engine_),
	      gs(gs_),
	      Eg(Eg_),
	      sites(sites_),
	      centralSite(centralSite_),
	      observable(what_)
	{}

	const EngineType& engine;
	const HilbertStateType& gs;
	RealType Eg;
	SizeType sites;
	SizeType centralSite;
	ObservableEnum observable;
}; // struct SqOmegaParams

class DynamicObservable {

public:

	DynamicObservable(ObservableEnum what,
	                  const EngineType& engine)
	    : what_(what),opLibFactory_(engine),opNormalFactory_(engine)
	{}

	void applyRightOperator(HilbertStateType& phiKet,
	                        SizeType site0,
	                        SizeType sigma0,
	                        SizeType dynType,
	                        SizeType threadNum)
	{
		if (what_ == OBS_SZ) {
			LibraryOperatorType& nOpI = opLibFactory_(LibraryOperatorType::N,
			                                          site0,
			                                          sigma0,
			                                          threadNum);
			nOpI.applyTo(phiKet);
			return;
		}

		assert(what_ == OBS_C);

		if (dynType == 0) {
			OperatorType& nOpI =  opNormalFactory_(LibraryOperatorType::DESTRUCTION,
			                                       site0,
			                                       sigma0,
			                                       threadNum);
			nOpI.applyTo(phiKet);
			return;
		}

		assert(dynType == 1);
		OperatorType& nOpI =  opNormalFactory_(LibraryOperatorType::CREATION,
		                                       site0,
		                                       sigma0,
		                                       threadNum);
		nOpI.applyTo(phiKet);
	}

	void applyLeftOperator(HilbertStateType& phiKet,
	                       SizeType site0,
	                       SizeType sigma0,
	                       SizeType dynType,
	                       SizeType threadNum)
	{
		if (what_ == OBS_SZ) {
			LibraryOperatorType& nOpI = opLibFactory_(LibraryOperatorType::N,
			                                          site0,
			                                          sigma0,
			                                          threadNum);
			nOpI.applyTo(phiKet);
			return;
		}

		assert(what_ == OBS_C);

		if (dynType == 0) {
			OperatorType& nOpI =  opNormalFactory_(LibraryOperatorType::CREATION,
			                                       site0,
			                                       sigma0,
			                                       threadNum);
			nOpI.applyTo(phiKet);
			return;
		}

		assert(dynType == 1);
		OperatorType& nOpI =  opNormalFactory_(LibraryOperatorType::DESTRUCTION,
		                                       site0,
		                                       sigma0,
		                                       threadNum);
		nOpI.applyTo(phiKet);
	}

	RealType signForWeight(SizeType sigma0,
	                       SizeType sigma1,
	                       SizeType dynType)
	{
		RealType sign = (sigma0 == sigma1) ? 1.0 : -1.0;
		sign *= (dynType == 0) ? 1 : -1;
		return (what_ == OBS_SZ) ? sign : -1.0;
	}

private:

	ObservableEnum what_;
	OpLibFactoryType opLibFactory_;
	OpNormalFactoryType opNormalFactory_;
}; // class DynamicObservable

class SqOmegaParallel {

public:

	SqOmegaParallel(const SqOmegaParams& params,
	                RealType total,
	                RealType step,
	                RealType offset)
	    : params_(params),
	      step_(step),
	      offset_(offset),
	      observable_(params.observable,params.engine),
	      result_(total,params.sites)
	{}

	void doTask(SizeType taskNumber,SizeType threadNum)
	{
		RealType omega = taskNumber*step_ + offset_;
		doOneOmega(taskNumber,omega, threadNum);
	}

	SizeType tasks() const { return result_.n_row(); }

	void print(std::ostream& os) const
	{
		for (SizeType it = 0; it< result_.n_row(); it++) {
			RealType omega = it * step_ + offset_;
			os<<omega<<" ";
			for (SizeType site1 = 0; site1 < result_.n_col(); ++site1) {
				ComplexType val = result_(it,site1);
				os<<PsimagLite::real(val)<<" "<<PsimagLite::imag(val)<<" ";
			}

			os<<"\n";
		}
	}

private:

	FieldType doOneOmegaOneSitePair(SizeType site0,
	                                SizeType site1,
	                                RealType omega,
	                                SizeType threadNum)
	{
		FieldType tmpC = 0.0;
		RealType epsilon = 0.1;
		SizeType sigma0 = 0;
		SizeType sigma1 = 0;
		for (SizeType dynType = 0; dynType < 2; ++dynType) {
			RealType sign = observable_.signForWeight(sigma0,sigma1,dynType);
			RealType signForDen = (dynType == 1) ? -1.0 : 1.0;
			HilbertStateType phiKet = params_.gs;
			observable_.applyRightOperator(phiKet,
			                               site0,
			                               sigma0,
			                               dynType,
			                               threadNum);

			HilbertStateType phiBra = params_.gs;
			observable_.applyLeftOperator(phiBra,
			                              site1,
			                              sigma1,
			                              dynType,
			                              threadNum);

			OpDiagonalFactoryType opDiagonalFactory(params_.engine);

			ComplexType z = ComplexType(omega,epsilon);
			OneOverZminusHType eih(z,signForDen,params_.Eg,params_.engine);
			DiagonalOperatorType& eihOp = opDiagonalFactory(eih);
			HilbertStateType phi3 = phiKet;
			eihOp.applyTo(phi3);

			tmpC += sign*scalarProduct(phiBra,phi3);
		}

		return tmpC;
	}

	void doOneOmega(SizeType it, RealType omega, SizeType threadNum)
	{
		SizeType site0 = params_.centralSite;
		for (SizeType site1 = 0; site1 < params_.sites; ++site1) {
			std::cerr<<"site1="<<site1<<"\n";
			result_(it,site1) = doOneOmegaOneSitePair(site0,
			                                          site1,
			                                          omega,
			                                          threadNum);
		}
	}

	const SqOmegaParams& params_;
	RealType step_;
	RealType offset_;
	DynamicObservable observable_;
	MatrixComplexType result_;
};

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	SizeType total=0;
	RealType offset = 0;
	RealType step = 0;
	SizeType centralSite =0;
	ObservableEnum what = OBS_SZ;

	while ((opt = getopt(argc, argv, "f:t:o:i:c:w:")) != -1) {
		switch (opt) {
		case 'f':
			file=optarg;
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
		case 'c':
			centralSite = atoi(optarg);
			break;
		case 'w':
			if (PsimagLite::String(optarg) == "c")
				what = OBS_C;
			else if (PsimagLite::String(optarg) == "sz")
				what = OBS_SZ;
			else
				throw std::runtime_error("observable" +
			                             PsimagLite::String(optarg) +
			                             " not supported\n");
			break;
		default: /* '?' */
			throw std::runtime_error("Wrong usage\n");
		}
	}

	if (file=="") {
		throw std::runtime_error("Wrong usage\n");
	}

	FreeFermions::InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(file,inputCheck);
	InputNgType::Readable io(ioWriteable);

	GeometryParamsType geometryParams(io);
	SizeType electronsUp = GeometryParamsType::readElectrons(io,geometryParams.sites);

	SizeType dof = 1;

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;

	SizeType npthreads = 1;
	io.readline(npthreads,"Threads=");
	ConcurrencyType concurrency(&argc,&argv,npthreads);
	EngineType engine(geometry.matrix(),
	                  geometryParams.outputFile,
	                  dof,
	                  EngineType::VERBOSE_YES);
	PsimagLite::Vector<SizeType>::Type ne(dof,electronsUp);
	bool debug = false;
	HilbertStateType gs(engine,ne,debug);

	RealType Eg = 0;
	for (SizeType i=0;i<ne[0];i++) Eg += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*Eg<<"\n";

	std::cout<<"#TotalNumberOfSites="<<geometryParams.sites<<"\n";
	std::cout<<"#OmegaTotal="<<total<<"\n";
	std::cout<<"#OmegaBegin="<<offset<<"\n";
	std::cout<<"#OmegaStep="<<step<<"\n";
	std::cout<<"#GeometryKind="<<geometryParams.geometry<<"\n";
	std::cout<<"#TSPSites 1 "<<centralSite<<"\n";
	std::cout<<"#Threads="<<PsimagLite::Concurrency::codeSectionParams.npthreads<<"\n";
	std::cout<<"#What="<<what<<"\n";
	std::cout<<"#############\n";

	typedef PsimagLite::Parallelizer<SqOmegaParallel> ParallelizerType;
	ParallelizerType threadObject(PsimagLite::Concurrency::codeSectionParams);

	SqOmegaParams params(engine,gs,Eg,geometryParams.sites,centralSite,what);
	SqOmegaParallel helperSqOmega(params,total,step,offset);

	threadObject.loopCreate(helperSqOmega);

	helperSqOmega.print(std::cout);
}
