// SAmple of how to use FreeFermions core engine to calculate
// <A_i 1/(z-H) A_j>
#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "OneOverZminusH.h"
#include "DiagonalOperator.h"
#include "Tokenizer.h"
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

enum {DYN_TYPE_0,DYN_TYPE_1};

void usage(const PsimagLite::String& thisFile)
{
	std::cout<<thisFile<<": USAGE IS "<<thisFile<<" ";
	std::cout<<" -n sites -e electronsUp -g geometry,[leg,filename]\n";
}

void setMyGeometry(GeometryParamsType& geometryParams,
                   const PsimagLite::Vector<PsimagLite::String>::Type& vstr)
{
	// default value
	geometryParams.type = GeometryLibraryType::CHAIN;

	if (vstr.size()<2) {
		// assume chain
		return;
	}

	PsimagLite::String gName = vstr[0];
	if (gName == "chain") {
		PsimagLite::String str("setMyGeometry: ");
		throw std::runtime_error(str + "-g chain takes no further arguments\n");
	}

	geometryParams.leg = atoi(vstr[1].c_str());

	if (gName == "ladder") {
		if (vstr.size()!=3) {
			usage("setMyGeometry");
			PsimagLite::String str("setMyGeometry: usage is: ");
			throw std::runtime_error(str + "-g ladder,leg,isPeriodic\n");
		}
		geometryParams.type = GeometryLibraryType::LADDER;
		geometryParams.hopping.resize(2);
		geometryParams.hopping[0] =  geometryParams.hopping[1]  = 1.0;
		geometryParams.isPeriodic[GeometryParamsType::DIRECTION_Y] =
		        (atoi(vstr[2].c_str())>0);
		return;
	}

	if (vstr.size()!=3) {
		usage("setMyGeometry");
		PsimagLite::String str("setMyGeometry: usage is: ");
		throw std::runtime_error(str + "-g {feas | ktwoniffour} leg filename\n");
	}

	geometryParams.filename = vstr[2];

	if (gName == "feas") {
		geometryParams.type = GeometryLibraryType::FEAS;
		return;
	}

	if (gName == "kniffour") {
		geometryParams.type = GeometryLibraryType::KTWONIFFOUR;
		return;
	}
}

struct SqOmegaParams {
	SqOmegaParams(const EngineType& engine_,
	              const HilbertStateType& gs_,
	              RealType Eg_,
	              SizeType sites_,
	              SizeType centralSite_)
	    : engine(engine_),
	      gs(gs_),
	      Eg(Eg_),
	      sites(sites_),
	      centralSite(centralSite_)
	{}

	const EngineType& engine;
	const HilbertStateType& gs;
	RealType Eg;
	SizeType sites;
	SizeType centralSite;
}; // struct SqOmegaParams

class SqOmegaParallel {

public:

	SqOmegaParallel(const SqOmegaParams& params,
	                RealType total,
	                RealType step,
	                RealType offset)
	    : params_(params),step_(step), offset_(offset),result_(total,params.sites)
	{}

	void thread_function_(SizeType threadNum,
	                      SizeType blockSize,
	                      SizeType total,
	                      typename ConcurrencyType::MutexType*)
	{
		SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
		SizeType npthreads = ConcurrencyType::npthreads;
		for (SizeType p=0;p<blockSize;p++) {
			SizeType it = (threadNum+npthreads*mpiRank)*blockSize + p;
			if (it>=total) continue;
			RealType omega = it * step_ + offset_;
			doOneOmega(it,omega);
		}
	}

	void print(std::ostream& os) const
	{
		for (SizeType it = 0; it< result_.n_row(); it++) {
			RealType omega = it * step_ + offset_;
			os<<omega<<" ";
			for (SizeType site1 = 0; site1 < result_.n_col(); ++site1) {
				ComplexType val = result_(it,site1);
				os<<std::real(val)<<" "<<std::imag(val)<<" ";
			}

			os<<"\n";
		}
	}

private:

	FieldType doOneOmegaOneSitePair(OpLibFactoryType& opLibFactory,
	                                SizeType site0,
	                                SizeType site1,
	                                RealType omega)
	{
		FieldType tmpC = 0.0;
		RealType epsilon = 0.1;
		SizeType sigma0 = 0;
		SizeType sigma1 = 0;
		for (SizeType dynType = 0; dynType < 2; ++dynType) {
			RealType sign = (sigma0 == sigma1) ? 1.0 : -1.0;
			sign *= (dynType == 0) ? 1 : -1;
			RealType signForDen = (dynType== DYN_TYPE_1) ? -1.0 : 1.0;
			LibraryOperatorType& nOpI = opLibFactory(LibraryOperatorType::N,
			                                         site0,
			                                         sigma0);

			HilbertStateType phiKet = params_.gs;
			nOpI.applyTo(phiKet);

			LibraryOperatorType& nOpJ = opLibFactory(LibraryOperatorType::N,
			                                         site1,
			                                         sigma1);
			HilbertStateType phiBra = params_.gs;
			nOpJ.applyTo(phiBra);

			//FieldType density = scalarProduct(phiBra,phiKet);
			//std::cerr<<"density="<<density<<"\n";

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

	void doOneOmega(SizeType it, RealType omega)
	{
		OpLibFactoryType opLibFactory(params_.engine);
		SizeType site0 = params_.centralSite;
		for (SizeType site1 = 0; site1 < params_.sites; ++site1) {
			result_(it,site1) = doOneOmegaOneSitePair(opLibFactory,
			                                          site0,
			                                          site1,
			                                          omega);
		}
	}

	const SqOmegaParams& params_;
	RealType step_;
	RealType offset_;
	MatrixComplexType result_;
};

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	SizeType total=0;
	RealType offset = 0;
	RealType step = 0;

	while ((opt = getopt(argc, argv, "f:t:o:i:")) != -1) {
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
	EngineType engine(geometry.matrix(),dof,true);

	PsimagLite::Vector<SizeType>::Type ne(dof,electronsUp);
	bool debug = false;
	HilbertStateType gs(engine,ne,debug);

	RealType Eg = 0;
	for (SizeType i=0;i<ne[0];i++) Eg += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*Eg<<"\n";

	SizeType centralSite = static_cast<SizeType>(geometryParams.sites/2)-1;
	std::cout<<"#TotalNumberOfSites="<<geometryParams.sites<<"\n";
	std::cout<<"#OmegaTotal="<<total<<"\n";
	std::cout<<"#OmegaBegin="<<offset<<"\n";
	std::cout<<"#OmegaStep="<<step<<"\n";
	std::cout<<"#GeometryKind="<<geometryParams.geometry<<"\n";
	std::cout<<"#TSPSites 1 "<<centralSite<<"\n";
	std::cout<<"#Threads="<<PsimagLite::Concurrency::npthreads<<"\n";
	std::cout<<"#############\n";

	typedef PsimagLite::Parallelizer<SqOmegaParallel> ParallelizerType;
	ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
	                              PsimagLite::MPI::COMM_WORLD);

	SqOmegaParams params(engine,gs,Eg,geometryParams.sites,centralSite);
	SqOmegaParallel helperSqOmega(params,total,step,offset);

	threadObject.loopCreate(total,helperSqOmega);

	helperSqOmega.print(std::cout);
}

