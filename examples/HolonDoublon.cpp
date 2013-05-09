

// Calculates <phi | phi>
// where 
// |phi> = c_{p up} exp(iHt) c_{i\sigma} nbar_{i \bar{sigma}} c^dagger_{j sigma'} n_{j \bar{sigma'}} |gs>

#include <cstdlib>
#include "unistd.h"
#include "Engine.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
#include "TypeToString.h"
#include "Tokenizer.h" // in PsimagLite
#include "GeometryParameters.h"
#include "ParallelHolonDoublon.h"
#ifdef USE_PTHREADS
#include "Pthreads.h"
#define PTHREADS_NAME PsimagLite::Pthreads
#else
#include "NoPthreads.h"
#define PTHREADS_NAME PsimagLite::NoPthreads
#endif

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef ComplexType FieldType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::GeometryParameters<RealType> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType,ConcurrencyType> EngineType;

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	size_t total=0;
	RealType offset = 0;
	size_t site3 = 0;

	RealType step = 0;
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

	GeometryParamsType geometryParams(file);
	size_t electronsUp = GeometryParamsType::readElectrons(file,geometryParams.sites);
	PsimagLite::Vector<size_t>::Type sites;
	GeometryParamsType::readVector(sites,file,"TSPSites");
	sites.resize(3);
	sites[2] = site3;
	size_t nthreads = 0;
	GeometryParamsType::readLabel(nthreads,file,"Threads=");
	assert(nthreads>0);

	size_t dof = 2; // spin up and down

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;
	
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(geometry,concurrency,dof,false);

	PsimagLite::Vector<size_t>::Type ne(dof,electronsUp); // 8 up and 8 down
	bool debug = false;
	bool verbose = false;

	
	RealType sum = 0;
	for (size_t i=0;i<electronsUp;i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";

	size_t sigma3 = 0;

	std::cout<<"#sites= ";
	for (size_t i=0;i<sites.size();i++) std::cout<<sites[i]<<" ";
	std::cout<<"\n";

	typedef FreeFermions::ParallelHolonDoublon<RealType,FieldType,EngineType> ParallelHolonDoublonType;
	PTHREADS_NAME<ParallelHolonDoublonType> threadedHolonDoublon;
	PTHREADS_NAME<ParallelHolonDoublonType>::setThreads(nthreads);

	ParallelHolonDoublonType::HolonDoublonParamsType params(ne,
															sites,
															sigma3,
															offset,
															step,
															debug,
															verbose);
	ParallelHolonDoublonType helperHolonDoublon(engine,params);

	FieldType superdensity = helperHolonDoublon.calcSuperDensity(sites[0],sites[1]);
	std::cout<<"#superdensity="<<superdensity<<"\n";

	threadedHolonDoublon.loopCreate(total,helperHolonDoublon,concurrency);
}
