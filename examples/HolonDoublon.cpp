

// Calculates <phi | phi>
// where 
// |phi> = c_{p up} exp(iHt) c_{i\sigma} nbar_{i \bar{sigma}}
//              c^dagger_{j sigma'} n_{j \bar{sigma'}} |gs>

#include <cstdlib>
#include <unistd.h>
#define USE_PTHREADS_OR_NOT_NG
#include "Engine.h"
#include "GeometryLibrary.h"
#include "TypeToString.h"
#include "GeometryParameters.h"
#include "ParallelHolonDoublon.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "InputNg.h"
#include "InputCheck.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef ComplexType FieldType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
typedef FreeFermions::GeometryParameters<RealType,InputNgType::Readable> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType> EngineType;

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	SizeType total=0;
	RealType offset = 0;
	SizeType site3 = 0;

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

	FreeFermions::InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(file,inputCheck);
	InputNgType::Readable io(ioWriteable);

	GeometryParamsType geometryParams(io);
	SizeType electronsUp = GeometryParamsType::readElectrons(io,geometryParams.sites);
	PsimagLite::Vector<SizeType>::Type sites;
	GeometryParamsType::readVector(sites,file,"TSPSites");
	sites.resize(3);
	sites[2] = site3;
	SizeType nthreads = 1;

	try {
		GeometryParamsType::readLabel(nthreads,file,"Threads=");
	} catch (std::exception& e) {}

	assert(nthreads>0);

	typedef PsimagLite::Concurrency ConcurrencyType;
	ConcurrencyType concurrency(&argc,&argv,nthreads);

	SizeType dof = 2; // spin up and down

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;
	
	EngineType engine(geometry.matrix(),
	                  geometryParams.outputFile,
	                  dof,
	                  EngineType::VerboseEnum::YES);

	PsimagLite::Vector<SizeType>::Type ne(dof,electronsUp); // 8 up and 8 down
	bool debug = false;
	bool verbose = false;

	
	RealType sum = 0;
	for (SizeType i=0;i<electronsUp;i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";

	SizeType sigma3 = 0;

	std::cout<<"#sites= ";
	for (SizeType i=0;i<sites.size();i++) std::cout<<sites[i]<<" ";
	std::cout<<"\n";

	typedef FreeFermions::ParallelHolonDoublon<RealType,
	        FieldType,
	        EngineType> ParallelHolonDoublonType;
	typedef PsimagLite::Parallelizer<ParallelHolonDoublonType> ParallelizerType;
	ParallelizerType threadedHolonDoublon(PsimagLite::Concurrency::codeSectionParams);

	ParallelHolonDoublonType::HolonDoublonParamsType params(ne,
	                                                        sites,
	                                                        sigma3,
	                                                        offset,
	                                                        step,
	                                                        debug,
	                                                        verbose);
	ParallelHolonDoublonType helperHolonDoublon(engine,params,total);

	FieldType superdensity = helperHolonDoublon.calcSuperDensity(sites[0],sites[1]);
	std::cout<<"#superdensity="<<superdensity<<"\n";

	threadedHolonDoublon.loopCreate(helperHolonDoublon);
	helperHolonDoublon.printTasks(std::cout);
}
