// Calculates <phi | phi>
// where 
// |phi> = decay |gs>

#include <cstdlib>
#include <unistd.h>
#define USE_PTHREADS_OR_NOT_NG
#include "Engine.h"
#include "GeometryLibrary.h"
#include "TypeToString.h"
#include "GeometryParameters.h"
#include "ParallelDecay.h"
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
	sites.resize(sites.size()+1);
	sites[sites.size()-1] = site3;
	SizeType nthreads = 1;

	try {
		GeometryParamsType::readLabel(nthreads,file,"Threads=");
	} catch (std::exception& e) {}
	assert(nthreads>0);

	typedef PsimagLite::Concurrency ConcurrencyType;
	ConcurrencyType concurrency(&argc,&argv,nthreads);

	SizeType dof = 2; // spin up and down

	GeometryLibraryType geometry(geometryParams,GeometryLibraryType::DECAY_0);

	std::cerr<<geometry;

	EngineType engine(geometry.matrix(),
	                  geometryParams.outputFile,
	                  dof,
	                  EngineType::VERBOSE_YES);
	PsimagLite::Vector<SizeType>::Type ne(dof,electronsUp);
	bool debug = false;
	bool verbose = false;

	RealType sum = 0;
	for (SizeType i=0;i<electronsUp;i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";

	typedef FreeFermions::ParallelDecay<RealType,FieldType,EngineType> ParallelDecayType;
	ParallelDecayType::HilbertStateType gs(engine,ne,debug);
	SizeType sigma3 = 0;

	std::cout<<"#sites= ";
	for (SizeType i=0;i<sites.size();i++) std::cout<<sites[i]<<" ";
	std::cout<<"\n";

	GeometryLibraryType geometry2(geometryParams,GeometryLibraryType::DECAY_1);
	std::cerr<<geometry2;
	EngineType engine2(geometry2.matrix(),
	                   geometryParams.outputFile,
	                   dof,
	                   EngineType::VERBOSE_YES);


	typedef PsimagLite::Parallelizer<ParallelDecayType> ParallelizerType;
	ParallelizerType threadedDecay(PsimagLite::Concurrency::codeSectionParams);

	ParallelDecayType::DecayParamsType params(geometryParams.orbitals,
	                                          geometryParams.sites,
	                                          sites,
	                                          sigma3,
	                                          offset,
	                                          step,
	                                          verbose);
	ParallelDecayType helperDecay(engine2,params,gs,total);

	FieldType superdensity = helperDecay.calcSuperDensity();
	std::cout<<"#superdensity="<<superdensity<<"\n";

	threadedDecay.loopCreate(helperDecay);
	helperDecay.printTasks(std::cout);
}

