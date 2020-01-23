// SAmple of how to use FreeFermions core engine to calculate
// <\sum_i W_i c^\dagger_i >
#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"
#include "GeometryParameters.h"
#include "Concurrency.h"
#include "InputNg.h"
#include "InputCheck.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::Matrix<ComplexType> MatrixType;
typedef PsimagLite::Matrix<MatrixType*> MatrixOfMatrixPointersType;
typedef PsimagLite::Vector<MatrixType*>::Type VectorMatrixType;
typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
typedef PsimagLite::Vector<ComplexType>::Type VectorType;
typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
typedef FreeFermions::GeometryParameters<RealType,InputNgType::Readable> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,ComplexType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;
typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
typedef PsimagLite::Vector<HilbertStateType*>::Type VectorHilbertStateType;
typedef DiagonalOperatorType::FactoryType OpDiagonalFactoryType;
typedef OperatorType::FactoryType OpNormalFactoryType;

void computeOneBucket(HilbertStateType& phi,
                      OpNormalFactoryType& opNormalFactory,
                      DiagonalOperatorType* eihOp,
                      SizeType physicalSite,
                      SizeType totalSites)
{
	SizeType sigma = 0;
	OperatorType& opC = opNormalFactory(OperatorType::DESTRUCTION,
	                                    physicalSite + 1*totalSites,
	                                    sigma);
	opC.applyTo(phi);

	OperatorType& opCdagger = opNormalFactory(OperatorType::CREATION,
	                                          physicalSite + 0*totalSites,
	                                          sigma);
	opCdagger.applyTo(phi);

	if (eihOp) eihOp->applyTo(phi);
}

ComplexType bucketFinal(const VectorHilbertStateType& buckets, const VectorType& weights)
{
	SizeType total = buckets.size();
	MatrixType m(total,total);
	for (SizeType i = 0; i < total; ++i) {
		for (SizeType j = i; j < total; ++j) {
			m(i, j) = scalarProduct(*buckets[i],*buckets[j]);
		}
	}

	for (SizeType i = 0; i < total; ++i) delete buckets[i];

	ComplexType result = 0;
	for (SizeType i = 0; i < total; ++i) {
		for (SizeType j = 0; j < total; ++j) {
			ComplexType value = (j > i) ? PsimagLite::conj(m(i,j)) : m(j,i);
			result += PsimagLite::conj(weights[i])*weights[j]*value;
		}
	}

	return result;
}

ComplexType doOneTime(RealType time,
                      const HilbertStateType& gs,
                      const EngineType& engine,
                      OpNormalFactoryType& opNormalFactory,
                      OpDiagonalFactoryType& opDiagonalFactory,
                      OperatorType& opCp,
                      const VectorSizeType& sites,
                      const VectorType& weights,
                      SizeType totalSites)
{
	EtoTheIhTimeType eih(time, engine, 0);
	DiagonalOperatorType& eihOp = opDiagonalFactory(eih);
	SizeType total = sites.size();
	VectorHilbertStateType buckets(total);
	for (SizeType i = 0; i < total; ++i) {
		const SizeType site = sites[i];
		buckets[i] = new HilbertStateType(gs);
		computeOneBucket(*buckets[i], opNormalFactory, &eihOp, site, totalSites);
		opCp.applyTo(*buckets[i]);
	}

	return bucketFinal(buckets, weights);
}

void doOneSite(MatrixType& values,
               SizeType site3,
               const HilbertStateType& gs,
               const EngineType& engine,
               const VectorSizeType& sites,
               const VectorType& weights,
               RealType offset,
               RealType step,
               SizeType totalSites)
{
	SizeType sigma = 0;
	OpNormalFactoryType opNormalFactory(engine);
	OperatorType& opCp = opNormalFactory(OperatorType::DESTRUCTION, site3, sigma);
	SizeType total = values.n_col();

	for (SizeType it = 0; it < total; ++it) {
		OpDiagonalFactoryType opDiagonalFactory(engine);
		const RealType time = it * step + offset;
		const ComplexType result =  doOneTime(time,
		                                      gs,
		                                      engine,
		                                      opNormalFactory,
		                                      opDiagonalFactory,
		                                      opCp,
		                                      sites,
		                                      weights,
		                                      totalSites);
		values(site3, it) = result;
	}
}

ComplexType superDensity(const HilbertStateType& gs,
                         const EngineType& engine,
                         const VectorSizeType& sites,
                         const VectorType& weights,
                         SizeType totalSites)
{
	SizeType total = sites.size();
	OpNormalFactoryType opNormalFactory(engine);
	VectorHilbertStateType buckets(total,0);
	for (SizeType i = 0; i < sites.size(); ++i) {
		SizeType site = sites[i];
		buckets[i] = new HilbertStateType(gs);
		computeOneBucket(*buckets[i], opNormalFactory, nullptr, site, totalSites);
	}

	return bucketFinal(buckets, weights);
}

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	SizeType total=0;
	RealType offset = 0;
	RealType step = 0;
	SizeType site3 = 0;
	bool allSites = true;

	while ((opt = getopt(argc, argv, "f:t:o:i:p:")) != -1) {
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
		case 'p':
			site3 = atoi(optarg);
			allSites = false;
			break;
		default: /* '?' */
			throw std::runtime_error("Wrong usage\n");
		}
	}

	if (file=="" || total == 0) {
		throw std::runtime_error("Wrong usage\n");
	}

	FreeFermions::InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(file,inputCheck);
	InputNgType::Readable io(ioWriteable);

	GeometryParamsType geometryParams(io);
	SizeType electronsUp = GeometryParamsType::readElectrons(io,geometryParams.sites);
	SizeType orbitals = 1;
	try {
		io.readline(orbitals,"Orbitals=");
	} catch (std::exception&) {}

	PsimagLite::Vector<SizeType>::Type sites;
	GeometryParamsType::readVector(sites,file,"TSPSites");
	VectorType weights;
	GeometryParamsType::readVector(weights,file,"TSPOperatorMultiplier");

	if (weights.size() != sites.size()) {
		throw PsimagLite::RuntimeError("sites.size() != weights.size()\n");
	}

	SizeType dof = 1; // spinless

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;

	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc, &argv, npthreads);
	EngineType engine(geometry.matrix(),
	                  geometryParams.outputFile,
	                  dof,
	                  EngineType::VerboseEnum::YES);

	PsimagLite::Vector<SizeType>::Type ne(dof,electronsUp);
	bool debug = false;
	HilbertStateType gs(engine,ne,debug);

	MatrixType values(geometryParams.sites*orbitals, total);
	SizeType sitesUpTo = values.n_row();

	ComplexType sd = superDensity(gs, engine, sites, weights, geometryParams.sites);

	std::cout<<"#Superdensity="<<sd<<"\n";

	if (!allSites) {
		doOneSite(values, site3, gs, engine, sites, weights, offset, step, geometryParams.sites);
		sitesUpTo = 1;
	} else {
		for (SizeType i = 0; i < orbitals*geometryParams.sites; ++i)
			doOneSite(values, i, gs, engine, sites, weights, offset, step, geometryParams.sites);
	}

	for (SizeType it = 0; it < total; ++it) {
		RealType time = it * step + offset;
		std::cout<<time<<"    ";
		for (SizeType i = 0; i < sitesUpTo; ++i) {
			SizeType thisSite = (allSites) ? i : site3;
			const RealType val = PsimagLite::real(values(thisSite, it)/sd);
			std::cout<<val<<" ";
		}

		std::cout<<"\n";
	}
}

