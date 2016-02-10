

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
#include "Tokenizer.h"
#include "GeometryParameters.h"
#include "Concurrency.h"
#include "InputNg.h"
#include "InputCheck.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::Matrix<ComplexType> MatrixType;
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
                      DiagonalOperatorType& eihOp,
                      OperatorType& opCp,
                      SizeType site)
{
	SizeType sigma = 0;
	OperatorType& opCsite = opNormalFactory(OperatorType::DESTRUCTION,site,sigma);
	opCsite.applyTo(phi);
	eihOp.applyTo(phi);
	opCp.applyTo(phi);
}

void doOneTime(RealType time,
               const HilbertStateType& gs,
               const EngineType& engine,
               OpNormalFactoryType& opNormalFactory,
               OpDiagonalFactoryType& opDiagonalFactory,
               OperatorType& opCp,
               const VectorSizeType& sites,
               const VectorType& weights)
{
	EtoTheIhTimeType eih(time,engine,0);
	DiagonalOperatorType& eihOp = opDiagonalFactory(eih);
	VectorHilbertStateType buckets(sites.size(),0);
	for (SizeType i = 0; i < sites.size(); ++i) {
		SizeType site = sites[i];
		buckets[i] = new HilbertStateType(gs);
		computeOneBucket(*buckets[i],opNormalFactory,eihOp,opCp,site);
	}

	MatrixType m(sites.size(),sites.size());
	for (SizeType i = 0; i < sites.size(); ++i) {
		for (SizeType j = i; j < sites.size(); ++j) {
			m(i,j) = scalarProduct(*buckets[i],*buckets[j]);
		}
	}

	for (SizeType i = 0; i < sites.size(); ++i) delete buckets[i];

	ComplexType sum = 0.0;
	for (SizeType i = 0; i < sites.size(); ++i) {
		for (SizeType j = 0; j < sites.size(); ++j) {
			ComplexType value = (j > i) ? std::conj(m(j,i)) : m(i,j);
			sum += std::conj(weights[i])*weights[j]*value;
		}
	}

	std::cout<<time<<" "<<sum<<"\n";
}

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	SizeType total=0;
	RealType offset = 0;
	RealType step = 0;
	SizeType site3 = 0;

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
	ConcurrencyType concurrency(&argc,&argv,npthreads);
	EngineType engine(geometry.matrix(),dof,true);

	PsimagLite::Vector<SizeType>::Type ne(dof,electronsUp);
	bool debug = false;
	HilbertStateType gs(engine,ne,debug);

	SizeType sigma =0;
	OpNormalFactoryType opNormalFactory(engine);
	OperatorType& opCp = opNormalFactory(OperatorType::DESTRUCTION,site3,sigma);

	for (SizeType it = 0; it<total; it++) {
		OpDiagonalFactoryType opDiagonalFactory(engine);
		RealType time = it * step + offset;
		doOneTime(time,gs,engine,opNormalFactory,opDiagonalFactory,opCp,sites,weights);
	}
}
