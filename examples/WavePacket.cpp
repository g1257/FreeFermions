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
                      SizeType siteMixed)
{
	SizeType sigma = 0;
	OperatorType& opCsite = opNormalFactory(OperatorType::DESTRUCTION,siteMixed,sigma);
	opCsite.applyTo(phi);
}

void computeOneBucket(HilbertStateType& phi,
                      OpNormalFactoryType& opNormalFactory,
                      DiagonalOperatorType& eihOp,
                      OperatorType& opCp,
                      SizeType siteMixed)
{
	SizeType sigma = 0;
	OperatorType& opCsite = opNormalFactory(OperatorType::DESTRUCTION,siteMixed,sigma);
	opCsite.applyTo(phi);
	eihOp.applyTo(phi);
	opCp.applyTo(phi);
}

void bucketFinal(MatrixType& result,
                 const VectorHilbertStateType& buckets,
                 const VectorType& weights,
                 SizeType orbitals)
{
	SizeType total = buckets.size();
	MatrixType m(total,total);
	for (SizeType i = 0; i < total; ++i) {
		for (SizeType j = i; j < total; ++j) {
			m(i,j) = scalarProduct(*buckets[i],*buckets[j]);
		}
	}

	for (SizeType i = 0; i < total; ++i) delete buckets[i];

	for (SizeType i = 0; i < total; ++i) {
		SizeType iSite = static_cast<SizeType>(i/orbitals);
		SizeType iOrb = i % orbitals;
		for (SizeType j = 0; j < total; ++j) {
			SizeType jOrb = j % orbitals;
			SizeType jSite = static_cast<SizeType>(j/orbitals);
			ComplexType value = (j > i) ? PsimagLite::conj(m(i,j)) : m(j,i);
			result(iOrb,jOrb) += PsimagLite::conj(weights[iSite])*weights[jSite]*value;
		}
	}
}

void doOneTime(MatrixType& result,
               RealType time,
               const HilbertStateType& gs,
               const EngineType& engine,
               OpNormalFactoryType& opNormalFactory,
               OpDiagonalFactoryType& opDiagonalFactory,
               OperatorType& opCp,
               const VectorSizeType& sites,
               const VectorType& weights,
               SizeType totalSites,
               SizeType orbitals)
{
	EtoTheIhTimeType eih(time,engine,0);
	DiagonalOperatorType& eihOp = opDiagonalFactory(eih);
	SizeType total = sites.size()*orbitals;
	VectorHilbertStateType buckets(total,0);
	for (SizeType i = 0; i < sites.size(); ++i) {
		SizeType site = sites[i];
		for (SizeType orb = 0; orb < orbitals; ++orb) {
			SizeType mixI = orb + i*orbitals;
			buckets[mixI] = new HilbertStateType(gs);
			SizeType siteMixed = site + orb*totalSites;
			computeOneBucket(*buckets[mixI],opNormalFactory,eihOp,opCp,siteMixed);
		}
	}

	return bucketFinal(result,buckets,weights,orbitals);
}

void doOneSite(MatrixOfMatrixPointersType& values,
               VectorMatrixType& garbage,
               SizeType site3,
               const HilbertStateType& gs,
               const EngineType& engine,
               const VectorSizeType& sites,
               const VectorType& weights,
               SizeType orbitals,
               RealType offset,
               RealType step)
{
	SizeType sigma =0;
	OpNormalFactoryType opNormalFactory(engine);
	OperatorType& opCp = opNormalFactory(OperatorType::DESTRUCTION,site3,sigma);
	SizeType totalSites = static_cast<SizeType>(values.n_row()/orbitals);
	SizeType total = values.n_col();

	for (SizeType it = 0; it < total; ++it) {
		OpDiagonalFactoryType opDiagonalFactory(engine);
		RealType time = it * step + offset;
		MatrixType *result = new MatrixType(orbitals,orbitals);
		garbage.push_back(result);
		doOneTime(*result,time,gs,engine,opNormalFactory,opDiagonalFactory,
		          opCp,sites,weights,totalSites,orbitals);
		values(site3,it) = result;
	}
}

void superDensity(MatrixType& sd,
                  const HilbertStateType& gs,
                  const EngineType& engine,
                  const VectorSizeType& sites,
                  const VectorType& weights,
                  SizeType totalSites,
                  SizeType orbitals)
{
	SizeType total = sites.size()*orbitals;
	OpNormalFactoryType opNormalFactory(engine);
	VectorHilbertStateType buckets(total,0);
	for (SizeType i = 0; i < sites.size(); ++i) {
		SizeType site = sites[i];
		for (SizeType orb = 0; orb < orbitals; ++orb) {
			SizeType mixI = orb + i*orbitals;
			buckets[mixI] = new HilbertStateType(gs);
			SizeType siteMixed = site + orb*totalSites;
			computeOneBucket(*buckets[mixI],opNormalFactory,siteMixed);
		}
	}

	return bucketFinal(sd,buckets,weights,orbitals);
}

void printValue(const MatrixOfMatrixPointersType& values,
                SizeType thisSite,
                SizeType it,
                int orb,
                SizeType orbitals)
{
	ComplexType sum = 0.0;
	if (orb < 0) {
		for (SizeType i = 0; i < orbitals; ++i)
			for (SizeType j = 0; j < orbitals; ++j)
				sum += values(thisSite,it)->operator()(i,j);
	} else {
		sum = values(thisSite,it)->operator()(orb,orb);
	}

	std::cout<<sum<<" ";

}

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	SizeType total=0;
	RealType offset = 0;
	RealType step = 0;
	SizeType site3 = 0;
	int option = 0;
	bool allSites = true;

	while ((opt = getopt(argc, argv, "f:t:o:i:p:O:")) != -1) {
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
		case 'O':
			option = atoi(optarg);
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
	ConcurrencyType concurrency(&argc,&argv,npthreads);
	EngineType engine(geometry.matrix(),dof,true);

	PsimagLite::Vector<SizeType>::Type ne(dof,electronsUp);
	bool debug = false;
	HilbertStateType gs(engine,ne,debug);

	VectorMatrixType garbage;
	MatrixOfMatrixPointersType values(orbitals*geometryParams.sites,total);
	SizeType sitesUpTo = values.n_row();

	MatrixType sd(orbitals,orbitals);
	superDensity(sd,gs,engine,sites,weights,geometryParams.sites,orbitals);

	std::cout<<"#Superdensity=\n";
	std::cout<<sd;

	ComplexType sum = 0.0;
	for (SizeType i = 0; i < sd.n_row(); ++i)
		for (SizeType j = 0; j < sd.n_col(); ++j)
			sum += sd(i,j);
	std::cout<<"#SuperdensityAllOrbitals="<<sum<<"\n";

	if (!allSites) {
		doOneSite(values,garbage,site3,gs,engine,sites,weights,orbitals,offset,step);
		sitesUpTo = 1;
	} else {
		for (SizeType i = 0; i < orbitals*geometryParams.sites; ++i)
			doOneSite(values,garbage,i,gs,engine,sites,weights,orbitals,offset,step);
	}

	for (SizeType it = 0; it < total; ++it) {
		RealType time = it * step + offset;
		std::cout<<time<<" ";
		for (SizeType i = 0; i < sitesUpTo; ++i) {
			SizeType thisSite = (allSites) ? i : site3;
			printValue(values,thisSite,it,option,orbitals);
		}

		std::cout<<"\n";
	}

	for (SizeType i = 0; i < garbage.size(); ++i) {
		delete garbage[i];
		garbage[i] = 0;
	}
}

