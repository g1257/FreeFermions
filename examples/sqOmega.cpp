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

FieldType doOneOmegaOneSitePair(const EngineType& engine,
                                const HilbertStateType& gs,
                                OpLibFactoryType& opLibFactory,
                                SizeType site0,
                                SizeType site1,
                                RealType Eg,
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

		HilbertStateType phiKet = gs;
		nOpI.applyTo(phiKet);

		LibraryOperatorType& nOpJ = opLibFactory(LibraryOperatorType::N,
		                                         site1,
		                                         sigma1);
		HilbertStateType phiBra = gs;
		nOpJ.applyTo(phiBra);

		//FieldType density = scalarProduct(phiBra,phiKet);
		//std::cerr<<"density="<<density<<"\n";

		OpDiagonalFactoryType opDiagonalFactory(engine);

		ComplexType z = ComplexType(omega,epsilon);
		OneOverZminusHType eih(z,signForDen,Eg,engine);
		DiagonalOperatorType& eihOp = opDiagonalFactory(eih);
		HilbertStateType phi3 = phiKet;
		eihOp.applyTo(phi3);

		tmpC += sign*scalarProduct(phiBra,phi3);
	}

	return tmpC;
}

void doOneOmega(const EngineType& engine,
                const HilbertStateType& gs,
                RealType Eg,
                SizeType totalSites,
                SizeType centralSite,
                RealType omega)
{
	OpLibFactoryType opLibFactory(engine);
	SizeType site0 = centralSite;
	std::cout<<omega<<" ";
	for (SizeType site1 = 0; site1 < totalSites; ++site1) {
		FieldType val = doOneOmegaOneSitePair(engine,
		                                      gs,
		                                      opLibFactory,
		                                      site0,
		                                      site1,
		                                      Eg,
		                                      omega);
		std::cout<<std::real(val)<<" "<<std::imag(val)<<" ";
	}

	std::cout<<"\n";
}

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
	std::cout<<"#############\n";

	for (SizeType it = 0; it<total; it++) {
		RealType omega = it * step + offset;
		doOneOmega(engine,
		           gs,
		           Eg,
		           geometryParams.sites,
		           centralSite,
		           omega);
	}
}

