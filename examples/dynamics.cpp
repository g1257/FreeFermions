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
#include "GeometryParameters.h"
#include "Concurrency.h"
#include "InputNg.h"
#include "InputCheck.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef ComplexType FieldType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
typedef FreeFermions::GeometryParameters<RealType,InputNgType::Readable> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::OneOverZminusH<EngineType> OneOverZminusHType;
typedef FreeFermions::DiagonalOperator<OneOverZminusHType> DiagonalOperatorType;
typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
typedef DiagonalOperatorType::FactoryType OpDiagonalFactoryType;
typedef OperatorType::FactoryType OpNormalFactoryType;

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
		throw std::runtime_error("setMyGeometry: -g chain takes no further arguments\n");
	}

	geometryParams.leg = atoi(vstr[1].c_str());

	if (gName == "ladder") {
		if (vstr.size()!=3) {
			usage("setMyGeometry");
			throw std::runtime_error("setMyGeometry: usage is: -g ladder,leg,isPeriodic \n");
		}
		geometryParams.type = GeometryLibraryType::LADDER;
		geometryParams.hopping.resize(2);
		geometryParams.hopping[0] =  geometryParams.hopping[1]  = 1.0;
		geometryParams.isPeriodic[GeometryParamsType::DIRECTION_Y] = (atoi(vstr[2].c_str())>0);
		return;
	}

	if (vstr.size()!=3) {
		usage("setMyGeometry");
		throw std::runtime_error("setMyGeometry: usage is: -g {feas | ktwoniffour} leg filename\n");
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

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	SizeType total=0;
	RealType offset = 0;
	RealType step = 0;
	SizeType dynType = DYN_TYPE_0;

	while ((opt = getopt(argc, argv, "f:t:o:i:d")) != -1) {
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
		case 'd':
			dynType = DYN_TYPE_1;
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

	SizeType dof = 1; // spinless

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;

	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc,&argv,npthreads);
	EngineType engine(geometry.matrix(),
	                  geometryParams.outputFile,
	                  dof,
	                  EngineType::VERBOSE_YES);

	PsimagLite::Vector<SizeType>::Type ne(dof,electronsUp); // 8 up and 8 down
	bool debug = false;
	HilbertStateType gs(engine,ne,debug);

	RealType Eg = 0;
	for (SizeType i=0;i<ne[0];i++) Eg += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*Eg<<"\n";

	SizeType sigma =0;
	OpNormalFactoryType opNormalFactory(engine);
	SizeType creatOrDest = (dynType == DYN_TYPE_1) ? OperatorType::CREATION :
	                                                 OperatorType::DESTRUCTION;

	OperatorType& opKet = opNormalFactory(creatOrDest,sites[0],sigma);
	HilbertStateType phiKet = gs;
	opKet.applyTo(phiKet);

	if (sites.size() == 1) {
		sites.resize(2);
		sites[1] = sites[0];
	}

	OperatorType& opBra = opNormalFactory(creatOrDest,sites[1],sigma);
	HilbertStateType phiBra = gs;
	opBra.applyTo(phiBra);

	FieldType density = scalarProduct(phiBra,phiKet);
	std::cerr<<"density="<<density<<"\n";

	std::cout<<"#site="<<sites[0]<<"\n";
	std::cout<<"#site2="<<sites[1]<<"\n";

	RealType epsilon = 1e-2;
	for (SizeType it = 0; it<total; it++) {
		OpDiagonalFactoryType opDiagonalFactory(engine);
		RealType omega = it * step + offset;
		ComplexType z = ComplexType(omega,epsilon);
		int sign = (dynType== DYN_TYPE_1) ? -1 : 1;
		OneOverZminusHType eih(z,sign,Eg,engine);
		DiagonalOperatorType& eihOp = opDiagonalFactory(eih);
		HilbertStateType phi3 = phiKet;
		eihOp.applyTo(phi3);
		FieldType tmpC = scalarProduct(phiBra,phi3);

		std::cout<<omega<<" "<<PsimagLite::imag(tmpC)<<" "<<PsimagLite::real(tmpC)<<"\n";
	}
}

