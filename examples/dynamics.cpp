

// SAmple of how to use FreeFermions core engine to calculate
// <A_i 1/(z-H) A_j>
#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "OneOverZminusH.h"
#include "DiagonalOperator.h"
#include "Tokenizer.h"
#include "GeometryParameters.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef ComplexType FieldType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::GeometryParameters<RealType> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType,ConcurrencyType> EngineType;
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

void setMyGeometry(GeometryParamsType& geometryParams,const PsimagLite::Vector<PsimagLite::String>::Type& vstr)
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
	size_t total=0;
	RealType offset = 0;
	RealType step = 0;
	size_t dynType = DYN_TYPE_0;

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

	GeometryParamsType geometryParams(file);
	size_t electronsUp = GeometryParamsType::readElectrons(file,geometryParams.sites);
	PsimagLite::Vector<size_t>::Type sites;
	GeometryParamsType::readVector(sites,file,"TSPSites");

	size_t dof = 1; // spinless

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;
	
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(geometry,concurrency,dof,true);

	PsimagLite::Vector<size_t>::Type ne(dof,electronsUp); // 8 up and 8 down
	bool debug = false;
	HilbertStateType gs(engine,ne,debug);

	RealType Eg = 0;
	for (size_t i=0;i<ne[0];i++) Eg += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*Eg<<"\n";

	size_t sigma =0;
	OpNormalFactoryType opNormalFactory(engine);
	size_t creatOrDest = (dynType == DYN_TYPE_1) ? OperatorType::CREATION : OperatorType::DESTRUCTION;
	OperatorType& myOp = opNormalFactory(creatOrDest,sites[0],sigma);

	HilbertStateType phi2 = gs;
	myOp.applyTo(phi2);
	
	FieldType density = scalarProduct(phi2,phi2);
	std::cerr<<"density="<<density<<"\n";
	
	std::cout<<"#site="<<sites[0]<<"\n";
	std::cout<<"#site2="<<sites[1]<<"\n";

	RealType epsilon = 1e-2;
	for (size_t it = 0; it<total; it++) {
		OpDiagonalFactoryType opDiagonalFactory(engine);
		RealType omega = it * step + offset;
		ComplexType z = ComplexType(omega,epsilon);
		int sign = (dynType== DYN_TYPE_1) ? -1 : 1;
		OneOverZminusHType eih(z,sign,Eg,engine);
		DiagonalOperatorType& eihOp = opDiagonalFactory(eih);
		HilbertStateType phi3 = phi2;
		eihOp.applyTo(phi3);
		FieldType tmpC = scalarProduct(phi2,phi3);

		std::cout<<omega<<" "<<std::imag(tmpC)<<" "<<std::real(tmpC)<<"\n";
	}
}
