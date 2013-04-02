

// Calculates <phi | phi>
// where 
// |phi> = c_{p up} exp(iHt) c_{i\sigma} nbar_{i \bar{sigma}} c^dagger_{j sigma'} n_{j \bar{sigma'}} |gs>

#include <cstdlib>
#include "unistd.h"
#include "Engine.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"
#include "LibraryOperator.h"
#include "Tokenizer.h" // in PsimagLite
#include "GeometryParameters.h"
#include "Range.h"
#include "DriverHelper.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef ComplexType FieldType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::GeometryParameters<RealType> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType,ConcurrencyType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;
typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
typedef FreeFermions::LibraryOperator<OperatorType> LibraryOperatorType;
typedef OperatorType::FactoryType OpNormalFactoryType;
typedef LibraryOperatorType::FactoryType OpLibFactoryType;
typedef DiagonalOperatorType::FactoryType OpDiagonalFactoryType;
typedef FreeFermions::DriverHelper<GeometryLibraryType> DriverHelperType;

FieldType calcSuperDensity(size_t site,
                             size_t site2,
                             const HilbertStateType& gs,
                             const EngineType& engine)
{
	HilbertStateType savedVector = gs;
	FieldType savedValue = 0;
	FieldType sum = 0;
	OpNormalFactoryType opNormalFactory(engine);
	OpLibFactoryType opLibFactory(engine);

	for (size_t sigma = 0;sigma<2;sigma++) {
		HilbertStateType phi = gs;

		LibraryOperatorType& myOp = opLibFactory(
		                         LibraryOperatorType::N,site,1-sigma);

		myOp.applyTo(phi);
		OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,
		                         site,sigma);

		myOp2.applyTo(phi);
		
		for (size_t sigma2 = 0;sigma2 < 2;sigma2++) {
			HilbertStateType phi3 = phi;
			LibraryOperatorType& myOp3 = opLibFactory(
			                         LibraryOperatorType::NBAR,site2,1-sigma2);
			myOp3.applyTo(phi3);

			OperatorType& myOp4 = opNormalFactory(
			                         OperatorType::DESTRUCTION,site2,sigma2);

			myOp4.applyTo(phi3);

			if (sigma ==0 && sigma2 ==0) savedVector = phi3;
			sum += scalarProduct(phi3,phi3);

			if (sigma ==1 && sigma2 ==1) {
				savedValue = scalarProduct(phi3,savedVector);
			}
		}
	}
	sum += 2*real(savedValue);
	//std::cerr<<"#sum2="<<scalarProduct(tmpV,tmpV)<<"\n";
	return sum;
}


int main(int argc,char *argv[])
{
	int opt;
	std::string file("");
	size_t total=0;
	RealType offset = 0;
	std::vector<size_t> sites;
	std::vector<std::string> str;

	RealType step = 0;
	while ((opt = getopt(argc, argv, "f:s:t:o:i:")) != -1) {
		switch (opt) {
		case 'f':
			file=optarg;
			break;
		case 's':
			str.clear();
			PsimagLite::tokenizer(optarg,str,",");
			for (size_t i=0;i<str.size();i++)
				sites.push_back(atoi(str[i].c_str()));
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
	size_t electronsUp = DriverHelperType::readLabel(file,"TargetElectronsUp=");
	size_t dof = 2; // spin up and down

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;
	
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(geometry,concurrency,dof,false);

	std::vector<size_t> ne(dof,electronsUp); // 8 up and 8 down
	bool debug = false;
	bool verbose = false;
	HilbertStateType gs(engine,ne,debug);
	
	RealType sum = 0;
	for (size_t i=0;i<electronsUp;i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";

	size_t sigma3 = 0;
	
	FieldType superdensity = calcSuperDensity(sites[0],sites[1],gs,engine);
	std::cout<<"#superdensity="<<superdensity<<"\n";

	std::cout<<"#sites= ";
	for (size_t i=0;i<sites.size();i++) std::cout<<sites[i]<<" ";
	std::cout<<"\n";

	PsimagLite::Range<ConcurrencyType> range(0,total,concurrency);
	
	while(!range.end()) {
		size_t it = range.index();

		OpNormalFactoryType opNormalFactory(engine);
		OpLibFactoryType opLibFactory(engine);
		OpDiagonalFactoryType opDiagonalFactory(engine);

		RealType time = it * step + offset;
		EtoTheIhTimeType eih(time,engine,0);
		DiagonalOperatorType& eihOp = opDiagonalFactory(eih);

		HilbertStateType savedVector = gs;
		FieldType savedValue = 0;
		FieldType sum = 0;

		for (size_t sigma = 0;sigma<2;sigma++) {
			HilbertStateType phi = gs;
			LibraryOperatorType& myOp = opLibFactory(
					                      LibraryOperatorType::N,sites[0],1-sigma);
			myOp.applyTo(phi);

			OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,
					                         sites[0],sigma);
			myOp2.applyTo(phi);

			for (size_t sigma2 = 0;sigma2 < 2;sigma2++) {
				HilbertStateType phi3 = phi;


				LibraryOperatorType& myOp3 = opLibFactory(
				                 LibraryOperatorType::NBAR,sites[1],1-sigma2);
				myOp3.applyTo(phi3);

				OperatorType& myOp4 = opNormalFactory(
				                 OperatorType::DESTRUCTION,sites[1],sigma2);
				myOp4.applyTo(phi3);

				if (verbose) std::cerr<<"Applying exp(iHt)\n";
				eihOp.applyTo(phi3);

				if (verbose) std::cerr<<"Applying c_p\n";
				OperatorType& myOp6 = opNormalFactory(
								  OperatorType::DESTRUCTION,sites[2],sigma3);
				myOp6.applyTo(phi3);

				if (verbose) std::cerr<<"Adding "<<sigma<<" "<<sigma2<<" "<<it<<"\n";

				if (sigma ==0 && sigma2 ==0) savedVector = phi3;
				if (sigma ==1 && sigma2 ==1) {
					savedValue = scalarProduct(phi3,savedVector);
				}
				sum += scalarProduct(phi3,phi3);
				if (verbose) std::cerr<<"Done with scalar product\n";
			}
		}
		sum += 2*real(savedValue);
		std::cout<<time<<" "<<real(sum)<<"\n";
		range.next();
	}
}
