

// SAmple of how to use FreeFermions core engine to calculate
// <c^\dagger_i c_j >
#include <cstdlib>
#include <unistd.h>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "GeometryParameters.h"
#include "Tokenizer.h"
#include "DriverHelper.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef RealType FieldType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::GeometryParameters<RealType> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType,ConcurrencyType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::HilbertState<OperatorType> HilbertStateType;
typedef FreeFermions::DriverHelper<GeometryLibraryType> DriverHelperType;
typedef OperatorType::FactoryType OpNormalFactoryType;


int main(int argc,char* argv[])
{
	int opt = 0;
	std::string file("");

	while ((opt = getopt(argc, argv, "f:")) != -1) {
		switch (opt) {
		case 'f':
			file=optarg;
			break;
		default: /* '?' */
			DriverHelperType::usage(argv[0],"-n sites -e electronsUp -g geometry,[leg,filename]");
			throw std::runtime_error("Wrong usage\n");
		}
	}

	if (file=="") {
		DriverHelperType::usage(argv[0],"-f file -g geometry,[leg,filename]");
		throw std::runtime_error("Wrong usage\n");
	}

	GeometryParamsType geometryParams(file);
	size_t electronsUp = 	DriverHelperType::readLabel(file,"TargetElectronsUp=");

	std::vector<RealType> v;
	DriverHelperType::readPotential(v,file);

	size_t n = geometryParams.sites;
	if (v.size()==4*n) {
		v.resize(2*n);
	}
	if (v.size()==2*n && geometryParams.type == GeometryLibraryType::LADDER) {
		v.resize(n);
	}

	size_t dof = 1; // spinless

	GeometryLibraryType geometry(geometryParams);

 	if (geometryParams.type!=GeometryLibraryType::KTWONIFFOUR) {
		geometry.addPotential(v);
	}

	std::cerr<<geometry;
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(geometry,concurrency,dof,true);
	std::vector<size_t> ne(dof,electronsUp); // n. of up (= n. of  down electrons)
	HilbertStateType gs(engine,ne);
	RealType sum = 0;
	for (size_t i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";	
	size_t sigma = 0;
	//MatrixType cicj(n,n);
	size_t norb = (geometryParams.type == GeometryLibraryType::FEAS || geometryParams.type == GeometryLibraryType::FEAS1D) ? geometryParams.orbitals : 1;
	for (size_t orbital=0; orbital<norb; orbital++) {
		for (size_t site = 0; site<n ; site++) {
			OpNormalFactoryType opNormalFactory(engine);
			OperatorType& myOp = opNormalFactory(OperatorType::DESTRUCTION,site+orbital*n,sigma);
			for (size_t site2=0; site2<n; site2++) {
				HilbertStateType phi = gs;
				myOp.applyTo(phi);
				OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,site2+orbital*n,sigma);
				myOp2.applyTo(phi);
				std::cout<<scalarProduct(gs,phi)<<" ";
				//cicj(site,site2) += scalarProduct(gs,phi);
			}
			std::cout<<"\n";
		}
		std::cout<<"-------------------------------------------\n";
	}
}

