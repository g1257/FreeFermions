

// SAmple of how to use FreeFermions core engine to calculate
// <c^\dagger_i c_j >
#include <cstdlib>
#include <unistd.h>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "GeometryParameters.h"
#include "Tokenizer.h"
#include "Concurrency.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef RealType FieldType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::GeometryParameters<RealType> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::HilbertState<OperatorType> HilbertStateType;
typedef OperatorType::FactoryType OpNormalFactoryType;


int main(int argc,char* argv[])
{
	int opt = 0;
	PsimagLite::String file("");

	while ((opt = getopt(argc, argv, "f:")) != -1) {
		switch (opt) {
		case 'f':
			file=optarg;
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

	size_t dof = 1; // spinless

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;
	size_t npthreads = 1;
	ConcurrencyType concurrency(&argc,&argv,npthreads);
	EngineType engine(geometry,dof,true);
	PsimagLite::Vector<size_t>::Type ne(dof,electronsUp); // n. of up (= n. of  down electrons)
	HilbertStateType gs(engine,ne);
	RealType sum = 0;
	for (size_t i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";	
	size_t sigma = 0;
	//MatrixType cicj(n,n);
	size_t norb = (geometryParams.type == GeometryLibraryType::FEAS || geometryParams.type == GeometryLibraryType::FEAS1D) ? geometryParams.orbitals : 1;
	for (size_t orbital=0; orbital<norb; orbital++) {
		for (size_t site = 0; site<geometryParams.sites ; site++) {
			OpNormalFactoryType opNormalFactory(engine);
			OperatorType& myOp = opNormalFactory(OperatorType::DESTRUCTION,site+orbital*geometryParams.sites,sigma);
			for (size_t site2=0; site2<geometryParams.sites; site2++) {
				HilbertStateType phi = gs;
				myOp.applyTo(phi);
				OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,site2+orbital*geometryParams.sites,sigma);
				myOp2.applyTo(phi);
				std::cout<<scalarProduct(gs,phi)<<" ";
				//cicj(site,site2) += scalarProduct(gs,phi);
			}
			std::cout<<"\n";
		}
		std::cout<<"-------------------------------------------\n";
	}
}

