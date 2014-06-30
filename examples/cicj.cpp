

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
#include "InputNg.h"
#include "InputCheck.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef ComplexType FieldType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::Matrix<FieldType> MatrixType;
typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
typedef FreeFermions::GeometryParameters<FieldType,InputNgType::Readable> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::HilbertState<OperatorType> HilbertStateType;
typedef OperatorType::FactoryType OpNormalFactoryType;


int main(int argc,char* argv[])
{
	int opt = 0;
	PsimagLite::String file("");
	bool energyOnly = false;

	while ((opt = getopt(argc, argv, "f:e")) != -1) {
		switch (opt) {
		case 'f':
			file=optarg;
			break;
		case 'e':
			energyOnly=true;
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

	SizeType dof = 1; // spinless

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;
	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc,&argv,npthreads);
	EngineType engine(geometry.matrix(),dof,true);
	PsimagLite::Vector<SizeType>::Type ne(dof,electronsUp); // n. of up (= n. of  down electrons)
	HilbertStateType gs(engine,ne);
	RealType sum = 0;
	for (SizeType i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";	
	
	if (energyOnly) return 0;

	SizeType sigma = 0;
	//MatrixType cicj(n,n);
	SizeType norb = 1;
	io.readline(norb,"Orbitals=");
	for (SizeType orbital=0; orbital<norb; orbital++) {
		for (SizeType site = 0; site<geometryParams.sites ; site++) {
			OpNormalFactoryType opNormalFactory(engine);
			OperatorType& myOp = opNormalFactory(OperatorType::DESTRUCTION,site+orbital*geometryParams.sites,sigma);
			for (SizeType site2=0; site2<geometryParams.sites; site2++) {
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

