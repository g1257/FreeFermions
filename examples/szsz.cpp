

// SAmple of how to use FreeFermions core engine to calculate
// <sz sz >
#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "GeometryParameters.h"
#include "Concurrency.h"
#include "InputNg.h"
#include "InputCheck.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef RealType FieldType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
typedef FreeFermions::GeometryParameters<RealType,InputNgType::Readable> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::HilbertState<OperatorType> HilbertStateType;
typedef OperatorType::FactoryType OpNormalFactoryType;

enum {SPIN_UP,SPIN_DOWN};

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

	FreeFermions::InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(file,inputCheck);
	InputNgType::Readable io(ioWriteable);

	GeometryParamsType geometryParams(io);
	SizeType electronsUp = GeometryParamsType::readElectrons(io, geometryParams.sites);

	SizeType dof = 2; // spin up and down

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;
	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc, &argv, npthreads);
	EngineType engine(geometry.matrix(),
	                  geometryParams.outputFile,
	                  dof,
	                  EngineType::VERBOSE_YES);
	PsimagLite::Vector<SizeType>::Type ne(dof, electronsUp);
	HilbertStateType gs(engine,ne);
	RealType sum = 0;
	for (SizeType i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";
	SizeType n = geometryParams.sites;

	std::cout.precision(20);
	//MatrixType cicj(n,n);
	SizeType norb = (geometryParams.type == GeometryLibraryType::FEAS ||
	                 geometryParams.type == GeometryLibraryType::FEAS1D) ?
	            geometryParams.orbitals : 1;
	for (SizeType orbital1=0; orbital1<norb; orbital1++) {
		for (SizeType orbital2=0; orbital2<norb; orbital2++) {
			for (SizeType site = 0; site<n ; site++) {
				OpNormalFactoryType opNormalFactory(engine);
				OperatorType& myOp1 = opNormalFactory(OperatorType::DESTRUCTION,
				                                      site + orbital1*n,
				                                      SPIN_UP);
				OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,
				                                      site + orbital1*n,
				                                      SPIN_UP);
				OperatorType& myOp3 = opNormalFactory(OperatorType::DESTRUCTION,
				                                      site + orbital1*n,
				                                      SPIN_DOWN);
				OperatorType& myOp4 = opNormalFactory(OperatorType::CREATION,
				                                      site + orbital1*n,
				                                      SPIN_DOWN);
				HilbertStateType phi1 = gs;
				myOp1.applyTo(phi1);
				myOp2.applyTo(phi1);
				HilbertStateType phi2 = gs;
				myOp3.applyTo(phi2);
				myOp4.applyTo(phi2);
				for (SizeType site2=0; site2<n; site2++) {
					OperatorType& myOp5 = opNormalFactory(OperatorType::DESTRUCTION,
					                                      site2 + orbital2*n,
					                                      SPIN_UP);
					OperatorType& myOp6 = opNormalFactory(OperatorType::CREATION,
					                                      site2 + orbital2*n,
					                                      SPIN_UP);
					OperatorType& myOp7 = opNormalFactory(OperatorType::DESTRUCTION,
					                                      site2 + orbital2*n,
					                                      SPIN_DOWN);
					OperatorType& myOp8 = opNormalFactory(OperatorType::CREATION,
					                                      site2 + orbital2*n,
					                                      SPIN_DOWN);
					HilbertStateType phi3 = gs;
					myOp5.applyTo(phi3);
					myOp6.applyTo(phi3);
					HilbertStateType phi4 = gs;
					myOp7.applyTo(phi4);
					myOp8.applyTo(phi4);
					RealType x13 = scalarProduct(phi1,phi3);
					RealType x24 = scalarProduct(phi2,phi4);
					RealType x14 = scalarProduct(phi1,phi4);
					RealType x23 = scalarProduct(phi2,phi3);
					std::cout<<(x13+x24-x14-x23)<<" ";
				}

				std::cout<<"\n";
			}

			std::cout<<"-------------------------------------------\n";
		}
	}
}

