

// SAmple of how to use FreeFermions core engine to calculate
// < \delta_{i,\gamma} \delta^\dagger_{j,\gamma'} >

#include <cstdlib>
#include <unistd.h>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "GeometryParameters.h"
#include "Tokenizer.h"
#include "LibraryOperator.h"
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
typedef FreeFermions::LibraryOperator<OperatorType> LibraryOperatorType;
typedef LibraryOperatorType::FactoryType OpLibFactoryType;

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

	size_t dof = 2; // spin up and down

	GeometryLibraryType geometry(geometryParams);

	std::cerr<<geometry;
	ConcurrencyType concurrency(argc,argv);

	EngineType engine(geometry,dof,true);
	PsimagLite::Vector<size_t>::Type ne(dof,electronsUp); // n. of up (= n. of  down electrons)
	HilbertStateType gs(engine,ne);
	RealType sum = 0;
	for (size_t i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";

	size_t n = geometryParams.sites;
	MatrixType m(n,n);
	FieldType sum2 = 0;
	size_t effectiveN  = n;
	size_t norb = (geometryParams.type == GeometryLibraryType::FEAS || geometryParams.type == GeometryLibraryType::FEAS1D) ? geometryParams.orbitals : 1;

	OpLibFactoryType opLibFactory(engine);
	for (size_t site = 0; site<effectiveN ; site++) {
		FieldType y = 0;
		for (size_t orb1 = 0;orb1<norb; orb1++) {
			HilbertStateType phi = gs;
			LibraryOperatorType& myOp = opLibFactory(LibraryOperatorType::DELTA,site+orb1*n,0);
			myOp.applyTo(phi);
			y += scalarProduct(phi,gs);
			for (size_t site2=0; site2<effectiveN; site2++) {
				for (size_t orb2=0;orb2<norb;orb2++) {
					HilbertStateType phi2 = gs;
					LibraryOperatorType& myOp2 = opLibFactory(LibraryOperatorType::DELTA,site2+orb2*n,0);
					myOp2.applyTo(phi2);
					FieldType x = scalarProduct(phi2,phi);
					std::cout<<x<<" ";
					std::cout.flush();
					m(site,site2) += x;
					sum2 += y*y;
				}
			}
			std::cout<<"\n";
		}
	}
}
