

// SAmple of how to use FreeFermions core engine to calculate
// <c^\dagger_i c_j >
#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "EtoTheBetaH.h"
#include "DiagonalOperator.h"
#include "Tokenizer.h"
#include "LibraryOperator.h"
#include "Combinations.h"
#include "GeometryParameters.h"
#include "Vector.h"
#include "Concurrency.h"

typedef double RealType;
typedef RealType FieldType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::GeometryParameters<RealType> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::EToTheBetaH<EngineType> EtoTheBetaHType;
typedef FreeFermions::DiagonalOperator<EtoTheBetaHType> DiagonalOperatorType;
typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
typedef DiagonalOperatorType::FactoryType OpDiagonalFactoryType;
typedef OperatorType::FactoryType OpNormalFactoryType;
typedef FreeFermions::LibraryOperator<OperatorType> LibraryOperatorType;
typedef LibraryOperatorType::FactoryType OpLibFactoryType;

RealType fermi(const RealType& x)
{
	return 1.0/(1.0+exp(x));
}

int main(int argc,char *argv[])
{
	int opt;
	PsimagLite::String file("");
	RealType beta = 0;
	RealType mu=0;

	while ((opt = getopt(argc, argv, "f:m:b:")) != -1) {
		switch (opt) {
		case 'f':
			file=optarg;
			break;
		case 'm':
			mu = atof(optarg);
			break;
		case 'b':
			beta = atof(optarg);
			break;
		default: /* '?' */
			throw std::runtime_error("Wrong usage\n");
		}
	}

	if (file=="") {
		throw std::runtime_error("Wrong usage\n");
	}

	GeometryParamsType geometryParams(file);

	size_t dof = 1; // spinless
	GeometryLibraryType geometry(geometryParams);

	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc,&argv,npthreads);
	EngineType engine(geometry,dof,true);

	size_t n = engine.size();

	std::cout<<geometry;
	std::cout<<"#beta="<<beta<<" mu="<<mu<<"\n";

	PsimagLite::Vector<RealType>::Type ni(n);
	RealType density = 0.0;
	RealType energy = 0.0;
	for (size_t i=0;i<n;i++) {
		for (size_t j=0;j<n;j++) {
			RealType value = 0;
			for (size_t k=0;k<n;k++) {
				value += std::conj(engine.eigenvector(i,k))*engine.eigenvector(j,k)*
						fermi(beta*(engine.eigenvalue(k)-mu));
				if (i == j && i == 0) {
					density += fermi(beta*(engine.eigenvalue(k)-mu));
					energy += engine.eigenvalue(k)*fermi(beta*(engine.eigenvalue(k)-mu));
				}
			}
			std::cout<<value<<" ";
			if (i == j) ni[i] = value;
		}
		std::cout<<"\n";
	}

	RealType sum1 = 0.0;
	RealType sum2 = 0.0;
	std::cout<<"#density = "<<density<<"\n";
	for (size_t i=0;i<ni.size();i++) {
		std::cout<<i<<" "<<ni[i]<<"\n";
		if (i<ni.size()/2) {
			sum1 += ni[i];
		} else {
			sum2 += ni[i];
		}
	}
	std::cout<<"#total density= "<<sum1<<" "<<sum2<<" "<<(sum1+sum2)<<"\n";
	std::cout<<"#energy= "<<energy<<"\n";
}
