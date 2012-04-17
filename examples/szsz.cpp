

// SAmple of how to use FreeFermions core engine to calculate
// <sz sz >
#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
#include "TypeToString.h"
#include "CreationOrDestructionOp.h"
#include "HilbertState.h"
#include "GeometryParameters.h"

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
typedef OperatorType::FactoryType OpNormalFactoryType;

enum {SPIN_UP,SPIN_DOWN};

int main(int argc,char* argv[])
{
	int argce = 3;
	size_t whatGeometry = GeometryLibraryType::FEAS; //CHAIN; // FEAS; //CHAIN; // KTWONIFFOUR;
	std::string s = "Needs " + ttos(argce) + " argument(s)\n";
	if (argc<argce) throw std::runtime_error(s.c_str());
	size_t n = atoi(argv[1]); // n. of  sites
	size_t dof = 2; // spin
	GeometryParamsType geometryParams;
	geometryParams.sites = n;
	geometryParams.type =whatGeometry;
	if (whatGeometry==GeometryLibraryType::LADDER || whatGeometry==GeometryLibraryType::FEAS)
		geometryParams.leg = 2;
	if (whatGeometry==GeometryLibraryType::FEAS || whatGeometry==GeometryLibraryType::KTWONIFFOUR)
		geometryParams.filename = argv[3];

// 	MatrixType t(4,4);
// 	std::vector<RealType> tHop(3);
// 	tHop[0] = 0.1;
// 	tHop[1] = 0.5; 
// 	tHop[2] = 1.3;
// 	t(0,1) = t(1,0) = tHop[0];
// 	t(1,2) = t(2,1) = tHop[1];
// 	t(2,3) = t(3,2) = tHop[2];
	GeometryLibraryType geometry(geometryParams);
	//GeometryLibraryType geometry(n,GeometryLibraryType::CHAIN);
	//geometry.setGeometry(t,GeometryLibraryType::OPTION_PERIODIC);
	

 	/* std::vector<RealType> w;
	PsimagLite::IoSimple::In io(argv[3]);
	try {
		io.read(w,"PotentialT");
	} catch (std::exception& e) {
		std::cerr<<"No PotentialT in file "<<argv[3]<<"\n";
	}
	io.rewind();
	std::vector<RealType> v;
	io.read(v,"potentialV");
	for (size_t i=0;i<v.size();i++) v[i] += w[i];

 	geometry.addPotential(v);*/
	std::cerr<<geometry;
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(geometry,concurrency,dof,true);
	std::vector<size_t> ne(dof,atoi(argv[2])); // n. of up (= n. of  down electrons)
	HilbertStateType gs(engine,ne);
	RealType sum = 0;
	for (size_t i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";	
	std::cout.precision(20);
	//MatrixType cicj(n,n);
	size_t norb = (whatGeometry == GeometryLibraryType::FEAS) ? 2 : 1;
	for (size_t orbital1=0; orbital1<norb; orbital1++) {
		for (size_t orbital2=0; orbital2<norb; orbital2++) {
			for (size_t site = 0; site<n ; site++) {
				OpNormalFactoryType opNormalFactory(engine);
				OperatorType& myOp1 = opNormalFactory(OperatorType::DESTRUCTION,site+orbital1*n,SPIN_UP);
				OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,site+orbital1*n,SPIN_UP);
				OperatorType& myOp3 = opNormalFactory(OperatorType::DESTRUCTION,site+orbital1*n,SPIN_DOWN);
				OperatorType& myOp4 = opNormalFactory(OperatorType::CREATION,site+orbital1*n,SPIN_DOWN);
				HilbertStateType phi1 = gs;
				myOp1.applyTo(phi1);
				myOp2.applyTo(phi1);
				HilbertStateType phi2 = gs;
				myOp3.applyTo(phi2);
				myOp4.applyTo(phi2);
				for (size_t site2=0; site2<n; site2++) {
					OperatorType& myOp5 = opNormalFactory(OperatorType::DESTRUCTION,site2+orbital2*n,SPIN_UP);
					OperatorType& myOp6 = opNormalFactory(OperatorType::CREATION,site2+orbital2*n,SPIN_UP);
					OperatorType& myOp7 = opNormalFactory(OperatorType::DESTRUCTION,site2+orbital2*n,SPIN_DOWN);
					OperatorType& myOp8 = opNormalFactory(OperatorType::CREATION,site2+orbital2*n,SPIN_DOWN);
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
					//cicj(site,site2) += scalarProduct(gs,phi);
				}
				std::cout<<"\n";
			}
			std::cout<<"-------------------------------------------\n";
		}
	}
}

