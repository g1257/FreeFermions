

// SAmple of how to use FreeFermions core engine to calculate
// < n_i n_j >

#include "Engine.h"
#include "ObservableLibrary.h"
#include "GeometryLibrary.h"
#include "ConcurrencySerial.h"
#include "TypeToString.h"
//#include "EtoTheIhTime.h"
//#include "DiagonalOperator.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef std::vector<bool> LevelsType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::Engine<RealType,RealType,LevelsType,ConcurrencyType> EngineType;
typedef EngineType::HilbertVectorType HilbertVectorType;
typedef EngineType::FreeOperatorType FreeOperatorType;
typedef FreeFermions::GeometryLibrary<MatrixType> GeometryLibraryType;
typedef FreeFermions::ObservableLibrary<EngineType> ObservableLibraryType;
//typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
//typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;

int main(int argc,char* argv[])
{
	int argce = 3;
	std::string s = "Needs " + ttos(argce) + " argument(s)\n";
	if (argc!=argce) throw std::runtime_error(s.c_str());
	size_t n = atoi(argv[1]); // 16 sites
	size_t dof = 2; // spin up and down
	MatrixType t(n,n);
	//std::vector<MatrixType> t;
	GeometryLibraryType geometry(n,GeometryLibraryType::LADDER);
	//GeometryLibraryType geometry(n,GeometryLibraryType::CHAIN);
	//geometry.setGeometry(t,GeometryLibraryType::OPTION_PERIODIC);
	
	//GeometryLibraryType geometry(n,GeometryLibraryType::FEAS);
	//size_t leg = atoi(argv[3]);
	geometry.setGeometry(t,2); //,argv[4],leg); //,GeometryLibraryType::OPTION_PERIODIC);
	std::vector<RealType> tb(2);
	tb[0] = -0.672705035; tb[1] = -0.672705035;	
	geometry.bathify(t,tb);
	std::vector<RealType> pot(t.n_row(),0);
	for (size_t i=0;i<4;i++) pot[i] = -4;
	for (size_t i=4;i<pot.size();i++) pot[i] = (i&1) ? 1.84716619 : -1.84716619;
	geometry.addPotential(t,pot);
	//std::cerr<<t;
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(t,concurrency,dof,true);
	std::vector<size_t> ne(dof,atoi(argv[2])); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	RealType sum = 0;
	for (size_t i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";	
	//ObservableLibraryType library(engine);
	size_t sigma = 0;
	MatrixType cicj(n,n);
	for (size_t orbital=0; orbital<2; orbital++) {
		for (size_t site = 0; site<n ; site++) {
			HilbertVectorType phi = engine.newState();
			FreeOperatorType myOp = engine.newSimpleOperator("destruction",site+orbital*n,sigma);
			myOp.apply(phi,gs,FreeOperatorType::SIMPLIFY);
			for (size_t site2=0; site2<n; site2++) {
				HilbertVectorType phi2 = engine.newState();
				FreeOperatorType myOp2 = engine.newSimpleOperator("creation",site2+orbital*n,sigma);
				myOp2.apply(phi2,phi,FreeOperatorType::SIMPLIFY);
				std::cout<<scalarProduct(phi2,gs)<<" ";
				cicj(site,site2) += scalarProduct(phi2,gs);
			}
			std::cout<<"\n";
		}
		std::cout<<"-------------------------------------------\n";
	}
	RealType charge = 0;
	for (size_t site = 0; site<n ; site++) {
		for (size_t site2=0; site2<n; site2++) {
			std::cout<<cicj(site,site2)<<" ";
			if (site2==site) charge += cicj(site,site);
		}
		std::cout<<"\n";	
	}
	std::cout<<"Total charge for spin "<<sigma<<" = "<<charge<<"\n";
}
