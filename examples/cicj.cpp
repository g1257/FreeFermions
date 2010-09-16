

// SAmple of how to use FreeFermions core engine to calculate
// < n_i n_j >

#include "Includes.h"

int main(int argc,char* argv[])
{
	if (argc!=3) throw std::runtime_error("Needs 2 argument(s)\n");
	size_t n = atoi(argv[1]); // 16 sites
	size_t dof = 2; // spin up and down
	MatrixType t(n,n);
	//GeometryLibraryType geometry(n,GeometryLibraryType::LADDER);
	GeometryLibraryType geometry(n,GeometryLibraryType::CHAIN);
	geometry.setGeometry(t,GeometryLibraryType::OPTION_PERIODIC);
	std::cerr<<t;	
	EngineType engine(t,dof,true);
	std::vector<size_t> ne(dof,atoi(argv[2])); // 8 up and 8 down
	HilbertVectorType gs = engine.newGroundState(ne);
	RealType sum = 0;
	for (size_t i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";	
	//ObservableLibraryType library(engine);
	size_t sigma = 1;
	for (size_t site = 0; site<n ; site++) {
		HilbertVectorType phi = engine.newState();
		FreeOperatorType myOp = engine.newSimpleOperator("destruction",site,sigma);
		myOp.apply(phi,gs,FreeOperatorType::SIMPLIFY);
		for (size_t site2=0; site2<n; site2++) {
			HilbertVectorType phi2 = engine.newState();
			FreeOperatorType myOp2 = engine.newSimpleOperator("creation",site2,sigma);
			myOp2.apply(phi2,phi,FreeOperatorType::SIMPLIFY);
			std::cout<<scalarProduct(phi2,gs)<<" ";
		}
		std::cout<<"\n";
	}
}
