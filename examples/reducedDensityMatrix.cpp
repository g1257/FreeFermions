

// TBW FIXME
#include "Includes.h"
#include "ReducedDensityMatrix.h"

int main(int argc,char* argv[])
{
	if (argc!=3) throw std::runtime_error("Needs 2 argument(s)\n");
        size_t n = atoi(argv[1]); 
	size_t nup =  atoi(argv[2]); 
	typedef FreeFermions::ReducedDensityMatrix<double,double> ReducedDensityMatrixType;
	
	MatrixType t(n,n);
	//GeometryLibraryType geometry(n,GeometryLibraryType::LADDER);
	GeometryLibraryType geometry(n,GeometryLibraryType::CHAIN);
	geometry.setGeometry(t); //,GeometryLibraryType::OPTION_PERIODIC);
	std::cerr<<t;	
	
	size_t dof = 1;
	EngineType engine(t,dof,false);
	
	ReducedDensityMatrixType reducedDensityMatrix(engine,n,nup);
	
	std::vector<double> e(reducedDensityMatrix.rank());
	reducedDensityMatrix.diagonalize(e);
	std::cout<<"DensityMatrixEigenvalues:\n";
	std::cout<<e;
}
