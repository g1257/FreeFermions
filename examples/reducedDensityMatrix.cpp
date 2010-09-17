

// TBW FIXME
#include "Includes.h"
#include "ReducedDensityMatrix.h"

int main(int argc,char* argv[])
{
	if (argc<3) throw std::runtime_error("Needs 2 or more argument(s)\n");
        size_t n = atoi(argv[1]); 
	size_t nup =  atoi(argv[2]); 
	typedef FreeFermions::ReducedDensityMatrix<double,double> ReducedDensityMatrixType;
	
	MatrixType t(2*n,2*n);
	size_t geometryType = GeometryLibraryType::CHAIN;
	size_t leg = GeometryLibraryType::OPTION_PERIODIC;
	if (argc>=4) {
		leg = atoi(argv[3]);
		if (leg>0) geometryType = GeometryLibraryType::LADDER;
	}
	GeometryLibraryType geometry(2*n,geometryType);
	geometry.setGeometry(t,leg);
	std::cerr<<t;
	
	size_t dof = 1;
	EngineType engine(t,dof,false);
	
	ReducedDensityMatrixType reducedDensityMatrix(engine,n,nup);
	
	std::vector<double> e(reducedDensityMatrix.rank());
	reducedDensityMatrix.diagonalize(e);
	std::cout<<"DensityMatrixEigenvalues:\n";
	std::cout<<e;
}
