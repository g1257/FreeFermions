

// TBW FIXME
#include "Includes.h"
#include "ReducedDensityMatrix.h"

int main(int argc,char* argv[])
{
	if (argc!=3) throw std::runtime_error("Needs 2 argument(s)\n");
        size_t n = atoi(argv[1]); 
	size_t nup =  atoi(argv[2]); 
	typedef FreeFermions::ReducedDensityMatrix<double,double> ReducedDensityMatrixType;
	
	ReducedDensityMatrixType reducedDensityMatrix(n,nup);
	
	std::vector<double> e(reducedDensityMatrix.rank());
	reducedDensityMatrix.diagonalize(e);
	std::cout<<e;
}
