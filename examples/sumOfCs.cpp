

// SAmple of how to use FreeFermions core engine to calculate
// <c^\dagger_i c_j >
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


// <phi| n_p | phi>
RealType phiNpPhi(OpNormalFactoryType& opNormalFactory,const HilbertStateType& gs,size_t siteP,const std::vector<size_t>& sites,size_t sigma)
{
	RealType sum = 0;
	OperatorType& cdaggerP = opNormalFactory(OperatorType::CREATION,siteP,sigma);
	OperatorType& cP = opNormalFactory(OperatorType::DESTRUCTION,siteP,sigma);
	for (size_t i = 0;i<sites.size();i++) {
		HilbertStateType backvector = gs;
		OperatorType& cdaggerI = opNormalFactory(OperatorType::CREATION,sites[i],sigma);
		cdaggerI.applyTo(backvector);
		cP.applyTo(backvector);
		cdaggerP.applyTo(backvector);
		for (size_t j = 0;j<sites.size();j++) {
				HilbertStateType tmp = gs; 
				OperatorType& cdaggerJ = opNormalFactory(OperatorType::CREATION,sites[j],sigma);
				cdaggerJ.applyTo(tmp);
				sum += scalarProduct(backvector,tmp);
		}
	}
	return sum;
}

// <phi| n_p | phi>
RealType phiPhi(OpNormalFactoryType& opNormalFactory,const HilbertStateType& gs,const std::vector<size_t>& sites,size_t sigma)
{
	RealType sum = 0;
	for (size_t i = 0;i<sites.size();i++) {
		HilbertStateType backvector = gs;
		OperatorType& cdaggerI = opNormalFactory(OperatorType::CREATION,sites[i],sigma);
		cdaggerI.applyTo(backvector);
		
		for (size_t j = 0;j<sites.size();j++) {
				HilbertStateType tmp = gs; 
				OperatorType& cdaggerJ = opNormalFactory(OperatorType::CREATION,sites[j],sigma);
				cdaggerJ.applyTo(tmp);
				sum += scalarProduct(backvector,tmp);
		}
	}
	return sum;
}

int main(int argc,char* argv[])
{
	int argce = 3;
	size_t whatGeometry = GeometryLibraryType::CHAIN; // FEAS; //CHAIN; // KTWONIFFOUR;
	std::string s = "Needs " + ttos(argce) + " argument(s)\n";
	if (argc<argce) throw std::runtime_error(s.c_str());
	size_t n = atoi(argv[1]); // n. of  sites
	size_t dof = 1; // spinless
	GeometryParamsType geometryParams;
	geometryParams.sites = n;
	geometryParams.type =whatGeometry;
	if (whatGeometry==GeometryLibraryType::LADDER || whatGeometry==GeometryLibraryType::FEAS)
		geometryParams.leg = 2;
	if (whatGeometry==GeometryLibraryType::FEAS || whatGeometry==GeometryLibraryType::KTWONIFFOUR)
		geometryParams.filename = argv[3];

	GeometryLibraryType geometry(geometryParams);
	//GeometryLibraryType geometry(n,GeometryLibraryType::CHAIN);
	//geometry.setGeometry(t,GeometryLibraryType::OPTION_PERIODIC);
	

//  	std::vector<RealType> w;
// 	PsimagLite::IoSimple::In io(argv[3]);
// 	try {
// 		io.read(w,"PotentialT");
// 	} catch (std::exception& e) {
// 		std::cerr<<"No PotentialT in file "<<argv[3]<<"\n";
// 	}
// 	io.rewind();
// 	std::vector<RealType> v;
// 	io.read(v,"potentialV");
// 	for (size_t i=0;i<v.size();i++) v[i] += w[i];
// 
//  	geometry.addPotential(v);
	std::cerr<<geometry;
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(geometry,concurrency,dof,true);
	std::vector<size_t> ne(dof,atoi(argv[2])); // n. of up (= n. of  down electrons)
	HilbertStateType gs(engine,ne);
	RealType sum = 0;
	for (size_t i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";	
	size_t sigma = 0;
	
	OpNormalFactoryType opNormalFactory(engine);
	
	//MatrixType cicj(n,n);
	size_t norb = (whatGeometry == GeometryLibraryType::FEAS) ? 2 : 1;
	std::vector<size_t> sites(3);
	sites[0]=3; sites[1]=4; sites[2]=5;
	RealType denominator = phiPhi(opNormalFactory,gs,sites,sigma);
	std::cout<<"site\tvalue\tnumerator\tdenominator\n";
	for (size_t orbital=0; orbital<norb; orbital++) {
		RealType total = 0;
		for (size_t site = 0; site<n ; site++) {
			RealType numerator = phiNpPhi(opNormalFactory,gs,site,sites,sigma);
			RealType value = numerator/denominator;
			std::cout<<site<<" "<<value<<" "<<numerator<<" "<<denominator<<"\n";
			total += value;
		}
		std::cout<<"total value="<<total<<"         -------------------------------------------\n";
	}
}

