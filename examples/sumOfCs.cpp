

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
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"

typedef double RealType;
typedef std::complex<double> ComplexType;
typedef ComplexType FieldType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef FreeFermions::GeometryParameters<RealType> GeometryParamsType;
typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
typedef FreeFermions::Engine<RealType,FieldType,ConcurrencyType> EngineType;
typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;
typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
typedef DiagonalOperatorType::FactoryType OpDiagonalFactoryType;
typedef OperatorType::FactoryType OpNormalFactoryType;


enum {SPIN_UP,SPIN_DOWN};

// <phi| n_p | phi>
FieldType phiNpPhi(OpNormalFactoryType& opNormalFactory,const HilbertStateType& gs,size_t siteP,const std::vector<size_t>& sites,size_t sigma,DiagonalOperatorType& eihOp,const std::vector<ComplexType>& weights)
{
	FieldType sum = 0;
	OperatorType& cdaggerP = opNormalFactory(OperatorType::CREATION,siteP,sigma);
	OperatorType& cP = opNormalFactory(OperatorType::DESTRUCTION,siteP,sigma);
	for (size_t i = 0;i<sites.size();i++) {
		HilbertStateType backvector = gs;
		OperatorType& cdaggerI = opNormalFactory(OperatorType::CREATION,sites[i],SPIN_UP);
		cdaggerI.applyTo(backvector);
		eihOp.applyTo(backvector);

		cP.applyTo(backvector);
		cdaggerP.applyTo(backvector);

		for (size_t j = 0;j<sites.size();j++) {
				HilbertStateType tmp = gs; 
				OperatorType& cdaggerJ= opNormalFactory(OperatorType::CREATION,sites[j],SPIN_UP);
				cdaggerJ.applyTo(tmp);
				eihOp.applyTo(tmp);

				sum += scalarProduct(backvector,tmp)*weights[i]*weights[j];
		}
	}
	return sum;
}

// <phi| | phi>
FieldType phiPhi(OpNormalFactoryType& opNormalFactory,const HilbertStateType& gs,const std::vector<size_t>& sites,size_t sigma,DiagonalOperatorType& eihOp,const std::vector<ComplexType>& weights)
{
	FieldType sum = 0;
	for (size_t i = 0;i<sites.size();i++) {
		HilbertStateType backvector = gs;
		OperatorType& cdaggerI = opNormalFactory(OperatorType::CREATION,sites[i],SPIN_UP);
		cdaggerI.applyTo(backvector);
		eihOp.applyTo(backvector);
		for (size_t j = 0;j<sites.size();j++) {
				HilbertStateType tmp = gs; 
				OperatorType& cdaggerJ = opNormalFactory(OperatorType::CREATION,sites[j],SPIN_UP);
				cdaggerJ.applyTo(tmp);
				eihOp.applyTo(tmp);
				sum += scalarProduct(backvector,tmp)*weights[i]*weights[j];
		}
	}
	return sum;
}

int main(int argc,char* argv[])
{
	int argce = 5;
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
	RealType timeMax =atof(argv[3]);
	RealType timeStep=atof(argv[4]);

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
	
	OpNormalFactoryType opNormalFactory(engine);
	
	//MatrixType cicj(n,n);
	//size_t norb = (whatGeometry == GeometryLibraryType::FEAS) ? 2 : 1;
	std::vector<size_t> sites(n-2);
	std::vector<ComplexType> weights(n-2);
	for (size_t i=1;i<n-1;i++) {
		sites[i-1]=i;
		RealType tmp123 = (i-n/2)*(i-n/2)/8.;
		weights[i-1] = exp(-tmp123);
		std::cout<<"WEIGHT["<<(i-1)<<"]="<<weights[i-1]<<" "<<tmp123<<" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
	}

//	std::cout<<"site\tvalue\tnumerator\tdenominator\n";
	std::cout<<"time ";
	for (size_t site = 0; site<n ; site++) std::cout<<site<<" ";
	std::cout<<"\n";
	for (RealType time=0;time<timeMax;time+=timeStep)  {
		FieldType total = 0;
		OpDiagonalFactoryType opDiagonalFactory(engine);
		EtoTheIhTimeType eih(time,engine,0);
		DiagonalOperatorType& eihOp = opDiagonalFactory(eih);
		FieldType denominator = phiPhi(opNormalFactory,gs,sites,SPIN_UP,eihOp,weights);
		std::cout<<time<<" ";
		for (size_t site = 0; site<n ; site++) {
			FieldType numerator = phiNpPhi(opNormalFactory,gs,site,sites,SPIN_UP,eihOp,weights);
			FieldType value = numerator/denominator;
			//std::cout<<site<<" "<<value<<" "<<numerator<<" "<<denominator<<"\n";
			RealType valueReal = std::real(value);
			assert(fabs(std::imag(value))<1e-6);
			std::cout<<valueReal<<" ";
			total += value;
		}
		RealType totalReal = std::real(total);
		std::cout<<totalReal<<"\n";
	}
}

