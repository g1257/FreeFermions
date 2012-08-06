

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
#include "Tokenizer.h"

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

void usage(const std::string& thisFile)
{
	std::cout<<thisFile<<": USAGE IS "<<thisFile<<" ";
	std::cout<<" -n sites -e electronsUp -g geometry,[leg,filename]\n";
}
	
void readPotential(std::vector<RealType>& v,const std::string& filename)
{
	std::vector<RealType> w;
	PsimagLite::IoSimple::In io(filename);
	try {
		io.read(w,"PotentialT");
	} catch (std::exception& e) {
		std::cerr<<"INFO: No PotentialT in file "<<filename<<"\n";
	}
	io.rewind();
	
	io.read(v,"potentialV");
	if (v.size()>w.size()) v.resize(w.size());
	for (size_t i=0;i<w.size();i++) v[i] += w[i];
}

void setMyGeometry(GeometryParamsType& geometryParams,const std::vector<std::string>& vstr)
{
	// default value
	geometryParams.type = GeometryLibraryType::CHAIN;
		
	if (vstr.size()<2) {
		// assume chain
		return;
	}

	std::string gName = vstr[0];
	if (gName == "chain") {
		throw std::runtime_error("setMyGeometry: -g chain takes no further arguments\n");
	}
	
	geometryParams.leg = atoi(vstr[1].c_str());
	
	if (gName == "ladder") {
		if (vstr.size()!=2) {
			usage("setMyGeometry");
			throw std::runtime_error("setMyGeometry: usage is: -g ladder leg \n");
		}
		geometryParams.type = GeometryLibraryType::LADDER;
		return;
	}
	
	if (vstr.size()!=3) {
			usage("setMyGeometry");
			throw std::runtime_error("setMyGeometry: usage is: -g {feas | ktwoniffour} leg filename\n");
	}

	geometryParams.filename = vstr[2];

	if (gName == "feas") {
		geometryParams.type = GeometryLibraryType::FEAS;
		return;
	}
	
	if (gName == "feas") {
		geometryParams.type = GeometryLibraryType::KTWONIFFOUR;
		return;
	}
}

int main(int argc,char* argv[])
{
	size_t n = 0;
	size_t electronsUp = 0;
	std::vector<RealType> v;
	GeometryParamsType geometryParams;
	std::vector<std::string> str;
	int opt = 0;
	
	geometryParams.type = GeometryLibraryType::CHAIN;
	
	while ((opt = getopt(argc, argv, "n:e:g:p")) != -1) {
		switch (opt) {
			case 'n':
				n = atoi(optarg);
				v.resize(n,0);
				geometryParams.sites = n;
				break;
			case 'e':
				electronsUp = atoi(optarg);
				break;
			case 'g':
				PsimagLite::tokenizer(optarg,str,",");
				setMyGeometry(geometryParams,str);
				break;
			case 'p':
				readPotential(v,optarg);
				break;
			default: /* '?' */
				usage("setMyGeometry");
				throw std::runtime_error("Wrong usage\n");
		}
	}
	if (n==0 || geometryParams.sites==0 || v.size()!=n) {
		usage("setMyGeometry");
		throw std::runtime_error("Wrong usage\n");
	}

	size_t dof = 1; // spinless

	GeometryLibraryType geometry(geometryParams);

 	geometry.addPotential(v);
	std::cerr<<geometry;
	ConcurrencyType concurrency(argc,argv);
	EngineType engine(geometry,concurrency,dof,true);
	std::vector<size_t> ne(dof,electronsUp); // n. of up (= n. of  down electrons)
	HilbertStateType gs(engine,ne);
	RealType sum = 0;
	for (size_t i=0;i<ne[0];i++) sum += engine.eigenvalue(i);
	std::cerr<<"Energy="<<dof*sum<<"\n";	
	size_t sigma = 0;
	//MatrixType cicj(n,n);
	size_t norb = (geometryParams.type == GeometryLibraryType::FEAS) ? 2 : 1;
	for (size_t orbital=0; orbital<norb; orbital++) {
		for (size_t site = 0; site<n ; site++) {
			OpNormalFactoryType opNormalFactory(engine);
			OperatorType& myOp = opNormalFactory(OperatorType::DESTRUCTION,site+orbital*n,sigma);
			for (size_t site2=0; site2<n; site2++) {
				HilbertStateType phi = gs;
				myOp.applyTo(phi);
				OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,site2+orbital*n,sigma);
				myOp2.applyTo(phi);
				std::cout<<scalarProduct(gs,phi)<<" ";
				//cicj(site,site2) += scalarProduct(gs,phi);
			}
			std::cout<<"\n";
		}
		std::cout<<"-------------------------------------------\n";
	}
}

