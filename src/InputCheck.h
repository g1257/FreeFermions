
#ifndef INPUT_CHECK_H
#define INPUT_CHECK_H
#include <vector>
#include <stdexcept>
#include "Options.h"

namespace FreeFermions {

class InputCheck {

	typedef PsimagLite::Options::Readable OptionsReadableType;

public:

	InputCheck() : optsReadable_(0) {}

	~InputCheck()
	{
		if (optsReadable_!=0) delete optsReadable_;
	}

	bool check(const PsimagLite::String& label,
	           const PsimagLite::Vector<PsimagLite::String>::Type& vec,
	           SizeType line) const
	{
		if (label=="JMVALUES") {
			if (vec.size()!=2) return error1("JMVALUES",line);
			return true;
		} else if (label=="RAW_MATRIX") {
			SizeType row = atoi(vec[0].c_str());
			SizeType col = atoi(vec[1].c_str());
			SizeType n = row*col;
			if (vec.size()!=n+2) return error1("RAW_MATRIX",line);
			return true;
		} else if (label=="Connectors") {
			return true;
		} else if (label=="MagneticField") {
			return true;
		} else if (label=="FiniteLoops") {
			SizeType n = atoi(vec[0].c_str());
			if (vec.size()!=3*n+1)  return error1("FiniteLoops",line);
			return true;
		}
		return false;
	}

	bool checkSimpleLabel(const PsimagLite::String& label,
	                      SizeType line) const
	{
		// FIXME: needs implementation
		return true;
	}

	void check(const PsimagLite::String& label,const PsimagLite::String& val,SizeType)
	{
		if (label!="SolverOptions") return;
		PsimagLite::Vector<PsimagLite::String>::Type registerOpts;

		registerOpts.push_back("restart");
		registerOpts.push_back("debugmatrix");
		registerOpts.push_back("test");
		registerOpts.push_back("useDavidson");
		registerOpts.push_back("verbose");
		registerOpts.push_back("nofiniteloops");
		registerOpts.push_back("nowft");
		registerOpts.push_back("inflate");
		registerOpts.push_back("none");
		registerOpts.push_back("twositedmrg");
		registerOpts.push_back("noloadwft");
		registerOpts.push_back("ChebyshevSolver");
		registerOpts.push_back("InternalProductStored");
		registerOpts.push_back("InternalProductKron");
		registerOpts.push_back("useSu2Symmetry");
		registerOpts.push_back("TimeStepTargetting");
		registerOpts.push_back("DynamicTargetting");
		registerOpts.push_back("AdaptiveDynamicTargetting");
		registerOpts.push_back("CorrectionVectorTargetting");
		registerOpts.push_back("CorrectionTargetting");
		registerOpts.push_back("MettsTargetting");

		PsimagLite::Options::Writeable optWriteable(registerOpts,
		                                            PsimagLite::Options::Writeable::PERMISSIVE);
		optsReadable_ = new  OptionsReadableType(optWriteable,val);
	}

	bool isSet(const PsimagLite::String& thisOption) const
	{
		return optsReadable_->isSet(thisOption);
	}

	void checkForThreads(SizeType nthreads) const
	{
		if (nthreads==1) return;

		PsimagLite::String message1(__FILE__);
		message1 += " FATAL: You are requesting nthreads>0 but you ";
		message1 += " did not compile with USE_PTHREADS enabled\n";
		message1 += " Either set Threads=1 in the input file ";
		message1 += "(you won't have threads though) or\n";
		message1 += " add -DUSE_PTHREADS to the CPP_FLAGS in your Makefile and recompile\n";
		throw PsimagLite::RuntimeError(message1.c_str());
	}

	void usageMain(const PsimagLite::String& name) const
	{
		std::cerr<<"USAGE is "<<name<<"\n";
	}

private:

	bool error1(const PsimagLite::String& message,SizeType line) const
	{
		PsimagLite::String s(__FILE__);
		s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
		throw PsimagLite::RuntimeError(s.c_str());

	}

	OptionsReadableType* optsReadable_;

}; // class InputCheck
} // namespace FreeFermions

/*@}*/
#endif
