/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[FreeFermions, Version 1.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
/** \ingroup DMRG */
/*@{*/

/*! \file GeometryParameters.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef GEOMETRY_PARAMS_H
#define GEOMETRY_PARAMS_H

#include <assert.h>
#include "Io/IoSimple.h"
#include "AnsiColors.h"

namespace FreeFermions {

template<typename FieldType,typename IoInputType>
struct GeometryParameters {

	typedef typename PsimagLite::Real<FieldType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	enum {CHAIN,
		  LADDER,
		  LADDERX,
		  FEAS,
		  KTWONIFFOUR,
		  FEAS1D,
		  CHAIN_EX,
		  LADDER_BATH,
		  STAR,
		  KANE_MELE_HUBBARD,
		  RAW};

	enum {DIRECTION_X=0,DIRECTION_Y=1,DIRECTION_XPY=2,DIRECTION_XMY=3};

	GeometryParameters(IoInputType& io,FieldType defaultHopping = 1.0)
	    : type(CHAIN),
	      sites(0),
	      leg(0),
	      isPeriodic(2,false),
	      hopping(1,defaultHopping),
	      filename(io.filename()),
	      orbitals(1),
	      omega(0.0),
	      phase(0.0)
	{
		outputFile = "bogus.txt";
		if (outputFile.find(".") == PsimagLite::String::npos)
			outputFile += ".txt";

		io.readline(sites,"TotalNumberOfSites=");

		// default value
		type = CHAIN;

		io.readline(geometry,"GeometryKind=");

		try {
			io.readline(omega,"omega=");
		} catch (std::exception&) {}

		try {
			io.readline(phase,"phase=");
		} catch (std::exception&) {}

		if (geometry == "LongRange" || geometry == "Raw") {
			type = RAW;
			return;
		}

		if (geometry == "chainEx") type = CHAIN_EX;
		if (geometry == "star") type = STAR;

		PsimagLite::String model("");
		io.readline(model,"Model=");

		PsimagLite::String hoppingOptions;
		io.readline(hoppingOptions,"GeometryOptions=");

		bool constantHoppings = (hoppingOptions == "ConstantValues") ? true : false;

		int x = 0;

		if (type != STAR) io.readline(x,"IsPeriodicX=");
		if (x<0) {
			PsimagLite::String str("GeometryParameters: IsPeriodicX must be 0 or 1\n");
			throw std::runtime_error(str);
		}

		isPeriodic[DIRECTION_X] = (x>0);

		try {
			VectorRealType tmpVector;
			io.read(tmpVector, "hubbardU");
			bool flag = false;
			for (SizeType i = 0; i <tmpVector.size(); ++i) {
				if (tmpVector[i] == 0) continue;
				flag = true;
				break;
			}

			if (flag) {
				std::cerr<<PsimagLite::AnsiColor::red;
				std::cerr<<"WARNING: HubbardU is non-zero!\n";
				std::cerr<<PsimagLite::AnsiColor::reset;
			}
		} catch (std::exception&) {}

		if (model=="HubbardOneBand" ||
		        model=="HubbardOneOrbital" ||
		        model == "SuperHubbardExtended")  {

			if (geometry == "chain" || geometry == "chainEx") return;
			bool ladderOrLadderX = (geometry == "ladder" || geometry == "ladderx");
			if (!ladderOrLadderX && geometry != "ladderbath") {
				PsimagLite::String str("GeometryParameters: This model supports");
				throw std::runtime_error(str + " chain or ladder or ladderbath only\n");
			}

			if (geometry == "ladderbath")  {
				io.read(bathSitesPerSite,"BathSitesPerSite=");
				if (sites % (bathSitesPerSite+1) != 0)
					throw std::runtime_error("GeometryParameters:\n");
				sites = static_cast<SizeType>(sites/(bathSitesPerSite+1));
				type = LADDER_BATH;
			} else {
				type = LADDER;
			}

			io.readline(x,"LadderLeg=");
			if (x < 0) {
				PsimagLite::String str("GeometryParameters: HubbardOneBand ladder leg");
				err(str + " must be non negative\n");
			}

			leg = x;

			hopping.resize(2*sites - leg, defaultHopping);
			if (geometry == "ladderx") {
				hopping.resize(4*sites - 3*leg, defaultHopping);
				type = LADDERX;
			}

			if (ladderOrLadderX && constantHoppings) {
				typename PsimagLite::Vector<FieldType>::Type v;
				io.read(v,"Connectors");
				if (v.size() != 1)
					err("Connectors must be vector of size 1\n");
				for (SizeType i = 0; i < sites - leg; ++i)
					hopping[i] = v[0];
				io.read(v,"Connectors");
				if (v.size() != 1)
					err("Connectors must be vector of size 1\n");
				for (SizeType i = 0; i < sites; ++i)
					hopping[i + sites - leg] = v[0];

				if (geometry == "ladderx") {
					typename PsimagLite::Vector<FieldType>::Type v;
					io.read(v,"Connectors");
					if (v.size() != 1)
						err("Connectors must be vector of size 1\n");
					for (SizeType i = 0; i < sites - leg; ++i)
						hopping[i + 2*sites - leg] = v[0];
					io.read(v,"Connectors");
					if (v.size() != 1)
						err("Connectors must be vector of size 1\n");
					for (SizeType i = 0; i < sites - leg; ++i)
						hopping[i + 3*sites - 2*leg] = v[0];
				}
			}

			if (!constantHoppings) {
				typename PsimagLite::Vector<FieldType>::Type v;
				io.read(v,"Connectors");
				assert(v.size() == sites - leg);
				for (SizeType i = 0; i < v.size(); ++i)
					hopping[i] = v[i];
				io.read(v,"Connectors");
				assert(v.size() == sites);
				for (SizeType i = 0; i < v.size(); ++i)
					hopping[i + sites - leg] = v[i];
				if (geometry == "ladderx") {
					io.read(v,"Connectors");
					assert(v.size() == sites - leg);
					for (SizeType i = 0; i < v.size(); ++i)
						hopping[i + 2*sites - leg] = v[i];
					io.read(v,"Connectors");
					assert(v.size() == sites - leg);
					for (SizeType i = 0; i < v.size(); ++i)
						hopping[i + 3*sites - 2*leg] = v[i];
				}
			}

			if (geometry == "ladderbath")  {
				io.read(tb,"Connectors");
				if (tb.size() != bathSitesPerSite * sites)
					err("ladderbath Connectors\n");
			}

			io.readline(x,"IsPeriodicY=");
			if (x<0) {
				PsimagLite::String str("GeometryParameters: ");
				err(str + "IsPeriodicY must be 0 or 1\n");
			}

			isPeriodic[DIRECTION_Y] = (x>0);


			return;
		} else if (model=="FeAsBasedSc" || model == "FeAsBasedScExtended") {

			io.readline(x,"Orbitals=");
			orbitals = x;

			if (geometry == "chain" || geometry == "chainEx") {
				leg = 1;
				type = FEAS1D;
				return;
			}

			x = 0;
			if (geometry == "ladder" || geometry == "ladderx")
				io.readline(x,"IsPeriodicY=");
			if (x<0) {
				PsimagLite::String str("GeometryParameters: ");
				throw std::runtime_error(str + "IsPeriodicY must be 0 or 1\n");
			}

			isPeriodic[DIRECTION_Y] = (x>0);

			try {
				io.readline(x,"LadderLeg=");
			} catch (std::exception&) {
				if (geometry == "ladder" || geometry == "ladderx")
					throw std::runtime_error("No LadderLeg=\n");
			}

			if (x!=1 && x!=2) {
				PsimagLite::String str("GeometryParameters: HubbardOneBand ladder leg");
				throw std::runtime_error(str + " must be 1 or 2 only\n");
			}

			leg = x;
			if (type != STAR)
				type = (leg==1) ? FEAS1D :  FEAS;

			if (x<1 || x>3) {
				PsimagLite::String str("GeometryParameters: HubbardOneBand ladder leg");
				throw std::runtime_error(str + " must be 1 or 2 or 3 only\n");
			}

		} else if (model=="Immm") {
			type = KTWONIFFOUR;
			io.readline(x,"IsPeriodicY=");
			if (x<0) {
				PsimagLite::String str("GeometryParameters: ");
				throw std::runtime_error(str + "IsPeriodicY must be 0 or 1\n");
			}

			isPeriodic[DIRECTION_Y] = (x>0);
		} else if (model=="KaneMeleHubbard") {
			type = KANE_MELE_HUBBARD;
			io.readline(x,"NumberOfTerms=");
			if (x != 2) {
				PsimagLite::String str("GeometryParameters: KaneMeleHubbard expects");
				throw std::runtime_error(str + " 2 terms\n");
			}

			if (geometry != "ladderx" || constantHoppings) {
				PsimagLite::String str("GeometryParameters: KaneMeleHubbard expects");
				throw std::runtime_error(str + " ladderx with non-constant\n");
			}

			if (constantHoppings) {
				PsimagLite::String str("GeometryParameters: KaneMeleHubbard expects");
				throw std::runtime_error(str + " non constant hoppings");
			}

			io.readline(geometry,"GeometryKind=");
			io.readline(hoppingOptions,"GeometryOptions=");

			if (geometry != "longchain" || hoppingOptions != "ConstantValues") {
				PsimagLite::String str("GeometryParameters: KaneMeleHubbard expects");
				throw std::runtime_error(str + " longchain with ConstantValues\n");
			}

		} else {
			PsimagLite::String str("GeometryParameters: unknown model ");
			str += model + "\n";
			throw std::runtime_error(str.c_str());
		}
	}

	static int readElectrons(IoInputType& io,SizeType nsites)
	{
		int x = -1;

		try {
			io.readline(x,"TargetElectronsUp=");
		} catch (std::exception& e) {
			io.readline(x,"TargetSzPlusConst=");
		}

		if (x >= 0) return x;

		PsimagLite::String str("Incorrect nup/sz line\n");
		throw std::runtime_error(str.c_str());
	}

	template<typename VectorLikeType>
	static void readVector(VectorLikeType& v,
	                       const PsimagLite::String& filename,
	                       const PsimagLite::String& label)
	{
		PsimagLite::IoSimple::In io(filename);
		io.read(v,label);
	}

	template<typename SomeFieldType>
	static void readLabel(SomeFieldType& v,
	                      const PsimagLite::String& filename,
	                      const PsimagLite::String& label)
	{
		PsimagLite::IoSimple::In io(filename);
		io.readline(v,label);
	}

	PsimagLite::String geometry;
	SizeType type;
	SizeType sites;
	SizeType leg;
	typename PsimagLite::Vector<bool>::Type isPeriodic;
	typename PsimagLite::Vector<FieldType>::Type hopping;
	PsimagLite::String outputFile;
	PsimagLite::String filename;
	SizeType orbitals;
	SizeType bathSitesPerSite;
	RealType omega;
	RealType phase;
	typename PsimagLite::Vector<FieldType>::Type tb;
}; // struct GeometryParameters
} // namespace Dmrg

/*@}*/
#endif //GEOMETRY_PARAMS_H

