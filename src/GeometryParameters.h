/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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
#include "String.h"
#include "IoSimple.h"

namespace FreeFermions {

	template<typename RealType_,typename IoInputType>
	struct GeometryParameters {

		typedef RealType_ RealType;

		enum {CHAIN,LADDER,FEAS,KTWONIFFOUR,FEAS1D};

		enum {DIRECTION_X=0,DIRECTION_Y=1,DIRECTION_XPY=2,DIRECTION_XMY=3};

		GeometryParameters(IoInputType& io)
		: type(CHAIN),
		  sites(0),
		  leg(0),
		  isPeriodic(2,false),
		  hopping(1,1.0),
		  orbitals(1)
		{
			io.readline(sites,"TotalNumberOfSites=");

			// default value
			type = CHAIN;

			PsimagLite::String geometry("");
			io.readline(geometry,"GeometryKind=");

			PsimagLite::String model("");
			io.readline(model,"Model=");

			int x = 0;

			io.readline(x,"IsPeriodicX=");
			if (x<0)
				throw std::runtime_error("GeometryParameters: IsPeriodicX must be non negative\n");
			isPeriodic[DIRECTION_X] = (x>0);

			if (model=="HubbardOneBand" || model=="HubbardOneOrbital")  {//    Immm  Tj1Orb

				if (geometry=="chain") return;
				if (geometry!="ladder")
					throw std::runtime_error("GeometryParameters: HubbardOneBand supports geometry chain or ladder only\n");

				type = LADDER;
				hopping.resize(2);
				hopping[0] =  hopping[1]  = 1.0;

				io.readline(x,"IsPeriodicY=");
				if (x<0)
					throw std::runtime_error("GeometryParameters: IsPeriodicY must be non negative\n");
				isPeriodic[DIRECTION_Y] = (x>0);

				io.readline(x,"LadderLeg=");
				if (x<0)
					throw std::runtime_error("GeometryParameters: HubbardOneBand ladder leg must be non negative\n");
				leg = x;
				return;
			} else if (model=="FeAsBasedSc") {

				io.readline(x,"IsPeriodicY=");
				if (x<0)
					throw std::runtime_error("GeometryParameters: IsPeriodicY must be non negative\n");
				isPeriodic[DIRECTION_Y] = (x>0);

				try {
					io.readline(x,"LadderLeg=");
				} catch (std::exception& e) {
					x=1;
				}
				if (x!=1 && x!=2)
					throw std::runtime_error("GeometryParameters: HubbardOneBand ladder leg must be 1 or 2 only\n");
				leg = x;
				type = (leg==1) ? FEAS1D :  FEAS;

				if (x<1 || x>3)
					throw std::runtime_error("GeometryParameters: HubbardOneBand ladder orbitals must be 1, 2 or 3 only\n");
				io.readline(x,"Orbitals=");
				orbitals = x;

			} else if (model=="Immm") {
				type = KTWONIFFOUR;
				io.readline(x,"IsPeriodicY=");
				if (x<0)
					throw std::runtime_error("GeometryParameters: IsPeriodicY must be non negative\n");
				isPeriodic[DIRECTION_Y] = (x>0);
			} else {
				PsimagLite::String str("GeometryParameters: unknown model ");
				str += model + "\n";
				throw std::runtime_error(str.c_str());
			}
		}

		static int readElectrons(IoInputType& io,size_t nsites)
		{
			int x = -1;
			try {
				io.readline(x,"TargetElectronsUp=");
			} catch (std::exception& e)
			{
				x=-1;
			}

			typename PsimagLite::Vector<RealType>::Type v;
			try {
				io.read(v,"TargetQuantumNumbers");
			} catch (std::exception& e)
			{
				v.resize(0);
			}

			if (x<0 && v.size()==0) {
				PsimagLite::String str("Either TargetElectronsUp or TargetQuantumNumbers is need\n");
				throw std::runtime_error(str.c_str());
			}

			if (x>=0 && v.size()>0) {
				PsimagLite::String str("Having both TargetElectronsUp and TargetQuantumNumbers is an error\n");
				throw std::runtime_error(str.c_str());
			}

			if (x>=0) return x;
			if (v.size()<1) {
				PsimagLite::String str("Incorrect TargetQuantumNumbers line\n");
				throw std::runtime_error(str.c_str());
			}
			return static_cast<size_t>(round(v[0]*nsites));
		}

		template<typename VectorLikeType>
		static void readVector(VectorLikeType& v,const PsimagLite::String& filename,const PsimagLite::String& label)
		{
			PsimagLite::IoSimple::In io(filename);
			io.read(v,label);
		}

		template<typename SomeFieldType>
		static void readLabel(SomeFieldType& v,const PsimagLite::String& filename,const PsimagLite::String& label)
		{
			PsimagLite::IoSimple::In io(filename);
			io.readline(v,label);
		}

		size_t type;
		size_t sites;
		size_t leg;
		typename PsimagLite::Vector<bool>::Type isPeriodic;
		typename PsimagLite::Vector<RealType>::Type hopping;
		PsimagLite::String filename;
		size_t orbitals;
	}; // struct GeometryParameters
} // namespace Dmrg 

/*@}*/
#endif //GEOMETRY_PARAMS_H
