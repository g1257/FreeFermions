// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file FermionFactor.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef FERMION_FACTOR_H
#define FERMION_FACTOR_H

#include "Complex.h" // in PsimagLite
#include "FreeOperators.h" // in PsimagLite


namespace FreeFermions {


template<typename OperatorType,typename OpPointerType>
class FermionFactor {
public:
	typedef typename OperatorType::RealType RealType;
	//typedef typename OperatorType::FieldType FieldType;
	typedef FreeOperators<OperatorType,OpPointerType> FreeOperatorsType;

	enum {CREATION = OperatorType::CREATION,
		  DESTRUCTION = OperatorType::DESTRUCTION
	     };

	FermionFactor(FreeOperatorsType& freeOps)
	    : value_(1)
	{
		if (freeOps()==0) {
			value_ = 0;
			return;
		}
		freeOps.removeNonCsOrDs();
		freeOps.reverse();

		pairUp(freeOps);
	}

	RealType operator()()
	{
		return value_;
	}

private:

	void pairUp(//const typename PsimagLite::Vector<SizeType>::Type& occupations,
	            FreeOperatorsType& freeOps)
	{
		// pair up daggers with non-daggers
		value_ = 1;
		while(freeOps.size()>0) {
			if (freeOps.notCreationOrDestruction(freeOps[0].type))
				continue;
			// take first operator's lambda:
			SizeType thisLambda = freeOps[0].lambda;
			// find next operator with same lambda:
			int x = freeOps.findOpGivenLambda(thisLambda,1);

			// if types are equal then result is zero, and we're done:
			if (x<0 || freeOps[0].type == freeOps[x].type) {
				value_ = 0;
				return;
			}
			// then we produce a sign, and either an occupation
			// or an anti-occupation
			// let's deal with the (anti)occupation first:

			bool nNormal = (freeOps[0].type == CREATION) ? true : false;
			SizeType occupationForThisLamda = 0; //occupations[thisLambda];
			if (nNormal && occupationForThisLamda==0) {
				value_=0;
				return;
			}

			if (!nNormal && occupationForThisLamda==1) {
				value_=0;
				return;
			}

			// now let's deail with the sign
			SizeType xSaved = x;
			x++;
			int f = (x&1) ? -1 : 1;
			value_ *= f;

			// finally, we remove the pair
			freeOps.removePair(xSaved);
		}
	}

	RealType value_;
}; // FermionFactor


} // namespace Dmrg 

/*@}*/
#endif // FERMION_FACTOR_H
