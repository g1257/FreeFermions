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
	typedef FreeOperators<OperatorType,OpPointerType> FreeOperatorsType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef std::pair<SizeType, SizeType> PairSizeType;

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

		freeOps.clear();
	}

	RealType operator()()
	{
		return value_;
	}

private:

	class WhoIsOn {

	public:

		WhoIsOn(SizeType n) : data_(n, true), beginEnd_(0, n) {}

		SizeType findActualIndex(SizeType ind) const
		{
			SizeType n = data_.size();
			PairSizeType bE = fromPair(ind, n);
			for (SizeType i = bE.first; i < bE.second; ++i)
				if (data_[i]) return i;

			throw std::runtime_error("FreeOperators::findActualIndex()\n");
		}

		SizeType count(SizeType from, SizeType to) const
		{
			PairSizeType bE = fromPair(from, to);
			SizeType counter = 0;
			for (SizeType i = bE.first; i < bE.second; ++i)
				if (data_[i]) ++counter;

			return counter;
		}

		PairSizeType fromPair(SizeType from, SizeType to) const
		{
			PairSizeType bE(from, to);
			if (beginEnd_.first > from) bE.first = beginEnd_.first;
			if (beginEnd_.second < to) bE.second = beginEnd_.second;
			return bE;
		}

		bool operator[](SizeType ind) const
		{
			assert(ind < data_.size());
			return data_[ind];
		}

		void set(SizeType ind, bool value)
		{
			assert(ind < data_.size());
			data_[ind] = value;
		}

		void updateBounds()
		{
			for (SizeType i = beginEnd_.first; i < beginEnd_.second; ++i) {
				if (!data_[i]) continue;
				if (i > beginEnd_.first) beginEnd_.first = i;
				break;
			}

			for (SizeType i = 0; i < beginEnd_.second; ++i) {
				SizeType j = beginEnd_.second - i - 1;
				if (!data_[j]) continue;
				if (j + 1 < beginEnd_.second)  beginEnd_.second = j + 1;
				break;
			}
		}

	private:

		VectorBoolType data_;
		PairSizeType beginEnd_;
	};

	void pairUp(FreeOperatorsType& freeOps)
	{
		// pair up daggers with non-daggers
		value_ = 1;
		WhoIsOn bits(freeOps.size());
		SizeType count = freeOps.size();

		while(count > 0) {
			SizeType firstIndex = bits.findActualIndex(0);

			if (FreeOperatorsType::notCreationOrDestruction(freeOps[firstIndex].type))
				continue;
			// take first operator's lambda:
			SizeType thisLambda = freeOps[firstIndex].lambda;
			// find next operator with same lambda:
			SizeType firstIndexPlusOne = bits.findActualIndex(firstIndex + 1);
			int x = findOpGivenLambda(thisLambda, firstIndexPlusOne, bits, freeOps);

			// if types are equal then result is zero, and we're done:
			if (x < 0 || freeOps[firstIndex].type == freeOps[x].type) {
				value_ = 0;
				return;
			}

			// then we produce a sign, and either an occupation
			// or an anti-occupation
			// let's deal with the (anti)occupation first:

			bool nNormal = (freeOps[firstIndex].type == CREATION) ? true : false;
			SizeType occupationForThisLamda = 0; //occupations[thisLambda];
			if (nNormal && occupationForThisLamda == 0) {
				value_ = 0;
				return;
			}

			if (!nNormal && occupationForThisLamda == 1) {
				value_ = 0;
				return;
			}

			// now let's deal with the sign
			const SizeType xSaved = x;
			SizeType counter = bits.count(firstIndex, xSaved + 1);

			int f = (counter & 1) ? -1 : 1;
			value_ *= f;

			// we remove the pair
			assert(x > 0);

			SizeType actualFirst = bits.findActualIndex(0);
			bits.set(actualFirst, false);

			SizeType second = bits.findActualIndex(x);
			bits.set(second, false);

			bits.updateBounds();
			count -= 2;
		}
	}

	static int findOpGivenLambda(SizeType thisLambda,
	                      SizeType start,
	                      const WhoIsOn& bits,
	                      const FreeOperatorsType& freeOps)
	{
		const SizeType n = freeOps.size();
		PairSizeType bE = bits.fromPair(start, n);
		for (SizeType  i = bE.first; i < bE.second; ++i) {
			if (!bits[i]) continue;
			if (FreeOperatorsType::notCreationOrDestruction(freeOps[i].type)) continue;
			if (freeOps[i].lambda == thisLambda) return i;
		}

		return -1;
	}

	RealType value_;
}; // FermionFactor


} // namespace Dmrg 

/*@}*/
#endif // FERMION_FACTOR_H
