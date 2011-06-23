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

/*! \file FreeOperators.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef FREE_OPERATORS_H
#define FREE_OPERATORS_H

#include "Complex.h" // in PsimagLite
#include "Sort.h" // in PsimagLite
#include "Permutations.h"
#include "IndexGenerator.h"

namespace FreeFermions {
	
	struct FreeOperator {
		size_t lambda;
		size_t type;
	};

	template<typename OperatorType,typename OpPointerType>
	class FreeOperators {
	public:
		typedef typename OperatorType::RealType RealType;
		typedef typename OperatorType::FieldType FieldType;
		typedef IndexGenerator IndexGeneratorType;
		typedef Permutations<IndexGeneratorType> PermutationsType;
		typedef std::vector<OpPointerType> OpPointersType;

		enum {CREATION = OperatorType::CREATION,
		       DESTRUCTION = OperatorType::DESTRUCTION
		};

		FreeOperators(const OpPointersType& opPointers,
		                const IndexGeneratorType& lambda,
		                const PermutationsType& lambda2)
		{
			size_t counter=0;
			size_t counter2=0;

			for (size_t i=0;i<opPointers.size();i++) {
				FreeOperator fo;
				fo.type = opPointers[i].type;
				if (fo.type==CREATION) {
					fo.lambda = lambda[counter++];
				} else {
					fo.lambda = lambda2[counter2++];
				}
				data_.push_back(fo);
			}
			// if daggers > non-daggers, result is zero
			if (counter!=counter2) {
				value_ = 0;
				return;
			}

		}

		size_t size() const { return data_.size(); }

		const FreeOperator& operator[](size_t i) const { return data_[i]; }

		void reverse()
		{
			// flip'em
			std::vector<FreeOperator> dataCopy = data_;
			size_t n = data_.size();
			for (size_t i=0;i<n;i++)
				data_[i] = dataCopy[n-i-1];
		}

		RealType operator()()
		{
			return value_;
		}

		void removePair(size_t thisLambda)
		{
			std::vector<FreeOperator>::iterator itp = data_.begin();
			data_.erase(itp);
			// find again because the erase changed data:
			int y = findOpGivenLambda(thisLambda,0);
			itp = data_.begin()+y;
			data_.erase(itp);
		}

		int findOpGivenLambda(size_t thisLambda,
		                         size_t start) const
		{
			for (size_t i=start;i<data_.size();i++) {
					if (data_[i].lambda==thisLambda) return i;
			}
			throw std::runtime_error("FreeOperators::findOpGivenLambda()\n");
		}

	private:
		std::vector<FreeOperator> data_;
		RealType value_;

//		void doSort(std::vector<size_t>& v)
//		{
//			Sort<std::vector<size_t> > mysort;
//			std::vector<size_t> iperm(v.size());
//			mysort.sort(v,iperm);
//			//return (repeatedIndices(v)) ? 0 : parityOfPerm(iperm);
//		}

		// note that v must come in ordered
//		bool repeatedIndices(std::vector<size_t>& v) const
//		{
//			size_t prev = v[0];
//			for (size_t i=1;i<v.size();i++) {
//				if (v[i]==prev) return true;
//				prev = v[i];
//			}
//			return false;
//		}

//
//		int parityOfPerm(std::vector<size_t>& iperm)
//		{
//			std::vector<bool> elementsSeen(iperm.size(),false);
//			size_t cycles = 0;
//			for (size_t i=0;i<elementsSeen.size();i++) {
//				if (elementsSeen[i]) continue;
//				cycles++;
//				size_t current = i;
//				while (!elementsSeen[current]) {
//					elementsSeen[current] = true;
//					current = iperm[current];
//				}
//			}
//			int x = iperm.size() - cycles;
//			return (x &1) ? -1 : 1;
//		}

		//const OpPointersType& opPointers_;
        //const IndexGeneratorType& lambda_;
	}; // FreeOperators
	

} // namespace Dmrg 

/*@}*/
#endif // FREE_OPERATORS_H
