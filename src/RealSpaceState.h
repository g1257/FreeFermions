// BEGIN LICENSE BLOCK
/*
Copyright (c) 2011 , UT-Battelle, LLC
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file RealSpaceState.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef REAL_SPACE_STATE_H
#define REAL_SPACE_STATE_H

#include <assert.h>
#include "Complex.h" // in PsimagLite
#include "TypeToString.h"
#include "FlavoredState.h"
#include "Permutations.h"
#include "ArrangementsWithoutRepetition.h"
#include "Sort.h"

namespace FreeFermions {

	template<typename CorDOperatorType_>
	class RealSpaceState {
		typedef typename CorDOperatorType_::EngineType EngineType;
		typedef typename CorDOperatorType_::RealType RealType;
		typedef typename CorDOperatorType_::FieldType FieldType;
		typedef  FlavoredState<CorDOperatorType_> FlavoredStateType;
		typedef RealSpaceState<CorDOperatorType_> ThisType;
		typedef ArrangementsWithoutRepetition<typename PsimagLite::Vector<SizeType>::Type >
						ArrangementsWithoutRepetitionType;
		typedef PsimagLite::Permutations<ArrangementsWithoutRepetitionType> PermutationsType;

		enum {CREATION = CorDOperatorType_::CREATION,
		       DESTRUCTION = CorDOperatorType_::DESTRUCTION,
		       DIAGONAL
		};

	public:
		typedef CorDOperatorType_ CorDOperatorType;

		// it's the g.s. for now, FIXME change it later to allow more flex.
		RealSpaceState(const EngineType& engine,
		               const typename PsimagLite::Vector<SizeType>::Type& ne,
		               SizeType threadNum,
		               bool debug)
		:  engine_(&engine),ne_(ne),debug_(debug),sorted_(false),zeroVals_(0)
		{
			for (SizeType i=0;i<engine_->dof();i++)
				initTerms(i,threadNum);
		}

		void pushInto(const CorDOperatorType& op)
		{
			for (SizeType i=0;i<terms_.size();i++) {
				if (fabs(values_[i])<1e-8) continue;
				int x = terms_[i].apply(op.type(),op.sigma(),op.index());
				values_[i] *= x;
				if (x==0) zeroVals_++;
			}
		}

		FieldType scalarProduct(ThisType& other)
		{
			simplify();
			other.simplify();
			FieldType sum = 0;
			SizeType j=0;
//			std::cout<<"size1="<<terms_.size()<<" size2=";
//			std::cout<<other.terms_.size()<<"\n";
			for (SizeType i=0;i<other.terms_.size();i++) {
				while (j<terms_.size() && terms_[j]<other.terms_[i]) j++;
				SizeType k = j;
				while(k<terms_.size() && terms_[k]==other.terms_[i]) {
					sum += PsimagLite::conj(values_[k]) * other.values_[i];
					k++;
				}
			}
			return sum;
		}

	private:

		void simplify()
		{
			if (terms_.size()==0) return;

			killZeroVals();
			sort();

			FlavoredStateType prev = terms_[0];
			typename PsimagLite::Vector<FlavoredStateType>::Type terms2 = terms_;
			terms_.clear();
			typename PsimagLite::Vector<FieldType>::Type values2 = values_;
			values_.clear();
			terms_.push_back(prev);
			values_.push_back(values_[0]);
			for (SizeType i=1;i<terms2.size();i++) {
				if (terms2[i] == prev) continue;
				terms_.push_back(terms2[i]);
				values_.push_back(values2[i]);
			}
		}

		void sort()
		{
			if (terms_.size()<2 || sorted_) return;
			typename PsimagLite::Vector<SizeType>::Type iperm(terms_.size());
			PsimagLite::Sort<typename PsimagLite::Vector<FlavoredStateType>::Type > mysort;
			mysort.sort(terms_,iperm);
			typename PsimagLite::Vector<FieldType>::Type valuesNew(values_.size());
			for (SizeType i=0;i<values_.size();i++)
				valuesNew[i]=values_[iperm[i]];
			values_=valuesNew;
			sorted_ = true;
		}

		// \sum_p (-1)^p \sum_{lambda}
		// U_{lambda(0),p(0)} U_{lambda(1),p(1)} U_{lambda(2),p(2)}
		// ... c^\dagger_{p(0)} c^\dagger_{p(1)} c^\dagger_{p(2)}...
		// where the sum is over all permutations p of 0,1,2 ... N-1
		void initTerms(SizeType sigma, SizeType threadNum)
		{
			assert(engine_->dof()==1);
			SizeType n = engine_->size();
			typename PsimagLite::Vector<bool>::Type v(n,false);
			if (ne_[sigma]==0) {
				FlavoredStateType fl(engine_->dof(),v.size(),threadNum);
				fl.pushInto(sigma,v);
				terms_.push_back(fl);
				values_.push_back(1.0);
				return;
			}

			ArrangementsWithoutRepetitionType ap(n,ne_[sigma]);
			//std::cerr<<"ap.size="<<ap.size()<<"\n";
			while (ap.increase()) {
				PermutationsType p(ap);
				//std::cerr<<"p.size="<<p.size()<<" terms="<<terms_.size()<<"\n";
				for (SizeType i=0;i<v.size();i++) v[i] = false;
				for (SizeType i=0;i<ne_[sigma];i++) v[p[i]] = true;
				RealType sum = 0;
				do {
//					std::cerr<<"--------> "<<p<<" <---------\n";
					RealType prod = (isArrangementOdd(p)) ? -1.0 : 1.0;
					for (SizeType i=0;i<ne_[sigma];i++) {
						prod *= engine_->eigenvector(p[i],i);
					}
					sum += prod;
				} while (p.increase());
				FlavoredStateType fl(engine_->dof(),v.size(),threadNum);
				fl.pushInto(sigma,v);
				terms_.push_back(fl);
				values_.push_back(sum);
			};
		}

		template<typename SomeVectorType>
		bool isArrangementOdd(const SomeVectorType& v)
		{
			// make a copy because sort will modify it:
			typedef typename SomeVectorType::value_type SomeElementType;
			typename PsimagLite::Vector<SomeElementType>::Type w(v.size());
			for (SizeType i=0;i<v.size();i++) w[i] = v[i];
			PsimagLite::Sort<typename PsimagLite::Vector<SomeElementType>::Type > mysort;
			typename PsimagLite::Vector<SizeType>::Type iperm(w.size());
			mysort.sort(w,iperm);
			return isOdd(iperm);
		}

		bool isOdd(const typename PsimagLite::Vector<SizeType>::Type& x)
		{
			//Return even parity for the permutation
			int temp =  (x.size() - ncycles(x));
			return temp & 1;
		}

		void killZeroVals()
		{
			FlavoredStateType prev = terms_[0];
			typename PsimagLite::Vector<FlavoredStateType>::Type  terms2 = terms_;
			terms_.clear();
			typename PsimagLite::Vector<FieldType>::Type values2 = values_;
			values_.clear();
			terms_.push_back(prev);
			values_.push_back(values_[0]);
			for (SizeType i=1;i<terms2.size();i++) {
				if (fabs(values2[i])<1e-8) continue;
				terms_.push_back(terms2[i]);
				values_.push_back(values2[i]);
			}
			zeroVals_=0;
		}

//		void killZeroVals()
//		{
//			for (SizeType i=1;i<terms_.size();i++) {
//				if (fabs(values_[i])<1e-8) {
//					terms_.erase(terms_.begin()+i);
//					values_.erase(values_.begin()+i);
//					i--;
//				}
//			}
//			zeroVals_=0;
//		}

		SizeType ncycles(const typename PsimagLite::Vector<SizeType>::Type& x)
		{
			SizeType ncycles = 0;
			typename PsimagLite::Vector<bool>::Type seen(x.size(),false);

			for (SizeType i=0;i<seen.size();i++) {
				if (seen[i]) continue;
				ncycles++;
				//mark indices that belong to the cycle
				SizeType j = i;
				while (!seen[j]) {
					seen[j] = true;
					j = x[j];
				}
			}
			return ncycles;
		}
		const EngineType* engine_;
		typename PsimagLite::Vector<SizeType>::Type ne_;
		bool debug_;
		bool sorted_;
		SizeType zeroVals_;
		typename PsimagLite::Vector<FlavoredStateType>::Type terms_;
		typename PsimagLite::Vector<FieldType>::Type values_;
	}; // RealSpaceState
	
	template<typename CorDOperatorType>
	typename CorDOperatorType::FieldType scalarProduct(
	      const RealSpaceState<CorDOperatorType>& s1,
	      const RealSpaceState<CorDOperatorType>& s2)
	{
		RealSpaceState<CorDOperatorType> s3 = s2;
		RealSpaceState<CorDOperatorType> s4 = s1;
		return s3.scalarProduct(s4);
	}

} // namespace Dmrg 

/*@}*/
#endif //REAL_SPACE_STATE_H
