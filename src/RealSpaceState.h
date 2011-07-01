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

/*! \file RealSpaceState.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef REAL_SPACE_STATE_H
#define REAL_SPACE_STATE_H

#include "Complex.h" // in PsimagLite
#include "TypeToString.h"
#include "FlavoredState.h"
#include "ArrangementsWithoutRepetition.h"
#include "Sort.h"

namespace FreeFermions {

	template<typename CorDOperatorType_>
	class RealSpaceState {
		typedef typename CorDOperatorType_::EngineType EngineType;
		typedef typename CorDOperatorType_::RealType RealType;
		typedef typename CorDOperatorType_::FieldType FieldType;
		typedef  FlavoredState<std::vector<bool>,CorDOperatorType_> FlavoredStateType;
		typedef RealSpaceState<CorDOperatorType_> ThisType;

		enum {CREATION = CorDOperatorType_::CREATION,
		       DESTRUCTION = CorDOperatorType_::DESTRUCTION,
		       DIAGONAL
		};

	public:
		typedef CorDOperatorType_ CorDOperatorType;

		// it's the g.s. for now, FIXME change it later to allow more flex.
		RealSpaceState(const EngineType& engine,
		              const std::vector<size_t>& ne,
		              bool debug = false)
		:  engine_(&engine),ne_(ne),debug_(debug),sorted_(false)
		{
			for (size_t i=0;i<engine_->dof();i++)
				initTerms(i);
		}

		void pushInto(const CorDOperatorType& op)
		{
			for (size_t i=0;i<terms_.size();i++) {
				int x = terms_[i].apply(op.type(),op.sigma(),op.index());
				values_[i] *= x;
			}
			FIXME CREATE NEW TERM IF NEEDED
			// remove 0 values FIXME
		}

		FieldType scalarProduct(ThisType& other)
		{
			sort();
			other.sort();
			FieldType sum = 0;
			size_t j=0;
			for (size_t i=0;i<other.terms_.size();i++) {
				while (j<terms_.size() && terms_[j]<other.terms_[i]) j++;
				size_t k = j;
				while(k<terms_.size() && terms_[k]==other.terms_[i]) {
					sum += std::conj(values_[k]) * other.values_[i];
					k++;
				}
			}
			return sum;
		}

	private:

		void sort()
		{
			if (terms_.size()<2 || sorted_) return;
			std::vector<size_t> iperm(terms_.size());
			Sort<std::vector<FlavoredStateType> > mysort;
			mysort.sort(terms_,iperm);
			std::vector<FieldType> valuesNew(values_.size());
			for (size_t i=0;i<values_.size();i++)
				valuesNew[i]=values_[iperm[i]];
			values_=valuesNew;
			sorted_ = true;
		}

		// \sum_p (-1)^p \sum_{lambda}
		// U_{lambda(0),p(0)} U_{lambda(1),p(1)} U_{lambda(2),p(2)}
		// ... c^\dagger_{p(0)} c^\dagger_{p(1)} c^\dagger_{p(2)}...
		// where the sum is over all permutations p of 0,1,2 ... N-1
		void initTerms(size_t sigma)
		{
			size_t n = engine_->sites();
			std::vector<bool> v(n);
			typedef ArrangementsWithoutRepetition<std::vector<size_t> >
				ArrangementsWithoutRepetitionType;
			ArrangementsWithoutRepetitionType p(n,ne_[sigma]);
			do {
				RealType prod = (isArrangementEven(p)) ? 1.0 : -1.0;
				for (size_t i=0;i<ne_[sigma];i++) {
					prod *= engine_->eigenvector(p[i],i);
				}

				for (size_t i=0;i<v.size();i++) v[i] = false;
				for (size_t i=0;i<p.size();i++) v[p[i]] = true;
				assert(engine_->dof()==1);
				FlavoredStateType fl(engine_->dof(),v.size());
				fl.pushInto(sigma,v);
				terms_.push_back(fl);
				values_.push_back(prod);
			} while (p.increase());
		}

		template<typename SomeVectorType>
		bool isArrangementEven(const SomeVectorType& v)
		{
			// make a copy because sort will modify it:
			typedef typename SomeVectorType::value_type SomeElementType;
			std::vector<SomeElementType> w(v.size());
			for (size_t i=0;i<v.size();i++) w[i] = v[i];
			Sort<std::vector<SomeElementType> > mysort;
			std::vector<size_t> iperm(w.size());
			mysort.sort(w,iperm);
			return even(iperm);
		}

		static bool even(const std::vector<size_t>& x)
		{
			//Return even parity for the permutation
			return (x.size() - ncycles(x)) % 2 == 0;
		}

		static size_t ncycles(const std::vector<size_t>& x)
		{
			size_t ncycles = 0;
			std::vector<bool> seen(x.size(),false);

			for (size_t i=0;i<seen.size();i++) {
				if (seen[i]) continue;
				ncycles++;
				//mark indices that belongs to the cycle
				size_t j = i;
				while (!seen[j]) {
					seen[j] = true;
					j = x[j];
				}
			}
			return ncycles;
		}
		const EngineType* engine_;
		std::vector<size_t> ne_;
		bool debug_;
		bool sorted_;
		std::vector<FlavoredStateType> terms_;
		std::vector<FieldType> values_;
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
