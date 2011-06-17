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

/*! \file HilbertVector.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef HILBERT_VECTOR_H
#define HILBERT_VECTOR_H

#include "Complex.h" // in PsimagLite
#include "Sort.h" // in PsimagLite

namespace FreeFermions {
	
	template<typename RealType_,typename HilbertTermType_>
	class HilbertVector {
			
			//static size_t const TERMS_MAX_SOFT = 500000;
			//static size_t const TERMS_MAX_HARD = 5000000;
			typedef HilbertVector<RealType_,HilbertTermType_> ThisType;
			
		public:
			typedef RealType_ RealType;
			typedef HilbertTermType_ HilbertTermType;
			typedef typename HilbertTermType::FieldType FieldType;
			typedef typename HilbertTermType::StateType StateType;
			typedef typename HilbertTermType::FlavorFactoryType
			                                    FlavorFactoryType;
			
			HilbertVector(size_t size,size_t dof,bool verbose=false) :
				size_(size),dof_(dof),verbose_(verbose),sorted_(false)
			{
			}
			
			void add(const ThisType& another)
			{
				if (verbose_) std::cerr<<"Adding "<<another.terms()<<" to "<<terms_.size()<<"\n";
				for (size_t i=0;i<another.terms();i++)
					add(another.term(i));
				
			}
			
			// No grouping here, since it's too expensive
			// Use simplify if you need to group
			void add(const HilbertTermType& term)
			{
				terms_.push_back(term);
				sorted_ = false;
			}
			
			const HilbertTermType& term(size_t i) const
			{
				return terms_[i];
			}
			
			size_t terms() const { return terms_.size(); }
			
			void fill(const std::vector<size_t>& ne)
			{
				if (ne.size()!=dof_)
					throw std::runtime_error("HilbertVector::fill()\n");
				
				StateType fstate(dof_);
				FlavorFactoryType flavorFactory(fstate,size_);
				flavorFactory.fill(ne);
				clear();
				HilbertTermType term(fstate,1.0);
				terms_.push_back(term);
				sorted_ = false;
			}
			
			void clear()
			{
				terms_.clear();
				sorted_ = false;
			}
			
			FieldType scalarProduct(ThisType& v)
			{
				if (size_!=v.size_ || dof_!=v.dof_)
					throw std::runtime_error("ScalarProduct\n");
				FieldType sum = 0;
				sort();
				v.sort();
				size_t j=0;
				for (size_t i=0;i<v.terms_.size();i++) {
					
					while (j<terms_.size() && terms_[j].state<v.terms_[i].state) j++;
					size_t k = j;
					while(k<terms_.size() && terms_[k].state==v.terms_[i].state) {
						sum += std::conj(terms_[k].value) * v.terms_[i].value;
						k++;
					}
				}
				
				return sum;
			}

			// sort and then simplify
			void simplify()
			{
				if (terms_.size()==0) return;
				std::vector<size_t> iperm(terms_.size());
				if (verbose_) std::cerr<<"Tring to sort "<<terms_.size()<<" items.\n";
				Sort<std::vector<HilbertTermType> > sortObject;
				std::vector<HilbertTermType> termsOld = terms_;
				sortObject.sort(terms_,iperm);

				std::vector<HilbertTermType> termsNew;

				size_t prevIndex=0;
				StateType prev = terms_[0].state;

				if (verbose_) std::cerr<<"Will now simplify...\n";
				for (size_t i=0;i<terms_.size();i++) {
					if (i==0 || !(terms_[i].state==prev)) {
						termsNew.push_back(terms_[i]);
						prevIndex = termsNew.size()-1;
						prev = terms_[i].state;
					} else {
						termsNew[prevIndex].value += termsOld[iperm[i]].value;
					}
				}
				termsOld.clear();
				terms_ = termsNew;
				sorted_ = true;
				if (verbose_) std::cerr<<"Simplification done "<<terms_.size()<<"\n";
			}
			
			// sort 
			void sort()
			{
				if (terms_.size()<2 || sorted_) return;
				std::vector<size_t> iperm(terms_.size());
				if (verbose_) std::cerr<<"Sorting "<<terms_.size()<<" items.\n";
				Sort<std::vector<HilbertTermType> > sortObject;
				sortObject.sort(terms_,iperm);
				sorted_ = true;
			}

			
			template<typename T,typename V>
			friend std::ostream& operator<<(std::ostream& os,const HilbertVector<T,V>& v);

		private:	
			size_t size_;
			size_t dof_;
			bool verbose_;
			bool sorted_;
			std::vector<HilbertTermType> terms_;
			
	}; // HilbertVector
	
	template<typename T,typename V>
	std::ostream& operator<<(std::ostream& os,const HilbertVector<T,V>& v)
	{
		
		os<<"size="<<v.size_<<"\n";
		os<<"dof="<<v.dof_<<"\n";
		os<<"dataSize="<<v.terms_.size()<<"\n";
		for (size_t i=0;i<v.terms_.size();i++) {
			os<<v.terms_[i];
		}
		return os;
	}
	
	template<typename T,typename V>
	typename V::FieldType scalarProduct(
	    HilbertVector<T,V>& v1,
	    HilbertVector<T,V>& v2)
	{
		return v1.scalarProduct(v2);
	}
} // namespace Dmrg 

/*@}*/
#endif
