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

#include "FlavoredState.h"
#include "Complex.h" // in PsimagLite
#include "Sort.h" // in PsimagLite

namespace FreeFermions {
	// All interactions == 0
	
	template<typename FieldType,typename FlavoredStateType>
	struct HilbertTerm {
		HilbertTerm(const FlavoredStateType& state1,const FieldType& value1) :
				state(state1),value(value1) { }
		FlavoredStateType state;
		FieldType value;
	};
	
	template<typename T,typename V>
	std::ostream& operator<<(std::ostream& os,const HilbertTerm<T,V>& v)
	{
		os<<v.state<<"\n";
		os<<"value="<<v.value<<"\n";
		return os;	
	}
	
	template<typename RealType,typename HilbertTermType>
	RealType energy(const std::vector<RealType>& eigs,HilbertTermType& term)
	{
		return term.value*energy(eigs,term.state);
	}
	
	template<typename RealType_,typename FieldType_,typename LevelsType>
	class HilbertVector {
			
			static size_t const TERMS_MAX_SOFT = 500000;
			static size_t const TERMS_MAX_HARD = 5000000;
			typedef HilbertVector<RealType_,FieldType_,LevelsType> ThisType;
			
		public:
			typedef RealType_ RealType;
			typedef FieldType_ FieldType;
			typedef FlavoredState<LevelsType> FlavoredStateType;
			typedef HilbertTerm<FieldType,FlavoredStateType> HilbertTermType;
			
			HilbertVector(size_t size,size_t dof,bool verbose=false) :
				size_(size),dof_(dof),verbose_(verbose),sorted_(false)
			{
			}
			
			void add(const ThisType& another)
			{
				if (verbose_) std::cerr<<"Adding "<<another.terms()<<" to "<<data_.size()<<"\n";
				for (size_t i=0;i<another.terms();i++)
					add(another.term(i));
				
			}
			
			// No grouping here, since it's too expensive
			// Use simplify if you need to group
			void add(const HilbertTermType& term)
			{
				const FlavoredStateType& state = term.state;
				const FieldType& value = term.value;
				data_.push_back(state);
				values_.push_back(value);
				sorted_ = false;
			}
			
			HilbertTermType term(size_t i) const
			{
				return HilbertTermType(data_[i],values_[i]);
			}
			
			size_t terms() const { return data_.size(); }
			
			void fill(const std::vector<size_t>& ne)
			{
				if (ne.size()!=dof_) throw std::runtime_error("HilbertVector::fill()\n");
				
				FlavoredStateType fstate(dof_,size_);
				fstate.fill(ne);
				clear();
				data_.push_back(fstate);
				values_.push_back(1.0);
				sorted_ = false;
			}
			
			void clear()
			{
				data_.clear();
				values_.clear();
				sorted_ = false;
			}
			
			// This function needs to be robust enough to handle the case
			// where neither this nor v are grouped
// 			FieldType scalarProduct1(const ThisType& v) const
// 			{
// 				if (size_!=v.size_ || dof_!=v.dof_) throw std::runtime_error("ScalarProduct\n");
// 				FieldType sum = 0;
// 				typedef typename std::vector<FlavoredStateType>::const_iterator MyIterator;
// 
// 				for (size_t i=0;i<v.data_.size();i++) {
// 					MyIterator start = data_.begin();
// 					while(start!=data_.end()) {
// 						MyIterator x = find(start,data_.end(),v.data_[i]);
// 						if (x==data_.end()) break;
// 						sum += std::conj(values_[x-data_.begin()]) * v.values_[i];
// 						start = x+1;
// 					}
// 				}
// 				return sum;
// 			}
			
			FieldType scalarProduct(ThisType& v)
			{
				if (size_!=v.size_ || dof_!=v.dof_) throw std::runtime_error("ScalarProduct\n");
				FieldType sum = 0;
				sort();
				v.sort();
				size_t j=0;
				for (size_t i=0;i<v.data_.size();i++) {
					
					while (j<data_.size() && data_[j]<v.data_[i]) j++;
					size_t k = j;
					while(k<data_.size() && data_[k]==v.data_[i]) {
						sum += std::conj(values_[k]) * v.values_[i];
						k++;
					}
				}
				
				return sum;
			}
			
			// this is an expensive operation due to the search:
// 			void simplifySlow()
// 			{
// 				std::vector<FlavoredStateType> dataNew;
// 				std::vector<FieldType> valuesNew;
// 				for (size_t i=0;i<data_.size();i++) {
// 					int x = utils::isInVector(dataNew,data_[i]);
// 					if (x<0) {
// 						dataNew.push_back(data_[i]);
// 						valuesNew.push_back(values_[i]);
// 					} else {
// 						valuesNew[x] += values_[i];
// 					}
// 				}
// 				data_ = dataNew;
// 				values_ = valuesNew;
// 				std::cerr<<"Simplification done\n";
// 			}
// 			
// 			// sort and then simplify
// 			void simplifyLowMemory()
// 			{
// 				if (data_.size()==0) return;
// 				std::vector<size_t> iperm(data_.size());
// 				std::cerr<<"Tring to sort "<<data_.size()<<" items.\n";
// 				utils::sort(data_,iperm);
// 				//throw std::runtime_error("Don't forget to reorder values\n");
// 				std::vector<FieldType> valuesNew(values_.size(),0);
// 				//for (size_t i=0;i<values_.size();i++) valuesNew[i]=values_[iperm[i]];
// 				//values_=valuesNew;
// 				//return;
// 				size_t prevIndex=0;
// 				FlavoredStateType prev = data_[0];
// 				//for (size_t i=0;i<values_.size();i++) valuesNew[i]=0;
// 				std::cerr<<"Will now simplify...\n";
// 				for (size_t i=0;i<data_.size();i++) {
// 					if (i==0 || !(data_[i]==prev)) {
// 						valuesNew[i] = values_[iperm[i]];
// 						prevIndex = i;
// 						prev = data_[i];
// 					} else {
// 						valuesNew[prevIndex] += values_[iperm[i]];
// 						valuesNew[i] = static_cast<FieldType>(0.0);
// 					}
// 				}
// 				sorted_ = true;
// 				std::cerr<<"Now erasing...\n";
// 				iperm.clear();
// 				values_=valuesNew;
// 				valuesNew.clear();
// 				typedef typename std::vector<FlavoredStateType>::iterator MyIterator;
// 				MyIterator iter = data_.begin();
// 				typedef typename std::vector<FieldType>::iterator MyIterator2;
// 				MyIterator2 iter2 = values_.begin();
// 				while(iter!=data_.end() && iter2!=values_.end()) {
// 					if (*iter2 == static_cast<FieldType>(0.0)) {
// 						iter = data_.erase(iter);
// 						iter2 = values_.erase(iter2);
// 					} else {
// 						iter++,iter2++;
// 					}
// 				}
// 				std::cerr<<"Simplification done "<<data_.size()<<"\n";
// 			}
// 			
			// sort and then simplify
			void simplify()
			{
				if (data_.size()==0) return;
				std::vector<size_t> iperm(data_.size());
				if (verbose_) std::cerr<<"Tring to sort "<<data_.size()<<" items.\n";
				Sort<std::vector<FlavoredStateType> > sortObject;
				sortObject.sort(data_,iperm);
				//throw std::runtime_error("Don't forget to reorder values\n");
				std::vector<FieldType> valuesNew;
				std::vector<FlavoredStateType> dataNew;
				//for (size_t i=0;i<values_.size();i++) valuesNew[i]=values_[iperm[i]];
				//values_=valuesNew;
				//return;
				size_t prevIndex=0;
				FlavoredStateType prev = data_[0];
				//for (size_t i=0;i<values_.size();i++) valuesNew[i]=0;
				if (verbose_) std::cerr<<"Will now simplify...\n";
				for (size_t i=0;i<data_.size();i++) {
					if (i==0 || !(data_[i]==prev)) {
						valuesNew.push_back(values_[iperm[i]]);
						dataNew.push_back(data_[i]);
						prevIndex = valuesNew.size()-1;
						prev = data_[i];
					} else {
						valuesNew[prevIndex] += values_[iperm[i]];
					}
				}
				values_ = valuesNew;
				valuesNew.clear();
				data_ = dataNew;
				sorted_ = true;
				if (verbose_) std::cerr<<"Simplification done "<<data_.size()<<"\n";
			}
			
			// sort 
			void sort()
			{
				if (data_.size()<2 || sorted_) return;
				std::vector<size_t> iperm(data_.size());
				if (verbose_) std::cerr<<"Sorting "<<data_.size()<<" items.\n";
				Sort<std::vector<FlavoredStateType> > sortObject;
				sortObject.sort(data_,iperm);
				std::vector<FieldType> valuesNew(values_.size());
				for (size_t i=0;i<values_.size();i++) valuesNew[i]=values_[iperm[i]];
				/*std::cerr<<values_;
				std::cerr<<"---------------\n";
				std::cerr<<valuesNew;
				std::cerr<<"-+-+-+-+-+-+-+-\n";
				std::cerr<<data_<<"\n";
				std::cerr<<"-+-+-+-+-+-+-+-\n";*/
				values_=valuesNew;
				sorted_ = true;
			}
// 			
// 			ThisType& operator*=(const FieldType &rhs)
// 			{
// 				throw std::runtime_error("Operator *= unimplemented in HilbertVerctor.h\n");
// 			}

			
			template<typename T,typename V, typename U>
			friend std::ostream& operator<<(std::ostream& os,const HilbertVector<T,V,U>& v);

		private:	
			size_t size_;
			size_t dof_;
			bool verbose_;
			bool sorted_;
			std::vector<FlavoredStateType> data_;
			std::vector<FieldType> values_;
			
	}; // HilbertVector
	
	template<typename T,typename V, typename U>
	std::ostream& operator<<(std::ostream& os,const HilbertVector<T,V,U>& v)
	{
		
		os<<"size="<<v.size_<<"\n";
		os<<"dof="<<v.dof_<<"\n";
		os<<"dataSize="<<v.data_.size()<<"\n";
		for (size_t i=0;i<v.data_.size();i++) {
			typename HilbertVector<T,V,U>::HilbertTermType term(v.data_[i],v.values_[i]);
			os<<term;
		}
		return os;
	}
	
	template<typename T,typename V, typename U>
	V scalarProduct(HilbertVector<T,V,U>& v1,HilbertVector<T,V,U>& v2)
	{
		return v1.scalarProduct(v2);
	}
} // namespace Dmrg 

/*@}*/
#endif
