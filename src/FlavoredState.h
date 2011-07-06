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

/*! \file FlavoredState.h
 *
 * 
 *
 */
#ifndef FLAVORED_STATE_H
#define FLAVORED_STATE_H

namespace FreeFermions {
	// All interactions == 0
	template<typename LevelsType,typename OperatorType>
	class FlavoredState {
			//static size_t const SPIN_UP=0,SPIN_DOWN=1;
			typedef FlavoredState<LevelsType,OperatorType> ThisType;
			static int const FERMION_SIGN = -1;
			enum {CREATION = OperatorType::CREATION,
			       DESTRUCTION = OperatorType::DESTRUCTION};

		public:
			FlavoredState(size_t dof,size_t size)
			{

				for (size_t i=0;i<dof;i++) {
					LevelsType tmpV(size,false);
					data_.push_back(tmpV);
				}
			}

			void pushInto(size_t sigma,const LevelsType& portion)
			{
				data_[sigma] = portion;
			}
			
//			void fill(const std::vector<size_t>& ne)
//			{
//				if (ne.size()!=data_.size()) throw std::runtime_error(
//						"FlavoredState::fill()\n");
//				for (size_t i=0;i<ne.size();i++) { // sum over spins
//					fillInternal(data_[i],ne[i]);
//				}
//			}
			
			int apply(size_t label,size_t flavor,size_t lambda)
			{
				if (flavor>=data_.size())
					throw std::runtime_error("FlavoredState::create()\n");
				if (lambda>=data_[0].size())
					throw std::runtime_error("FlavoredState::create()\n");
				int interSign =
						(calcInterElectrons(flavor) %2) ? 1 : FERMION_SIGN;
				return applyInternal(label,data_[flavor],lambda)*interSign;
			}
			
			
//			void occupations(std::vector<size_t>& ns,size_t flavor) const
//			{
//				ns.resize(size_);
//				for (size_t i = 0; i < size_; i++) ns[i] = 0;
//
//				for (size_t counter=0;counter<data_[flavor].size();counter++) {
//					ns[counter] = (data_[flavor][counter]) ? 1 : 0;
//				}
//			}
			
			size_t flavors() const { return data_.size(); }
			
			template<typename T,typename U>
			friend std::ostream& operator<<(std::ostream& os,
			                                  const FlavoredState<T,U>& v);
			
			template<typename T,typename U>
			friend bool operator==(const FlavoredState<T,U>& v1,
			                         const FlavoredState<T,U>& v2);
			
			template<typename T,typename U>
			friend bool operator<(const FlavoredState<T,U>& v1,
			                         const FlavoredState<T,U>& v2);
					
			template<typename T,typename U>
			friend bool operator>(const FlavoredState<T,U>& v1,
			                         const FlavoredState<T,U>& v2);
					
			template<typename T,typename U>
			friend bool operator<=(const FlavoredState<T,U>& v1,
			                         const FlavoredState<T,U>& v2);
			
		private:
			
			void fillInternal(LevelsType& x,size_t ne)
			{
				if (ne>x.size())
					throw std::runtime_error("FlavoredState::fillInternal\n");
				for (size_t i=0;i<x.size();i++) x[i] = (i<ne) ? true : false;
			}
			
			int applyInternal(size_t label,LevelsType& x,size_t lambda)
			{
				size_t nflips = statesBetween(x,lambda);
				if (label == CREATION) {
					if (x[lambda]) return 0; // can't create, there's already one
					x[lambda] = true;
				} else if (label == DESTRUCTION) {
					if (!x[lambda]) return 0; // can't destroy, there's nothing
					x[lambda] =false;
				} else {
					throw std::runtime_error("FlavoredState::applyInternal()\n");
				}
				if (nflips ==0 || nflips % 2 ==0) return 1;
				return FERMION_SIGN;
			}
			
			size_t statesBetween(LevelsType x,size_t lambda) const
			{
				size_t sum = 0;
				for (size_t counter = 0;counter < lambda ; counter++)
					if (x[counter]) sum ++;
				return sum;
			}
			
			size_t calcInterElectrons(size_t flavor)
			{
				size_t sum = 0;
				for (size_t flavor2 = 0; flavor2 < flavor; flavor2++) {
					sum += numberOfDigits(data_[flavor2]);
				}
				return sum;
			}
			
			size_t numberOfDigits(const LevelsType& x)
			{
				size_t sum = 0;
				for (size_t i=0;i<x.size();i++) {
					if (x[i]) sum ++;
				}
				return sum;
			}

			std::vector<LevelsType> data_;
	}; // FlavoredState
	
	template<typename T,typename U>
	std::ostream& operator<<(std::ostream& os,const FlavoredState<T,U>& v)
	{
		os<<"size="<<v.data_.size()<<"\n";
		for (size_t i=0;i<v.data_.size();i++) {
			os<<v.data_[i]<<" ";
		}
		return os;
	}
	
	template<typename T,typename U>
	inline bool operator==(const FlavoredState<T,U>& v1,const FlavoredState<T,U>& v2)
	{
		// eliminated due to performance reasons:
		//if (size_!=b.size_ || data_.size()!=b.data_.size()) return false;

		for (size_t i=0;i<v1.data_.size();i++)
			if (v1.data_[i]!=v2.data_[i]) return false;
		return true;
	}

	inline bool operator<(const std::vector<bool>& v1,const std::vector<bool>& v2)
	{
		for (size_t i=0;i<v1.size();i++) {
			if (v1[i] && !v2[i]) return false;
			if (!v1[i] && v2[i]) return true;
		}
		return false;
	}

	template<typename T,typename U>
	inline bool operator<(const FlavoredState<T,U>& v1,const FlavoredState<T,U>& v2)
	{
		// eliminated due to performance reasons:
		//if (size_!=b.size_ || data_.size()!=b.data_.size()) return false;
		
		for (size_t i=0;i<v1.data_.size();i++) {
			if (v1.data_[i]>v2.data_[i]) return false;
			if (v1.data_[i]<v2.data_[i]) return true;
		}
		return false;
	}
	
//	template<typename T,typename U>
//	inline bool operator>(const FlavoredState<T,U>& v1,const FlavoredState<T,U>& v2)
//	{
//		// eliminated due to performance reasons:
//		//if (size_!=b.size_ || data_.size()!=b.data_.size()) return false;
//
//		for (size_t i=0;i<v1.data_.size();i++) {
//			if (!v1.data_[i] && v2.data_[i]) return false;
//			if (v1.data_[i]&& !v2.data_[i]) return true;
//		}
//		return false;
//	}
//
//	template<typename T,typename U>
//	inline bool operator<=(const FlavoredState<T,U>& v1,const FlavoredState<T,U>& v2)
//	{
//		// eliminated due to performance reasons:
//		//if (size_!=b.size_ || data_.size()!=b.data_.size()) return false;
//
//		for (size_t i=0;i<v1.data_.size();i++) {
//			if (v1.data_[i] && !v2.data_[i]) return false;
//			if (!v1.data_[i] && v2.data_[i]) return true;
//		}
//		return true;
//	}
} // namespace Dmrg 

/*@}*/
#endif
