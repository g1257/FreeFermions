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

/*! \file HilbertState.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef HILBERT_STATE_H
#define HILBERT_STATE_H

#include "Complex.h" // in PsimagLite
#include "Permutations.h"
#include "IndexGenerator.h"
#include "FermionFactor.h"

namespace FreeFermions {
	struct OperatorPointer {
		OperatorPointer(size_t t,size_t i) : type(t),index(i) {}
		size_t type;
		size_t index;
	};
	template<typename OperatorType>
	class HilbertState {
		typedef typename OperatorType::RealType RealType;
		typedef typename OperatorType::FieldType FieldType;
		typedef Permutations PermutationsType;
		typedef IndexGenerator IndexGeneratorType;
		typedef FermionFactor<OperatorType,OperatorPointer> FermionFactorType;
		typedef HilbertState<OperatorType> ThisType;

		enum {CREATION = OperatorType::CREATION,
		       DESTRUCTION = OperatorType::DESTRUCTION
		};

	public:
		// it's the g.s. for now, FIXME change it later to allow more flex.
		HilbertState(size_t hilbertSize,size_t ne,bool debug = false)
		:  hilbertSize_(hilbertSize),ne_(ne),debug_(debug) {}

		~HilbertState()
		{
			// get out the garbage:
			for (size_t i=0;i<garbage_.size();i++)
				delete garbage_[i];
		}

		void pushInto(const OperatorType& op)
		{
			if (op.type()==CREATION) {
				OperatorPointer opPointer(op.type(),operatorsCreation_.size());
				operatorsCreation_.push_back(&op);
				opPointers_.push_back(opPointer);
			} else {
				OperatorPointer opPointer(op.type(),operatorsDestruction_.size());
				operatorsDestruction_.push_back(&op);
				opPointers_.push_back(opPointer);
			}
		}

		FieldType close()
		{
			size_t m = operatorsCreation_.size();
			if (operatorsDestruction_.size()!=m) return 0;
			PermutationsType perm(m);
			FieldType sum = 0;
			size_t countFIXME = 0;
			do  {
				sum += compute(perm,countFIXME++);
				if (debug_) std::cerr<<"**********************\n";
			} while (perm.increase());
			return sum;
		}

		void pour(const ThisType& hs)
		{
			if (hs.hilbertSize_!=hilbertSize_ || hs.ne_!=ne_)
				throw std::runtime_error("HilbertState::pour(...)\n");

			size_t counter = 0;
			size_t counter2 = 0;
			size_t test1 = hs.opPointers_.size();
			for (size_t i=0;i<hs.opPointers_.size();i++) {
				if (test1!=hs.opPointers_.size())
					throw std::runtime_error("Changed\n");
				if (hs.opPointers_[i].type==CREATION) {
					const OperatorType* op = hs.operatorsCreation_[counter++];
					OperatorType *opCopy = new OperatorType(op);
					garbage_.push_back(opCopy);
					opCopy->transpose();
					pushInto(*opCopy);
				} else {
					const OperatorType* op = hs.operatorsDestruction_[counter2++];
					OperatorType *opCopy = new OperatorType(op);
					garbage_.push_back(opCopy);
					opCopy->transpose();
					pushInto(*opCopy);
				}
			}

		}

	private:



		FieldType compute(const PermutationsType& perm,size_t countFIXME)
		{
			IndexGeneratorType lambda(perm.size(),hilbertSize_);

			FieldType sum = 0;
			do  {
				FieldType prod = 1;
				FermionFactorType fermionFactor(opPointers_,lambda,perm,ne_,countFIXME);
				FieldType ff = fermionFactor();
				if (ff==0) continue;
				for (size_t i=0;i<lambda.size();i++) {
					prod *= operatorsCreation_[i]->operator()(lambda[i]);
				}
				for (size_t i=0;i<lambda.size();i++) {
					prod *= operatorsDestruction_[i]->operator()(lambda[perm[i]]);
				}

				sum += prod*ff;
				if (debug_) {
					std::cerr<<"sum ="<<sum<<" prod="<<prod<<" ff="<<fermionFactor();
					std::cerr<<" l="<<lambda<<"\n";
				}
			} while (lambda.increase());
			return sum;
		}

		size_t hilbertSize_,ne_;
		bool debug_;
		std::vector<const OperatorType*> operatorsCreation_,operatorsDestruction_;
		std::vector<OperatorPointer> opPointers_;
		std::vector<const OperatorType*> garbage_;
	}; // HilbertState
	
	template<typename OperatorType>
	typename OperatorType::FieldType scalarProduct(
	                                     const HilbertState<OperatorType>& s1,
	                                     const HilbertState<OperatorType>& s2)
	{
		HilbertState<OperatorType> s3 = s2;
		s3.pour(s1); // s1 --> s3
		return s3.close();
	}

} // namespace Dmrg 

/*@}*/
#endif //HILBERT_STATE_H
