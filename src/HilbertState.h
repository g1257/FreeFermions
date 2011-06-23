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
#include "FermionFactor.h"

namespace FreeFermions {
	struct OperatorPointer {
		OperatorPointer(size_t t,size_t i) : type(t),index(i) {}
		size_t type;
		size_t index;
	};

	template<typename FieldType>
	class DummyOperator {
	public:
		DummyOperator(const DummyOperator* x) {}
		template<typename T1,typename T2>
		FieldType operator()(const T1& l1,const T2& l2) const { return 1; }
		size_t type() const { return 0; }
		void transpose() {}
	};

	template<typename CorDOperatorType,
	          typename DiagonalOperatorType=
	                    DummyOperator<typename CorDOperatorType::RealType> >
	class HilbertState {
		typedef typename CorDOperatorType::RealType RealType;
		typedef typename CorDOperatorType::FieldType FieldType;
		typedef FermionFactor<CorDOperatorType,OperatorPointer> FermionFactorType;
		typedef typename FermionFactorType::PermutationsType PermutationsType;
		typedef typename FermionFactorType::IndexGeneratorType IndexGeneratorType;

		typedef HilbertState<CorDOperatorType,DiagonalOperatorType> ThisType;

		enum {CREATION = CorDOperatorType::CREATION,
		       DESTRUCTION = CorDOperatorType::DESTRUCTION,
		       DIAGONAL
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
			for (size_t i=0;i<garbage2_.size();i++)
				delete garbage2_[i];
		}

		void pushInto(const CorDOperatorType& op)
		{
			if (op.type()==CREATION) {
				OperatorPointer opPointer(op.type(),operatorsCreation_.size());
				operatorsCreation_.push_back(&op);
				opPointers_.push_back(opPointer);
			} else if (op.type()==DESTRUCTION) {
				OperatorPointer opPointer(op.type(),operatorsDestruction_.size());
				operatorsDestruction_.push_back(&op);
				opPointers_.push_back(opPointer);
			}
		}

		void pushInto(const DiagonalOperatorType& op)
		{
				OperatorPointer opPointer(op.type(),operatorsDiagonal_.size());
				operatorsDiagonal_.push_back(&op);
				opPointers_.push_back(opPointer);
		}

		void pour(const ThisType& hs)
		{
			if (hs.hilbertSize_!=hilbertSize_ || hs.ne_!=ne_)
				throw std::runtime_error("HilbertState::pour(...)\n");

			size_t counter = 0;
			size_t counter2 = 0;
			size_t counter3 = 0;
			size_t n1 = hs.opPointers_.size();
			// pour them in reverse order:
			for (size_t i=0;i<hs.opPointers_.size();i++) {
				if (hs.opPointers_[n1-i-1].type==CREATION) {
					size_t x1 = hs.operatorsCreation_.size() - 1 - counter;
					const CorDOperatorType* op = hs.operatorsCreation_[x1];
					counter++;
					CorDOperatorType *opCopy = new CorDOperatorType(op);
					garbage_.push_back(opCopy);
					opCopy->transpose();
					pushInto(*opCopy);
				} else if (hs.opPointers_[n1-i-1].type==DESTRUCTION) {
					size_t x2 = hs.operatorsDestruction_.size() - 1 - counter2;
					const CorDOperatorType* op = hs.operatorsDestruction_[x2];
					counter2++;
					CorDOperatorType *opCopy = new CorDOperatorType(op);
					garbage_.push_back(opCopy);
					opCopy->transpose();
					pushInto(*opCopy);
				} else {
					size_t x3 = hs.operatorsDiagonal_.size() - 1 - counter3;
					const DiagonalOperatorType* op = hs.operatorsDiagonal_[x3];
					counter3++;
					DiagonalOperatorType *opCopy = new DiagonalOperatorType(op);
					garbage2_.push_back(opCopy);
					opCopy->transpose();
					pushInto(*opCopy);
				}
			}

		}

		FieldType close()
		{
			size_t m = operatorsCreation_.size();
			IndexGeneratorType lambda(m,hilbertSize_);
			FieldType sum  = 0;
			do {
				sum += compute(lambda);
			} while (lambda.increase());
			return sum;
		}

	private:



		FieldType compute(const IndexGeneratorType& lambda)
		{

			PermutationsType lambda2(lambda);
			FieldType sum = 0;
			do  {
				FieldType prod = 1;
				FermionFactorType fermionFactor(opPointers_,lambda,lambda2,ne_);
				FieldType ff = fermionFactor();
				if (debug_) {
					std::cerr<<"sum ="<<sum<<" prod="<<prod<<" ff="<<fermionFactor();
					std::cerr<<" l="<<lambda<<" l2="<<lambda2<<"\n";
				}
				if (ff==0) continue;
				for (size_t i=0;i<lambda.size();i++) {
					prod *= operatorsCreation_[i]->operator()(lambda[i]);
				}
				for (size_t i=0;i<lambda2.size();i++) {
					prod *= operatorsDestruction_[i]->operator()(lambda2[i]);
				}
				FieldType dd = 1.0;
				for (size_t i=0;i<operatorsDiagonal_.size();i++)
					dd *= operatorsDiagonal_[i]->operator()(lambda,lambda2);

				sum += prod*ff*dd;

			} while(lambda2.increase());
			return sum;
		}

		size_t hilbertSize_,ne_;
		bool debug_;
		std::vector<const CorDOperatorType*> operatorsCreation_,operatorsDestruction_;
		std::vector<const DiagonalOperatorType*> operatorsDiagonal_;
		std::vector<OperatorPointer> opPointers_;
		std::vector<const CorDOperatorType*> garbage_;
		std::vector<const DiagonalOperatorType*> garbage2_;
	}; // HilbertState
	
	template<typename CorDOperatorType,typename DiagonalOperatorType>
	typename CorDOperatorType::FieldType scalarProduct(
	      const HilbertState<CorDOperatorType,DiagonalOperatorType>& s1,
	      const HilbertState<CorDOperatorType,DiagonalOperatorType>& s2)
	{
		HilbertState<CorDOperatorType,DiagonalOperatorType> s3 = s2;
		s3.pour(s1); // s1 --> s3
		return s3.close();
	}

} // namespace Dmrg 

/*@}*/
#endif //HILBERT_STATE_H
