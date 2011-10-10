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
#include "TypeToString.h"

namespace FreeFermions {
	struct OperatorPointer {
		OperatorPointer(size_t t,size_t s,size_t i)
		: type(t),sigma(s),index(i)
		{}

		size_t type;
		size_t sigma;
		size_t index;
	};

	template<typename FieldType>
	class DummyOperator {
	public:
		class DummyFactory {
		public:
			template<typename X>
			DummyFactory(const X& x) {}
			DummyOperator& operator()(const DummyOperator* op) {
				throw std::runtime_error("DummyOperatorFactory\n");
			}
		};
		typedef DummyFactory FactoryType;
		DummyOperator(const DummyOperator* x) {}
		template<typename T1>
		FieldType operator()(const T1& l1,size_t loc) const { return 1; }
		size_t sigma() const { return 0; }
		void transpose() {}
	};

	template<typename CorDOperatorType_,
	          typename DiagonalOperatorType_=
	                    DummyOperator<typename CorDOperatorType_::FieldType> >
	class HilbertState {
		typedef typename CorDOperatorType_::EngineType EngineType;
		typedef typename CorDOperatorType_::RealType RealType;
		typedef typename CorDOperatorType_::FieldType FieldType;
		typedef FermionFactor<CorDOperatorType_,OperatorPointer> FermionFactorType;
		typedef typename FermionFactorType::FreeOperatorsType FreeOperatorsType;
		typedef typename FreeOperatorsType::PermutationsType PermutationsType;
		typedef typename FreeOperatorsType::IndexGeneratorType IndexGeneratorType;

		typedef HilbertState<CorDOperatorType_,DiagonalOperatorType_> ThisType;

		enum {CREATION = CorDOperatorType_::CREATION,
		       DESTRUCTION = CorDOperatorType_::DESTRUCTION,
		       DIAGONAL
		};

	public:
		typedef CorDOperatorType_ CorDOperatorType;
		typedef DiagonalOperatorType_ DiagonalOperatorType;
		typedef typename CorDOperatorType::FactoryType OpNormalFactoryType;
		typedef typename DiagonalOperatorType::FactoryType OpDiagonalFactoryType;

		// it's the g.s. for now, FIXME change it later to allow more flex.
		HilbertState(const EngineType& engine,
		              const std::vector<size_t>& ne,
		              bool debug = false)
		: engine_(&engine),
		  debug_(debug),
		  occupations_(ne.size()),
		  opNormalFactory_(engine),
		  opDiagonalFactory_(engine)
		{
			   for (size_t i=0;i<occupations_.size();++i) {
				   occupations_[i].resize(engine.size(),0);
				   for (size_t j=0;j<ne[i];++j) {
					   occupations_[i][j] = 1;
				   }
			   }
		}

		HilbertState(const EngineType& engine,
		             const std::vector<std::vector<size_t> >& occupations,
		             bool debug = false)
		: engine_(&engine),
		  debug_(debug),
		  occupations_(occupations),
		  opNormalFactory_(engine),
		  opDiagonalFactory_(engine)
		{
		}

		void pushInto(const CorDOperatorType& op)
		{
			if (op.type()==CREATION) {
				OperatorPointer opPointer(op.type(),op.sigma(),
						operatorsCreation_.size());
				operatorsCreation_.push_back(&op);
				opPointers_.push_back(opPointer);
			} else if (op.type()==DESTRUCTION) {
				OperatorPointer opPointer(op.type(),op.sigma(),
						operatorsDestruction_.size());
				operatorsDestruction_.push_back(&op);
				opPointers_.push_back(opPointer);
			}
		}

		void pushInto(const DiagonalOperatorType& op)
		{
				OperatorPointer opPointer(DIAGONAL,0,
						operatorsDiagonal_.size());
				operatorsDiagonal_.push_back(&op);
				opPointers_.push_back(opPointer);
		}

		FieldType pourAndClose(const ThisType& hs)
		{
			pour(hs);
			return close(hs.occupations_);
		}
		
	private:
		void pour(const ThisType& hs)
		{
			if (hs.engine_->size()!=engine_->size()) {
				std::string s = "HilbertState::pour(...)  size1=" +
				                 ttos(engine_->size()) +
				                " size2=" + ttos(hs.engine_->size()) + "\n";
				throw std::runtime_error(s.c_str());
			}
// 			if (hs.occupations_!=occupations_ && !equalZero(occupations_)) {
// 				std::string s(__FILE__);
// 				s += std::string(" ") + ttos(__LINE__) + " ";
// 				s += std::string(__FUNCTION__) + "\n";
// 				throw std::runtime_error(s.c_str());
// 			}

			pourInternal(hs);
		}

		FieldType close(const std::vector<std::vector<size_t> >& occupations2) const
		{
			//std::cerr<<"DEBUG: closing with weight="<<opPointers_.size()<<"\n";
			FieldType prod = 1.0;
			if (occupations_.size()!=occupations2.size())
				throw std::runtime_error("HilbertState::close()\n");

			for (size_t i=0;i<occupations_.size();i++) {
				prod *= close(i,occupations2[i]);
			}
			return prod; // FIXME: NEEDS FERMION SIGN
		}

		bool equalZero(const std::vector<std::vector<size_t> >& v) const
		{
			for (size_t i=0;i<v.size();i++)
				if (!equalZero(v[i])) return false;
			return true;
		}
		
		bool equalZero(const std::vector<size_t>& v) const
		{
			for (size_t i=0;i<v.size();i++) if (v[i]!=0) return false;
			return true;
		}

		FieldType close(size_t sigma,const std::vector<size_t>& occupations2) const
		{
			size_t m = findCreationGivenSpin(sigma);
			IndexGeneratorType lambda(m,engine_->size());
			FieldType sum  = 0;
			do {
				sum += compute(lambda,sigma,occupations2);
			} while (lambda.increase());
			return sum;
		}

		size_t findCreationGivenSpin(size_t sigma) const
		{
			size_t counter = 0;
			for (size_t i=0;i<opPointers_.size();i++) {
				if (opPointers_[i].type == CREATION &&
				    opPointers_[i].sigma == sigma) counter++;
			}
			return counter;
		}

		FieldType compute(const IndexGeneratorType& lambda,
		                  size_t sigma,
		                  const std::vector<size_t>& occupations2) const
		{

			PermutationsType lambda2(lambda);
			FieldType sum = 0;
			do  {
				FieldType prod = 1;
				FreeOperatorsType lambdaOperators(opPointers_,lambda,lambda2,
				                   sigma,occupations_[sigma],occupations2);
				// diag. part need to be done here, because...
				FieldType dd = 1.0;
				for (size_t i=0;i<operatorsDiagonal_.size();i++) {
					size_t loc = lambdaOperators.findLocOfDiagOp(i);
					dd *= operatorsDiagonal_[i]->operator()(lambdaOperators,loc);
				}
				// ... fermionFactor ctor will modify lambdaOperators

				FermionFactorType fermionFactor(lambdaOperators);
				RealType ff = fermionFactor();

				if (ff==0) continue;
				for (size_t i=0;i<lambda.size();i++) {
					int loc = findLocOf(operatorsCreation_,i,sigma);
					if (loc<0) continue;
					prod *= operatorsCreation_[loc]->operator()(lambda[i]);
				}
				for (size_t i=0;i<lambda2.size();i++) {
					int loc = findLocOf(operatorsDestruction_,i,sigma);
					if (loc<0) continue;
					prod *= operatorsDestruction_[loc]->operator()(lambda2[i]);
				}
				if (debug_) {
					std::cerr<<" lambda="<<lambda;
					std::cerr<<" lambda2="<<lambda2;
					std::cerr<<" ff="<<ff<<" dd="<<dd<<" prod="<<prod;
					std::cerr<<"sum ="<<sum<<"\n";
				}
				sum += prod*ff*dd;

			} while(lambda2.increase());
			return sum;
		}

		int findLocOf(
		                  const std::vector<const CorDOperatorType*>& v,
		                  size_t ind,
		                  size_t sigma) const
		{
			size_t counter = 0;
			for (size_t i=0;i<v.size();i++) {
				if (v[i]->sigma() != sigma) continue;
				if (counter==ind) return i;
				counter++;
			}
			return -1;
		}

		void pourInternal(const ThisType& hs)
		{
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
					CorDOperatorType& opCopy = opNormalFactory_(op);
					opCopy.transpose();
					pushInto(opCopy);
				} else if (hs.opPointers_[n1-i-1].type==DESTRUCTION) {
					size_t x2 = hs.operatorsDestruction_.size() - 1 - counter2;
					const CorDOperatorType* op = hs.operatorsDestruction_[x2];
					counter2++;
					CorDOperatorType& opCopy = opNormalFactory_(op);
					opCopy.transpose();
					pushInto(opCopy);
				} else {
					size_t x3 = hs.operatorsDiagonal_.size() - 1 - counter3;
					const DiagonalOperatorType* op = hs.operatorsDiagonal_[x3];
					counter3++;
					DiagonalOperatorType& opCopy = opDiagonalFactory_(op);
					opCopy.transpose();
					pushInto(opCopy);
				}
			}
		}

		const EngineType* engine_;
		bool debug_;
		std::vector<std::vector<size_t> > occupations_;
// 		std::vector<size_t> ne_;
// 		std::vector<size_t> ne2_;
		std::vector<const CorDOperatorType*> operatorsCreation_,operatorsDestruction_;
		std::vector<const DiagonalOperatorType*> operatorsDiagonal_;
		std::vector<OperatorPointer> opPointers_;
		OpNormalFactoryType opNormalFactory_;
		OpDiagonalFactoryType opDiagonalFactory_;
	}; // HilbertState
	
	template<typename CorDOperatorType,typename DiagonalOperatorType>
	typename CorDOperatorType::FieldType scalarProduct(
	      const HilbertState<CorDOperatorType,DiagonalOperatorType>& s1,
	      const HilbertState<CorDOperatorType,DiagonalOperatorType>& s2)
	{
		HilbertState<CorDOperatorType,DiagonalOperatorType> s3 = s2;
		return s3.pourAndClose(s1);
	}

} // namespace Dmrg 

/*@}*/
#endif //HILBERT_STATE_H
