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

namespace FreeFermions {
	
	template<typename EngineType,typename OperatorType>
	class HilbertState {
		typedef typename EngineType::RealType RealType;
		typedef typename EngineType::FieldType FieldType;
		typedef Permutations PermutationsType;

	public:
		// it's the g.s. for now, FIXME change it later to allow more flex.
		HilbertState(const EngineType& engine,size_t ne)
		: engine_(engine), ne_(ne) {}

		void applyOperator(const OperatorType& op)
		{
			operators_.push_back(op);
		}

		FieldType close()
		{
			size_t n = operators_.size();
			if (n&1) throw std::runtime_error(
				"HilbertState::close(): n. of operators is odd\n");
			size_t m = n/2;
			//orderOperators(); // so that operators appear
			                  // dagger and non dagger
			PermutationsType perm(m,m);
			FieldType sum = 0;
			do  {
				sum += compute(perm);
			} while (perm.increase());
			return sum;
		}

	private:

		FieldType compute(const PermutationsType& perm)
		{
			PermutationsType lambda(perm.size(),ne_);
			FieldType sum = 0;
			do  {
				FieldType prod = 1;
				for (size_t i=0;i<lambda.size();i++) {
					prod *=
					  std::conj(engine_.eigenvector(operators_[i].site(),lambda[i]));
				}
				for (size_t i=0;i<lambda.size();i++) {
					prod *=
					  engine_.eigenvector(operators_[perm[i]].site(),lambda[perm[i]]);
				}
				sum += prod;
			} while (lambda.increase());
			return sum;
		}

		const EngineType& engine_;
		size_t ne_;
		std::vector<OperatorType> operators_;

	}; // HilbertState
	
	template<typename EngineType,typename OperatorType>
	typename EngineType::FieldType scalarProduct(HilbertState<EngineType,OperatorType>& s1,
	                           HilbertState<EngineType,OperatorType>& s2)
	{
		// s2.pour(s1); // s1 --> s2
		return s2.close();
	}

} // namespace Dmrg 

/*@}*/
#endif //HILBERT_STATE_H
