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

/*! \file ArrangementsWithoutRepetition.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef ARRANGEMENTS_NOR_H
#define ARRANGEMENTS_NOR_H

#include "Complex.h" // in PsimagLite
#include "Sort.h"

namespace FreeFermions {
	template<typename ContainerType>
	class ArrangementsWithoutRepetition {
		typedef typename ContainerType::value_type FieldType;
	public:
		typedef FieldType value_type;

		ArrangementsWithoutRepetition(size_t n,size_t k)
		: data_(n),k_(k)
		{
			if (k==0 || n<k) throw std::runtime_error(
			  "ArrangementsWithoutRepetition\n");
		}

		/*
		 * Generates the next combination of n elements as k after comb
		 * comb => the previous combination ( use (0, 1, 2, ..., k) for first)
		 * k => the size of the subsets to generate
		 *  n => the size of the original set
		 *  Returns:
		 *  true if a valid combination was found
		 *  false, otherwise
		 */
		bool increase()
		{
			if (data_.size()==0) return false;
			size_t n = data_.size();
			size_t i = k_ - 1;
			++data_[i];
			while (data_[i] >= n - k_ + 1 + i) {
				if (i==0) break;
				--i;
				++data_[i];
			}
			if (data_[0] > n - k_) /* Combination (n-k, n-k+1, ..., n) reached */
				return false; /* No more combinations can be generated */

			/* comb now looks like (..., x, n, n, n, ..., n).
			 * Turn it into (..., x, x + 1, x + 2, ...) */
			for (i = i + 1; i < k_; ++i)
				data_[i] = data_[i - 1] + 1;
			return true;
		}


		size_t operator[](size_t i) const
		{
			return data_[i];
		}

		size_t size() const { return k_; }

	private:

		typename PsimagLite::Vector<size_t>::Type data_;
		size_t k_;
	}; // ArrangementsWithoutRepetition
	
	template<typename T>
	std::ostream& operator<<(std::ostream& os,
	                          const ArrangementsWithoutRepetition<T>& ig)
	{
		for (size_t i=0;i<ig.size();i++) os<<ig[i]<<" ";
		return os;
	}
} // namespace Dmrg 

/*@}*/
#endif // ARRANGEMENTS_NOR_H
