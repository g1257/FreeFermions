/*
Copyright (c) 2009-2017, UT-Battelle, LLC
All rights reserved

[FreeFermions, Version 1.]
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
/** \ingroup DMRG */
/*@{*/

/*! \file Combinations.h
 *
 *  Algorithm by Donald Knuth.
 *
 */
#ifndef COMBINATIONS_H_H
#define COMBINATIONS_H_H

#include <vector>

namespace FreeFermions {
class Combinations {
public:

	Combinations(SizeType n,SizeType k)
	    : k_(k)
	{
		if (k==0 || n<k) throw std::runtime_error(
		            "Combinations::ctor\n");
		mainLoop(n);
	}

	const PsimagLite::Vector<SizeType>::Type& operator()(SizeType i) const
	{
		return data_[i];
	}

	SizeType size() const { return data_.size(); }

private:
	void mainLoop(SizeType n)
	{
		PsimagLite::Vector<int>::Type c(k_+3);
		PsimagLite::Vector<SizeType>::Type tmpVec(k_);
		int x = 0;
		for (int i=1; i <= int(k_); i++) c[i] = i;
		c[k_+1] = n+1;
		c[k_+2] = 0;
		int j = k_;
visit:
		SizeType counter = 0;
		for (int i=k_; i >= 1; i--)
			tmpVec[counter++] = c[i]-1;
		data_.push_back(tmpVec);

		if (j > 0) {x = j+1; goto incr;}
		if (c[1] + 1 < c[2]) {
			c[1] += 1;
			goto visit;
		}
		j = 2;
do_more:
		c[j-1] = j-1;
		x = c[j] + 1;
		if (x == c[j+1]) {j++; goto do_more;}
		if (j > int(k_)) return;
incr:
		c[j] = x;
		j--;
		goto visit;
	}

	SizeType k_;
	PsimagLite::Vector<PsimagLite::Vector<SizeType>::Type>::Type data_;
}; // Combinations

std::ostream& operator<<(std::ostream& os,const Combinations& ig)
{
	for (SizeType i=0;i<ig.size();++i) {
		for (SizeType j=0;j<ig(i).size();++j)
			os<<ig(i)[j]<<" ";
		os<<"\n";
	}
	return os;
}
} // namespace Dmrg 

/*@}*/
#endif // COMBINATIONS_H_H
