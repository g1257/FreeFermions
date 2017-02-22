/*
Copyright (c) 2011-2017, UT-Battelle, LLC
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

/*! \file FreeOperators.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef FREE_OPERATORS_H
#define FREE_OPERATORS_H

#include "Complex.h" // in PsimagLite
#include "Sort.h" // in PsimagLite
#include "Permutations.h"
#include "IndexGenerator.h"
#include <cassert>
#include <algorithm>

namespace FreeFermions {

struct FreeOperator {
	FreeOperator() : lambda(0),type(0) {}
	SizeType lambda;
	SizeType type;
};

template<typename OperatorType,typename OpPointerType>
class FreeOperators {

	enum {DRY_RUN,NORMAL_RUN};

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename OperatorType::RealType RealType;
	typedef typename OperatorType::FieldType FieldType;
	typedef IndexGenerator IndexGeneratorType;
	typedef PsimagLite::Permutations<IndexGeneratorType> PermutationsType;
	typedef typename PsimagLite::Vector<OpPointerType>::Type OpPointersType;

	enum {CREATION = OperatorType::CREATION,
		  DESTRUCTION = OperatorType::DESTRUCTION,
		  DIAGONAL};

	FreeOperators(const OpPointersType& opPointers,
	              const IndexGeneratorType& lambda,
	              const PermutationsType& lambda2,
	              SizeType sigma,
	              const VectorSizeType& occupations,
	              const VectorSizeType& occupations2)
	    : value_(1),loc_(0)
	{
		SizeType counter3 = dryRun(occupations);
		SizeType counter2 = 0;
		SizeType counter = 0;
		addAtTheMiddle(counter,counter2,opPointers,lambda,lambda2,sigma,DRY_RUN);
		counter += counter3;

		counter2 += dryRun(occupations2);

		data_.resize(loc_);
		loc_=0;

		addAtTheFront(occupations);
		SizeType counter4=0;
		SizeType counter5=0;
		addAtTheMiddle(counter4,counter5,opPointers,lambda,lambda2,sigma,NORMAL_RUN);
		addAtTheBack(occupations2);

		// if daggers > non-daggers, result is zero
		if (counter!=counter2) {
			value_ = 0;
			return;
		}
	}

	SizeType findLocOfDiagOp(SizeType ind) const
	{
		SizeType counter = 0;
		SizeType j = 0;
		for (SizeType i=0;i<data_.size();i++) {
			if (data_[i].type == DIAGONAL) {
				if (counter==ind) return j;
				counter++;
				j++;
			}
			j++;
		}
		assert(false);
		return 0;
	}

	void removeNonCsOrDs()
	{
		typename PsimagLite::Vector<FreeOperator>::Type dataOld = data_;
		data_.clear();
		for (SizeType i=0;i<dataOld.size();i++) {
			SizeType type1 = dataOld[i].type;
			if (type1 == CREATION || type1 == DESTRUCTION)
				data_.push_back(dataOld[i]);
		}
	}

	SizeType size() const { return data_.size(); }

	const FreeOperator& operator[](SizeType i) const { return data_[i]; }

	void reverse()
	{
		std::reverse(data_.begin(), data_.end());
	}

	RealType operator()()
	{
		return value_;
	}

	void removePair(SizeType loc)
	{
		typename PsimagLite::Vector<FreeOperator>::Type::iterator itp = data_.begin();
		data_.erase(itp);
		if (loc == 0) throw PsimagLite::RuntimeError("removePair failed\n");
		--loc;
		itp = data_.begin()+loc;
		assert(itp < data_.end());
		data_.erase(itp);
	}

	bool notCreationOrDestruction(SizeType type1) const
	{
		return (type1 != CREATION && type1 != DESTRUCTION);
	}

private:

	SizeType dryRun(const VectorSizeType& occupations2)
	{
		SizeType counter = 0;
		for (SizeType i = 0; i < occupations2.size(); ++i) {
			SizeType x = (occupations2[i] == 0) ? 0 : 1;
			counter += x;
			loc_ += x;
		}

		return counter;
	}

	void addAtTheBack(const VectorSizeType& occupations2)
	{
		for (int i=occupations2.size()-1;i>=0;i--) {
			if (occupations2[i]==0) continue;

			FreeOperator fo;
			fo.lambda = i;
			fo.type = DESTRUCTION;
			data_[loc_++]=fo;
		}
	}

	void addAtTheFront(const VectorSizeType&  occupations)
	{
		for (SizeType i=0;i<occupations.size();++i) {
			if (occupations[i]==0) continue;
			FreeOperator fo;
			fo.lambda = i;
			fo.type = CREATION;
			data_[loc_++]=fo;
		}
	}

	void addAtTheMiddle(SizeType& counter,
	                    SizeType& counter2,
	                    const OpPointersType& opPointers,
	                    const IndexGeneratorType& lambda,
	                    const PermutationsType& lambda2,
	                    SizeType sigma,
	                    SizeType typeOfRun)
	{
		for (SizeType i=0;i<opPointers.size();i++) {
			FreeOperator fo;
			fo.type = opPointers[i].type;
			if (notCreationOrDestruction(fo.type)) {
				fo.lambda = 0;
				if (typeOfRun==NORMAL_RUN) data_[loc_++]=fo;
				else loc_++;
				continue;
			}
			if (opPointers[i].sigma!=sigma) continue;

			if (typeOfRun==DRY_RUN) {
				loc_++;
				continue;
			}
			if (fo.type==CREATION) {
				if (counter<lambda.size()) fo.lambda = lambda[counter];
				counter++;
			} else if (fo.type==DESTRUCTION) {
				if (counter2<lambda2.size()) fo.lambda = lambda2[counter2];
				counter2++;
			} else {
				fo.lambda = 0;
			}
			data_[loc_++]=fo;
		}
	}

	typename PsimagLite::Vector<FreeOperator>::Type data_;
	RealType value_;
	SizeType loc_;
}; // FreeOperators
} // namespace Dmrg 

/*@}*/
#endif // FREE_OPERATORS_H
