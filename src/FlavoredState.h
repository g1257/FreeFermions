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
#include "Vector.h"
#include "BitManip.h"

namespace FreeFermions {

class Tstorage {

public:

	typedef PsimagLite::Vector<bool>::Type LevelsType;
	typedef PsimagLite::Vector<unsigned char>::Type VectorUcharType;
	typedef PsimagLite::Vector<VectorUcharType>::Type VectorVectorUcharType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	static void init(SizeType dof, SizeType blockSize, SizeType threadNum)
	{
		if (dof_ == dof) return;

		for (SizeType i = 0; i < x_.size(); ++i)
			x_[i].resize(100);

		dof_ = dof;
		blockSize_ = blockSize;
		assert(blockSize_ < x_[threadNum].size());
	}

	static void set(const LevelsType& portion, SizeType threadNum)
	{
		SizeType c = 0;
		SizeType j = 0;
		bytes_[threadNum] = 0;
		for (SizeType i = 0; i < portion.size(); ++i) {
			bool b = portion[i];
			SizeType mask = (1<<j);
			if (b)
				c |= mask;
			else
				c &= (~mask);
			j++;
			if (j == 8) {
				assert(bytes_[threadNum] < x_[threadNum].size());
				SizeType offset = bytes_[threadNum];
				x_[threadNum][offset] = c;
				bytes_[threadNum]++;
				j = 0;
				c = 0;
			}
		}

		if (j > 0) {
			assert(bytes_[threadNum] < x_[threadNum].size());
			SizeType offset = bytes_[threadNum];
			x_[threadNum][offset] = c;
			bytes_[threadNum]++;
		}
	}

	static void set(const VectorUcharType& d,
	                SizeType flavor,
	                SizeType threadNum)
	{
		SizeType index = flavor * dof_;
		for (SizeType i = 0; i < blockSize_; ++i)
			x_[threadNum][i] = d[index++];
		bytes_[threadNum] = blockSize_;
	}

	static SizeType numberOfDigits(SizeType threadNum)
	{
		SizeType sum = 0;
		for (SizeType i = 0; i < blockSize_; ++i)
			sum += PsimagLite::BitManip::count(x_[threadNum][i]);
		return sum;
	}

	static SizeType size(SizeType threadNum) {return bytes_[threadNum]; }

	static unsigned char byte(SizeType i, SizeType threadNum)
	{
		assert(i < bytes_[threadNum]);
		assert(i < x_[threadNum].size());
		return  x_[threadNum][i];
	}

	static SizeType statesBetween(SizeType lambda, SizeType threadNum)
	{
		SizeType index = getIndex(0,lambda);
		SizeType sum = 0;
		for (SizeType i = 0; i < index; ++i)
			sum += PsimagLite::BitManip::count(x_[threadNum][i]);
		SizeType offset = lambda % 8;
		offset++;
		SizeType mask = (1<<offset);
		mask--;
		SizeType x = x_[threadNum][index] & mask;
		sum += PsimagLite::BitManip::count(x);
		return sum;
	}

	static bool bitAt(SizeType lambda, SizeType threadNum)
	{
		SizeType index = getIndex(0,lambda);
		SizeType offset = lambda % 8;
		SizeType mask = (1<<offset);
		return ((x_[threadNum][index] & mask) > 0);
	}

	static SizeType getIndex(SizeType flavor, SizeType lambda)
	{
		SizeType index = flavor * dof_;
		return index + static_cast<SizeType>(lambda/8);
	}

private:

	static SizeType dof_;
	static SizeType blockSize_;
	static VectorSizeType bytes_;
	static VectorVectorUcharType x_;
};

// All interactions == 0
template<typename OperatorType>
class FlavoredState {

	//static SizeType const SPIN_UP=0,SPIN_DOWN=1;
	typedef FlavoredState<OperatorType> ThisType;
	static int const FERMION_SIGN = -1;
	enum {CREATION = OperatorType::CREATION,
	       DESTRUCTION = OperatorType::DESTRUCTION};

	typedef Tstorage TstorageType;
	typedef Tstorage::VectorUcharType VectorUcharType;

public:

	typedef Tstorage::LevelsType LevelsType;

	FlavoredState(SizeType dof,SizeType size, SizeType threadNum)
	: dof_(dof),threadNum_(threadNum)
	{
		blockSize_ = static_cast<SizeType>(size/8);
		if (size % 8 != 0) blockSize_++;
		data_.resize(blockSize_*dof,0);
		TstorageType::init(dof_,blockSize_,threadNum_);
	}

	void pushInto(SizeType sigma,const LevelsType& portion)
	{
		TstorageType::set(portion,threadNum_);
		SizeType index = sigma*blockSize_;
		for (SizeType i = 0; i < TstorageType::size(threadNum_); ++i) {
			data_[index++] = TstorageType::byte(i,threadNum_);
		}
	}

	int apply(SizeType label,SizeType flavor,SizeType lambda)
	{
		assert(flavor < dof_);
		int interSign = (calcInterElectrons(flavor) %2) ? 1 : FERMION_SIGN;
		bool b = false;
		int s = applyInternal(b,label,flavor,lambda)*interSign;
		setData(flavor,lambda,b);
		return s;
	}

	SizeType flavors() const { return dof_; }

	bool equalEqual(const ThisType& other) const
	{
		assert(data_.size() == other.data_.size());
		for (SizeType i = 0; i < data_.size(); ++i) {
			SizeType a = data_[i];
			SizeType b = other.data_[i];
			if (a != b) return false;
		}
		return true;
	}

	bool lessThan(const ThisType& other) const
	{
		assert(data_.size() == other.data_.size());
		for (SizeType i = 0; i < data_.size(); ++i) {
			SizeType  a = data_[i];
			SizeType b = other.data_[i];
			if (a > b) return false;
			if (a < b) return true;
		}

		return false;
	}

private:

	void setData(SizeType flavor,SizeType lambda,bool b)
	{
		SizeType index = Tstorage::getIndex(flavor,lambda);
		SizeType offset = lambda % 8;
		SizeType mask = (1<<offset);
		assert(index < data_.size());
		if (b)
			data_[index] |= mask;
		else
			data_[index] &= (~mask);
	}

	int applyInternal(bool& result,SizeType label,SizeType flavor,SizeType lambda)
	{
		TstorageType::set(data_,flavor,threadNum_);
		SizeType nflips = TstorageType::statesBetween(lambda,threadNum_);
		if (label == CREATION) {
			if (TstorageType::bitAt(lambda,threadNum_)) return 0; // can't create, there's already one
			result = true;
		} else if (label == DESTRUCTION) {
			if (!TstorageType::bitAt(lambda,threadNum_)) return 0; // can't destroy, there's nothing
			result =false;
		} else {
			throw std::runtime_error("FlavoredState::applyInternal()\n");
		}

		if (nflips ==0 || nflips % 2 ==0) return 1;
		return FERMION_SIGN;
	}

	SizeType calcInterElectrons(SizeType flavor)
	{
		SizeType sum = 0;
		for (SizeType flavor2 = 0; flavor2 < flavor; flavor2++) {
			TstorageType::set(data_,flavor2,threadNum_);
			sum += TstorageType::numberOfDigits(threadNum_);
		}

		return sum;
	}

	SizeType dof_;
	SizeType threadNum_;
	SizeType blockSize_;
	VectorUcharType data_;
}; // class FlavoredState

template<typename U>
inline bool operator==(const FlavoredState<U>& v1,const FlavoredState<U>& v2)
{
	return v1.equalEqual(v2);
}

template<typename U>
inline bool operator<(const FlavoredState<U>& v1,const FlavoredState<U>& v2)
{
	return v1.lessThan(v2);
}

Tstorage::VectorVectorUcharType Tstorage::x_(100);

Tstorage::VectorSizeType Tstorage::bytes_(100);

SizeType Tstorage::dof_ = 0;

SizeType Tstorage::blockSize_ = 0;
} // namespace Dmrg 

/*@}*/
#endif

