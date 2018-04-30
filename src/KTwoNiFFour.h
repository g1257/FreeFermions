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

/*! \file KTwoNiFFour.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef KTWONIFFOUR_H
#define KTWONIFFOUR_H
#include <stdexcept>
#include <cassert>
#include "Io/IoSimple.h"

namespace FreeFermions {

template<typename GeometryParamsType,typename MatrixType>
class KTwoNiFFour  {

	typedef std::pair<int,int> PairType;
	typedef typename GeometryParamsType::RealType RealType;
	typedef typename MatrixType::value_type FieldType;

	enum {TYPE_O,TYPE_C};
	enum {SUBTYPE_X,SUBTYPE_Y};
	enum {DIR_X,DIR_Y,DIR_XPY,DIR_XMY};

public:

	KTwoNiFFour(const GeometryParamsType& geometryParams)
	    : geometryParams_(geometryParams),signChange_(1)
	{
		PsimagLite::IoSimple::In io(geometryParams.filename);
		try {
			io.readline(signChange_,"SignChange=");
		} catch (std::exception& e) {
			io.rewind();
		}

		io.read(coHoppingsX_, "Connectors");
		io.read(coHoppingsY_, "Connectors");
		io.read(ooHoppingsXPY_, "Connectors");
		io.read(ooHoppingsXMY_, "Connectors");
	}

	void fillMatrix(MatrixType& t) const
	{
		SizeType sites = geometryParams_.sites;

		resizeAndZeroOutMatrix(t);
		for (SizeType i=0;i<sites;i++) {
			SizeType type1 = findTypeOfSite(i).first;
			for (SizeType j=0;j<sites;j++) {
				if (!connected(i,j)) continue;
				SizeType type2 = findTypeOfSite(j).first;
				if (type1==type2 && type1==TYPE_C) continue;
				if (type1==type2 && type1==TYPE_O) {
					orbitalsForO(t,i,j);
					continue;
				}
				SizeType newi = (type1==TYPE_C) ? i : j;
				SizeType newj = (type1==TYPE_C) ? j : i;
				orbitalsForCO(t,newi,newj);
			}
		}
		addPeriodicConnections(t);
	}

private:

	void addPeriodicConnections(MatrixType& t) const
	{
		if (!geometryParams_.isPeriodic[GeometryParamsType::DIRECTION_Y])
			return;
		
		SizeType n = geometryParams_.sites;

		// 0 --> N-1 Cu-O
		for (SizeType orb=0;orb<2;orb++)
			t(orb*n,n-1) = t(n-1,orb*n) = coOrbitals(DIR_X,orb);

		// 0 --> N-2 O-O
		for (SizeType orb1=0;orb1<2;orb1++) {
			for (SizeType orb2=0;orb2<2;orb2++) {
				t(index(0,orb1),index(n-2,orb2))=ooOrbitals(DIR_XPY,orb1,orb2);
				t(index(n-2,orb2),index(0,orb1))=ooOrbitals(DIR_XPY,orb1,orb2);
			}
		}

		// 0 --> N-3 O-O
		for (SizeType orb1=0;orb1<2;orb1++) {
			for (SizeType orb2=0;orb2<2;orb2++) {
				t(index(0,orb1),index(n-3,orb2))=ooOrbitals(DIR_XMY,orb1,orb2);
				t(index(n-3,orb2),index(0,orb1))=ooOrbitals(DIR_XMY,orb1,orb2);
			}
		}
	}

	void resizeAndZeroOutMatrix(MatrixType& t) const
	{
		SizeType rank = matrixRank();
		resizeAndZeroOut(t,rank,rank);
	}

	void resizeAndZeroOut(MatrixType& t,SizeType nrow,SizeType ncol) const
	{
		t.resize(nrow,ncol);
		for (SizeType i=0;i<nrow;i++)
			for(SizeType j=0;j<ncol;j++)
				t(i,j)=0;
	}

	SizeType matrixRank() const
	{
		SizeType sites = geometryParams_.sites;
		SizeType no = 0;
		SizeType nc = 0;
		for (SizeType i=0;i<sites;i++) {
			SizeType type1 = findTypeOfSite(i).first;
			if (type1==TYPE_C) nc++;
			else no++;
		}
		return 2*no+nc;
	}

	bool connected(SizeType i1,SizeType i2) const
	{
		if (i1==i2) return false;
		PairType type1 = findTypeOfSite(i1);
		PairType type2 = findTypeOfSite(i2);
		//! 4 possibilities
		//! c-c
		if (type1.first == type2.first && type1.first == TYPE_C)
			return false;
		//! o-o
		if (type1.first == type2.first) {
			assert(type1.first == TYPE_O);
			if (type1.second==type2.second) return false;
			SizeType newi1 = (type1.second==SUBTYPE_X) ? i1 : i2;
			SizeType newi2 = (type1.second==SUBTYPE_X) ? i2 : i1;
			if (newi1>newi2) {
				assert(newi1>=4);
				if (newi1-2==newi2 || newi1-3==newi2) return true;
				return false;
			}
			if (newi2-1==newi1) return true;
			if (newi2>=2 && newi2-2==newi1) return true;
			return false;
		}
		//! o-c or c-o
		SizeType newi1 = (type1.first==TYPE_O) ? i1 : i2;
		SizeType newi2 = (type1.first==TYPE_O) ? i2 : i1;
		assert(newi2>=3);
		if (newi2-1==newi1 || newi2-2==newi1 || newi2-3==newi1 || newi2+1==newi1)
			return true;
		return false;
	}

	void orbitalsForO(MatrixType& t,SizeType i1,SizeType i2) const
	{
		SizeType dir = calcDir(i1,i2);
		RealType sign = signChange(i1,i2);
		for (SizeType orb1=0;orb1<2;orb1++) {
			for (SizeType orb2=0;orb2<2;orb2++) {
				t(index(i1,orb1),index(i2,orb2))=ooOrbitals(dir,orb1,orb2)*sign;
			}
		}
	}

	void orbitalsForCO(MatrixType& t,SizeType i1,SizeType i2) const
	{
		SizeType dir = calcDir(i1,i2);
		RealType sign = signChange(i1,i2);
		for (SizeType orb=0;orb<2;orb++)
			t(index(i2,orb),index(i1)) = t(index(i1),index(i2,orb)) =
			        coOrbitals(dir,orb)*sign;
	}

	SizeType index(SizeType i,SizeType orb) const
	{
		if (orb==0) return i;
		SizeType sites = geometryParams_.sites;
		SizeType tmp = (i+1)/4;
		assert(sites+i>=tmp);
		return sites+i-tmp;
	}

	SizeType index(SizeType i) const
	{
		return i;
	}

	int signChange(SizeType i1,SizeType i2) const
	{
		// change the sign of Cu-O
		SizeType newi1 = std::min(i1,i2);
		SizeType newi2 = std::max(i1,i2);
		PairType type1 = findTypeOfSite(newi1);
		PairType type2 = findTypeOfSite(newi2);
		int sign1 = 1;
		if (type1.first!=type2.first) {

			int diff =  newi2-newi1;
			assert(diff==1 || diff==2 || diff==3);
			if (diff<2) sign1 = -1;
		}

		return (isInverted(i1) || isInverted(i2)) ? signChange_*sign1 : sign1;
	}

	bool isInverted(SizeType i) const
	{
		SizeType j = i+4;
		return ((j%8)==0);
	}

	FieldType ooOrbitals(SizeType dir,SizeType orb1,SizeType orb2) const
	{
		return (dir==DIR_XPY) ? ooHoppingsXPY_(orb1,orb2)
		                      : ooHoppingsXMY_(orb1,orb2);
	}

	FieldType coOrbitals(SizeType dir,SizeType orb) const
	{
		return (dir==DIR_X) ? coHoppingsX_(orb,0) : coHoppingsY_(orb,0);
	}

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType i1,SizeType i2) const
	{
		PairType type1 = findTypeOfSite(i1);
		PairType type2 = findTypeOfSite(i2);
		//! o-o
		if (type1.first == type2.first && type1.first == TYPE_O) {
			assert(type1.second!=type2.second);
			SizeType newi1 = (type1.second==SUBTYPE_X) ? i1 : i2;
			SizeType newi2 = (type1.second==SUBTYPE_X) ? i2 : i1;
			if (newi1>newi2) {
				assert(newi1>=4);
				SizeType distance = newi1-newi2;
				if (distance==2) return DIR_XPY;
				assert(distance==3);
				return DIR_XMY;
			}
			SizeType distance = newi2-newi1;
			if (distance==1)  return DIR_XPY;
			assert(distance==2);
			return DIR_XMY;
		}
		//! o-c or c-o
		SizeType newi1 = (type1.first==TYPE_O) ? i1 : i2;
		type1 = findTypeOfSite(newi1);
		return (type1.second==SUBTYPE_X) ? DIR_X : DIR_Y;
	}

	//! Given as a pair, first number is the type,
	//! If first number == TYPE_C then second number is bogus
	//! If first number == TYPE_O then second number is the subtype
	PairType findTypeOfSite(SizeType site) const
	{
		SizeType sitePlusOne = site + 1;
		SizeType r = sitePlusOne%4;
		if (r==0) return PairType(TYPE_C,0);

		if (r==1) return PairType(TYPE_O,SUBTYPE_X);
		return PairType(TYPE_O,SUBTYPE_Y);
	}

	const GeometryParamsType& geometryParams_;
	int signChange_;
	MatrixType ooHoppingsXPY_,ooHoppingsXMY_;
	MatrixType coHoppingsX_,coHoppingsY_;
}; // class KTwoNiFFour
} // namespace FreeFermions 

/*@}*/
#endif // KTWONIFFOUR_H
