/*
Copyright  2009-2012, UT-Battelle, LLC
All rights reserved

[FreeFermions, Version 1.0.0]
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

/*! \file GeometryLibrary.h
 *
 * handy geometries
 *
 */
#ifndef GEOMETRY_LIB_H
#define GEOMETRY_LIB_H
#include "IoSimple.h" // in psimaglite
#include "Matrix.h" // in psimaglite
#include <cassert>
#include "KTwoNiFFour.h"

namespace FreeFermions {

template<typename MatrixType,typename GeometryParamsType_>
class GeometryLibrary {

public:

	enum DecayEnum {DECAY_NONE, DECAY_0, DECAY_1};

	typedef GeometryParamsType_ GeometryParamsType;
	typedef typename MatrixType::value_type FieldType;
	typedef typename PsimagLite::Real<FieldType>::Type RealType;

	enum {CHAIN=GeometryParamsType::CHAIN,
	      LADDER=GeometryParamsType::LADDER,
	      FEAS=GeometryParamsType::FEAS,
	      KTWONIFFOUR=GeometryParamsType::KTWONIFFOUR,
	      FEAS1D=GeometryParamsType::FEAS1D,
	      CHAIN_EX = GeometryParamsType::CHAIN_EX,
	      LADDER_BATH = GeometryParamsType::LADDER_BATH,
	      STAR = GeometryParamsType::STAR
	     };

	enum {DIRECTION_X   = GeometryParamsType::DIRECTION_X,
	      DIRECTION_Y   = GeometryParamsType::DIRECTION_Y,
	      DIRECTION_XPY = GeometryParamsType::DIRECTION_XPY,
	      DIRECTION_XMY = GeometryParamsType::DIRECTION_XMY
	     };

	GeometryLibrary(GeometryParamsType& geometryParams,DecayEnum decay = DECAY_NONE)
	: geometryParams_(geometryParams),decay_(decay)
	{
		switch (geometryParams.type) {
		case CHAIN:
			setGeometryChain();
			break;
		case CHAIN_EX:
			setGeometryChainEx();
			break;
		case LADDER:
			setGeometryLadder();
			break;
		case LADDER_BATH:
			setGeometryLadder();
			bathify();
			break;
		case FEAS:
		case FEAS1D:
		case STAR:
			setGeometryFeAs();
			break;
		case KTWONIFFOUR:
			setGeometryKniffour();
			break;
		default:
			assert(false);
		}
		typename PsimagLite::Vector<RealType>::Type v;
		readPotential(v,geometryParams_.filename);
		addPotential(v);
	}

	template<typename ComplexType>
	void fourierTransform(typename PsimagLite::Vector<ComplexType>::Type& dest,
	                      const MatrixType& src,
	                      SizeType leg) const
	{
		SizeType n = src.n_row();
		if (n!=geometryParams_.sites)
			throw std::runtime_error("src must have the same number of sites as lattice\n");
		PsimagLite::Matrix<ComplexType> B(n,n);
		getFourierMatrix(B,leg);
		for (SizeType k=0; k<n; k++) {
			ComplexType sum = 0.0;
			for (SizeType i=0; i<n; i++) {
				for (SizeType j=0; j<n; j++) {
					sum += std::conj(B(i,k)) * src(i,j) * B(j,k);
				}
			}
			dest[k] = sum;
		}
	}

	const MatrixType& matrix() const
	{
		return t_;
	}

	PsimagLite::String name() const {

		switch (geometryParams_.type) {
		case CHAIN:
			return "chain";
			break;
		case CHAIN_EX:
			return "chainEx";
			break;
		case LADDER:
			return "ladder";
			break;
		case LADDER_BATH:
			return "ladderbath";
			break;
		case FEAS:
			return "feas";
			break;
		case KTWONIFFOUR:
			return "kniffour";
			break;
		case FEAS1D:
			return "feas1d";
			break;
		case STAR:
			return "star";
			break;
		default:
			assert(false);
		}
		return "unknown";
	}

	void addPotential(const RealType& value)
	{
		SizeType n = t_.n_row();

		for (SizeType i=0; i<n; i++) t_(i,i) += value;
	}

	template<typename MType,typename PType>
	friend std::ostream& operator<<(std::ostream& os,
	                                const GeometryLibrary<MType,PType>& gl);

private:

	template<typename VectorLikeType>
	void addPotential(const VectorLikeType& p2)
	{
		if (p2.size()==0) return;

		VectorLikeType p = p2;
		if (p.size()!=t_.n_row()) {
			SizeType halfSize = SizeType(p2.size()/2);
			p.resize(halfSize);
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "addPotential(...): resizing potential to " + ttos(t_.n_row());
			str += " from " + ttos(p2.size()) + "\n";
			std::cerr<<str;
		}
		if (p.size()!=t_.n_row()) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "addPotential(...): expecting " + ttos(t_.n_row());
			str += " numbers but " + ttos(p.size()) + " found instead.\n";
			throw std::runtime_error(str.c_str());
		}
		for (SizeType i=0; i<p.size(); i++) t_(i,i) = p[i];
	}

	void setGeometryFeAs()
	{
		typename PsimagLite::Vector<MatrixType>::Type t;
		SizeType edof = geometryParams_.orbitals;
		typename PsimagLite::Vector<FieldType>::Type oneSiteHoppings;
		SizeType dirs = (geometryParams_.type==FEAS1D || geometryParams_.type==STAR) ? 1 : 4;
		readOneSiteHoppings(oneSiteHoppings,geometryParams_.filename,dirs);
		if (oneSiteHoppings.size()!=dirs*edof*edof) {
			throw std::runtime_error("Wrong number of hoppings\n");
		}
		SizeType sites = geometryParams_.sites;
		for (SizeType i=0; i<edof*edof; i++) { // orbital*orbital cases: aa ab ba bb etc
			MatrixType oneT(sites,sites);
			if (geometryParams_.type == STAR)
				setGeometryFeAsStar(oneT,i,oneSiteHoppings);
			else
				setGeometryFeAs(oneT,i,oneSiteHoppings);
			reorderLadderX(oneT,geometryParams_.leg);
			assert(isHermitian(oneT,true));
			t.push_back(oneT);
		}
		resizeAndZeroOut(edof*sites,edof*sites);
		for (SizeType orbitalPair=0; orbitalPair<edof*edof; orbitalPair++) {
			SizeType orb1 = (orbitalPair % geometryParams_.orbitals);
			SizeType orb2 = orbitalPair/geometryParams_.orbitals;
			for (SizeType i=0; i<sites; i++)
				for (SizeType j=0; j<sites; j++)
					t_(i+orb1*sites,j+orb2*sites) += t[orbitalPair](i,j);
		}
	}

	void setGeometryKniffour()
	{
		KTwoNiFFour<GeometryParamsType,MatrixType> ktwoniffour(geometryParams_);
		ktwoniffour.fillMatrix(t_);
	}

	void setGeometryChain()
	{
		SizeType sites = geometryParams_.sites;
		resizeAndZeroOut(sites,sites);
		for (SizeType i=0; i<sites; i++) {
			for (SizeType j=0; j<sites; j++) {
				if (i-j==1 || j-i==1) t_(i,j) = geometryParams_.hopping[0];
			}
		}
		if (geometryParams_.isPeriodic[DIRECTION_X])
			t_(0,sites-1) = t_(sites-1,0) = geometryParams_.hopping[0];
	}

	void setGeometryChainEx()
	{
		typename PsimagLite::Vector<MatrixType>::Type t;
		SizeType edof = geometryParams_.orbitals;
		typename PsimagLite::Vector<FieldType>::Type oneSiteHoppings;
		SizeType dirs = 2;
		readOneSiteHoppings(oneSiteHoppings,geometryParams_.filename,dirs);
		if (oneSiteHoppings.size()!=dirs*edof*edof) {
			throw std::runtime_error("Wrong number of hoppings\n");
		}
		SizeType sites = geometryParams_.sites;
		for (SizeType i=0; i<edof*edof; i++) {
			MatrixType oneT(sites,sites);
			setGeometryChainEx(oneT,i,oneSiteHoppings);
			assert(isHermitian(oneT,true));
			t.push_back(oneT);
		}

		resizeAndZeroOut(edof*sites,edof*sites);
		for (SizeType orbitalPair=0; orbitalPair<edof*edof; orbitalPair++) {
			SizeType orb1 = (orbitalPair % geometryParams_.orbitals);
			SizeType orb2 = orbitalPair/geometryParams_.orbitals;
			for (SizeType i=0; i<sites; i++)
				for (SizeType j=0; j<sites; j++)
					t_(i+orb1*sites,j+orb2*sites) += t[orbitalPair](i,j);
		}
	}

	void setGeometryLadder()
	{
		SizeType leg = geometryParams_.leg;
		if (leg<2)
			throw std::runtime_error("GeometryLibrary:: ladder must have leg>1\n");
		assert(!geometryParams_.isPeriodic[DIRECTION_Y] || leg>2);
		SizeType sites = geometryParams_.sites;
		assert(geometryParams_.hopping.size()==sites/2 + sites -2);
		resizeAndZeroOut(sites,sites);
		for (SizeType i=0; i<sites; i++) {
			PsimagLite::Vector<SizeType>::Type v;
			ladderFindNeighbors(v,i,leg,geometryParams_.isPeriodic);
			for (SizeType k=0; k<v.size(); k++) {
				SizeType j = v[k];
				SizeType ijmin = (i < j) ? i : j;
				SizeType ijminOver2 = static_cast<SizeType>(ijmin/2);

				if (ijmin & 1) {
					ijmin = (ijmin-1)/2  + sites/2 -1  ;
				} else {
					ijmin = ijminOver2;
				}

				FieldType tmp = (ladderSameColumn(i,j,leg))?
				                geometryParams_.hopping[ijminOver2 + sites-2] :
				                geometryParams_.hopping[ijmin];
				t_(i,j) = tmp;
				t_(j,i) = std::conj(t_(i,j));
			}
		}
	}

	void ladderFindNeighbors(PsimagLite::Vector<SizeType>::Type& v,
	                         SizeType i,
	                         SizeType leg,
	                         const PsimagLite::Vector<bool>::Type& isPeriodic)
	{
		SizeType sites = geometryParams_.sites;
		SizeType lengthX = static_cast<SizeType>(sites/leg);

		v.clear();
		SizeType k = i + 1;
		if (k<sites && ladderSameColumn(k,i,leg)) v.push_back(k);
		k = i + leg;
		if (k<sites)  v.push_back(k);

		if (leg>2 && isPeriodic[DIRECTION_Y]) {
			if (i < lengthX) {
				k = i + (leg-1)*lengthX;
				if (k<sites) v.push_back(k);
			}
		}

		if (isPeriodic[DIRECTION_X]) {
			if (i == 0 || i % lengthX == 0) {
				k = i +lengthX -1;
				if (k<sites) v.push_back(k);
			}
		}

		// careful with substracting because i and k are unsigned!
		if (i==0) return;
		k = i-1;
		if (ladderSameColumn(i,k,leg)) v.push_back(k);

		if (i<leg) return;
		k= i-leg;
		v.push_back(k);
	}

	bool ladderSameColumn(SizeType i,SizeType k,SizeType leg)
	{
		SizeType i1 = i/leg;
		SizeType k1 = k/leg;
		if (i1 == k1) return true;
		return false;
	}

	// only 2 orbitals supported
	void setGeometryFeAs(MatrixType& t,
	                     SizeType orborb,
	                     const typename PsimagLite::Vector<FieldType>::Type& oneSiteHoppings)
	{
		SizeType sites = geometryParams_.sites;
		SizeType leg = geometryParams_.leg;
		SizeType orbitalsSquared = geometryParams_.orbitals * geometryParams_.orbitals;
		FieldType tx = oneSiteHoppings[orborb+DIRECTION_X*orbitalsSquared];
		SizeType lengthx  = sites/leg;
		if (sites%leg!=0) throw std::runtime_error("Leg must divide number of sites.\n");
		for (SizeType j=0; j<leg; j++) {
			for (SizeType i=0; i<lengthx; i++) {
				if (i+1<lengthx) t(i+1+j*lengthx,i+j*lengthx) = t(i+j*lengthx,i+1+j*lengthx) = tx;
				if (i>0) t(i-1+j*lengthx,i+j*lengthx) = t(i+j*lengthx,i-1+j*lengthx) = tx;
			}
			if (geometryParams_.isPeriodic[GeometryParamsType::DIRECTION_X])
				t(j*lengthx,lengthx-1+j*lengthx) = t(lengthx-1+j*lengthx,j*lengthx) = tx;
		}

		if (geometryParams_.type==FEAS1D) return;

		FieldType ty = oneSiteHoppings[orborb+DIRECTION_Y*orbitalsSquared];
		for (SizeType i=0; i<lengthx; i++) {
			for (SizeType j=0; j<leg; j++) {
				if (j>0) t(i+(j-1)*lengthx,i+j*lengthx) = t(i+j*lengthx,i+(j-1)*lengthx) = ty;
				if (j+1<leg) t(i+(j+1)*lengthx,i+j*lengthx) = t(i+j*lengthx,i+(j+1)*lengthx) = ty;
			}
			if (geometryParams_.isPeriodic[GeometryParamsType::DIRECTION_Y])
				t(i,i+(leg-1)*lengthx) = t(i+(leg-1)*lengthx,i) = ty;
		}
		FieldType txpy = oneSiteHoppings[orborb+DIRECTION_XPY*orbitalsSquared];
		FieldType txmy = oneSiteHoppings[orborb+DIRECTION_XMY*orbitalsSquared];
		for (SizeType i=0; i<lengthx; i++) {
			for (SizeType j=0; j<leg; j++) {
				if (j+1<leg && i+1<lengthx)
					t(i+1+(j+1)*lengthx,i+j*lengthx) = t(i+j*lengthx,i+1+(j+1)*lengthx) = txpy;
				if (i+1<lengthx && j>0)
					t(i+1+(j-1)*lengthx,i+j*lengthx) = t(i+j*lengthx,i+1+(j-1)*lengthx) = txmy;
				if (!geometryParams_.isPeriodic[GeometryParamsType::DIRECTION_X] || i>0)
					continue;
				if (j+1<leg)
					t((j+1)*lengthx,lengthx-1+j*lengthx) = t(lengthx-1+j*lengthx,(j+1)*lengthx) = txpy;
				if (j>0)
					t((j-1)*lengthx,lengthx-1+j*lengthx) = t(lengthx-1+j*lengthx,(j-1)*lengthx) = txmy;
			}
			if (!geometryParams_.isPeriodic[GeometryParamsType::DIRECTION_X])
				continue;
			if (i+1<lengthx) {
				t(i+1,i+(leg-1)*lengthx) = t(i+(leg-1)*lengthx,i+1) = txpy;
				t(i+1+(leg-1)*lengthx,i) = t(i,i+1+(leg-1)*lengthx) = txmy;
			}
			if (i>0) continue;
			t(0,lengthx-1+(leg-1)*lengthx) = t(lengthx-1+(leg-1)*lengthx,0) = txpy;
			t(0+(leg-1)*lengthx,lengthx-1) = t(lengthx-1,0+(leg-1)*lengthx) = txmy;
		}

	}

	void setGeometryFeAsStar(MatrixType& t,
	                         SizeType orborb,
	                         const typename PsimagLite::Vector<FieldType>::Type& oneSiteHoppings)
	{
		assert(geometryParams_.type == STAR);
		SizeType sites = geometryParams_.sites;
		SizeType orbitalsSquared = geometryParams_.orbitals * geometryParams_.orbitals;
		FieldType tx = oneSiteHoppings[orborb+DIRECTION_X*orbitalsSquared];

		for (SizeType i=1; i<sites; ++i) t(0,i) = t(i,0) = tx;
	}

	void setGeometryChainEx(MatrixType& t,
	                        SizeType orborb,
	                        const typename PsimagLite::Vector<FieldType>::Type& oneSiteHoppings)
	{
		SizeType sites = geometryParams_.sites;
		SizeType orbitalsSquared = geometryParams_.orbitals * geometryParams_.orbitals;
		for (SizeType i=0; i<sites; i++) {
			for (SizeType j=0; j<sites; j++) {
				if (i-j==1 || j-i==1) t(i,j) = oneSiteHoppings[orborb+0*orbitalsSquared];
				if (i-j==2 || j-i==2) {
					if ((i & 1) == 0) t(i,j) = oneSiteHoppings[orborb+1*orbitalsSquared];
				}
			}
		}
	}

	// from 0--1--2--...
	//      N-N+1-N+2--..
	// into:
	//      0--2--4--
	//      1--3--5--
	//
	void reorderLadderX(MatrixType& told,SizeType leg)
	{
		SizeType sites = geometryParams_.sites;
		MatrixType tnew(told.n_row(),told.n_col());
		for (SizeType i=0; i<sites; i++) {
			for (SizeType j=0; j<sites; j++) {
				SizeType i2 = reorderLadderX(i,leg);
				SizeType j2 = reorderLadderX(j,leg);
				tnew(i2,j2) = told(i,j);
			}
		}
		told = tnew;
	}

	SizeType reorderLadderX(SizeType i,SizeType leg)
	{
		SizeType sites = geometryParams_.sites;
		SizeType lengthx  = sites/leg;
		// i = x + y*lengthx;
		SizeType x = i%lengthx;
		SizeType y = i/lengthx;
		return y + x*leg;
	}

	void readOneSiteHoppings(typename PsimagLite::Vector<FieldType>::Type& v,
	                         const PsimagLite::String& filename,
	                         SizeType dirs)
	{
		typename PsimagLite::IoSimple::In io(filename);
		v.clear();
		for (SizeType i=0; i<dirs; i++) {
			MatrixType m;
			io.readMatrix(m,"Connectors");
			if (m.n_row()!=m.n_col() || m.n_row()!=geometryParams_.orbitals) {
				PsimagLite::String str(__FILE__);
				str += " " + ttos(__LINE__) + "\n";
				str += "Error in input file "+filename+" label: Connectors\n";
				throw std::runtime_error(str.c_str());
			}
			appendToVector(v,m);
		}
	}

	void appendToVector(typename PsimagLite::Vector<FieldType>::Type& v,
	                    const MatrixType& m) const
	{
		for (SizeType i=0; i<m.n_row(); i++) {
			for (SizeType j=0; j<m.n_col(); j++) {
				v.push_back(m(i,j));
			}
		}
	}

	template<typename MatrixComplexType>
	void getFourierMatrix(MatrixComplexType& B,SizeType leg) const
	{
		typedef typename MatrixComplexType::value_type ComplexType;
		if (geometryParams_.geometryType!=FEAS)
			throw std::runtime_error("getFourierMatrix: unsupported\n");
		SizeType n = B.n_row();
		SizeType lengthx = n/leg;
		if (n%leg !=0) throw std::runtime_error("Leg must divide number of sites for FEAS\n");
		for (SizeType i=0; i<n; i++) {
			SizeType rx = i % lengthx;
			SizeType ry = i/lengthx;
			for (SizeType k=0; k<n; k++) {
				SizeType kx = k % lengthx;
				SizeType ky = k / lengthx;
				RealType tmpx = rx*kx*2.0*M_PI/lengthx;
				RealType tmpy = ry*ky*2.0*M_PI/leg;
				RealType tmp1 = cos(tmpx+tmpy);
				RealType tmp2 = sin(tmpx+tmpy);
				B(i,k) = ComplexType(tmp1,tmp2);
			}
		}
	}

	void resizeAndZeroOut(SizeType nrow,SizeType ncol)
	{
		t_.resize(nrow,ncol);
		for (SizeType i=0; i<nrow; i++)
			for (SizeType j=0; j<ncol; j++)
				t_(i,j)=0;
	}

	void readPotential(typename PsimagLite::Vector<RealType>::Type& v,
	                   const PsimagLite::String& filename)
	{
		typename PsimagLite::Vector<RealType>::Type w;
		PsimagLite::IoSimple::In io(filename);

		try {
			io.read(w,"potentialT");
		} catch (std::exception& e) {
			std::cerr<<"INFO: No PotentialT in file "<<filename<<"\n";
		}
		io.rewind();

		io.read(v,"potentialV");
		if (decay_ == DECAY_0) {
			return;
		}

		if (w.size()==0) return;

		if (decay_ == DECAY_1) {
			v = w;
			return;
		}

		if (v.size()>w.size()) v.resize(w.size());
		for (SizeType i=0; i<w.size(); i++)
			v[i] += w[i];
	}

	void bathify()
	{
		SizeType sites = t_.n_row();
		SizeType nb =geometryParams_.bathSitesPerSite;
		SizeType nnew = sites*(1+nb);
		MatrixType tnew(nnew,nnew);
		SizeType sitesOver2 = static_cast<SizeType>(sites*0.5);
		SizeType firstC = sitesOver2*nb;
		for (SizeType i=0; i<t_.n_row(); i++) {
			for (SizeType j=0; j<t_.n_col(); j++)
				tnew(i+firstC,j+firstC) = t_(i,j);
			for (SizeType j=0; j<nb; j++) {
				SizeType k = j*sitesOver2 + i;
				SizeType k2 = j*sites + i;
				if (i >= sitesOver2) k += sitesOver2 * (nb +1);
				assert(k2 < geometryParams_.tb.size());
				tnew(i+firstC,k) = tnew(k,i+firstC) =geometryParams_.tb[k2];
			}
		}
		t_ = tnew;
	}

	const GeometryParamsType& geometryParams_;
	DecayEnum decay_;
	MatrixType t_;

}; // GeometryLibrary

template<typename MatrixType,typename ParamsType>
std::ostream& operator<<(std::ostream& os,
                         const GeometryLibrary<MatrixType,ParamsType>& gl)
{
	os<<gl.t_;
	os<<"GeometryName="<<gl.name()<<"\n";
	return os;
}
} // namespace FreeFermions

/*@}*/
#endif //GEOMETRY_LIB_H

