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

		typedef GeometryParamsType_ GeometryParamsType;
		typedef typename MatrixType::value_type RealType;

		enum {CHAIN=GeometryParamsType::CHAIN,
			  LADDER=GeometryParamsType::LADDER,
			  FEAS=GeometryParamsType::FEAS,
			  KTWONIFFOUR=GeometryParamsType::KTWONIFFOUR,
			  FEAS1D=GeometryParamsType::FEAS1D};

		enum {DIRECTION_X   = GeometryParamsType::DIRECTION_X,
			  DIRECTION_Y   = GeometryParamsType::DIRECTION_Y,
			  DIRECTION_XPY = GeometryParamsType::DIRECTION_XPY,
			  DIRECTION_XMY = GeometryParamsType::DIRECTION_XMY};
		
		GeometryLibrary(GeometryParamsType& geometryParams) 
		: geometryParams_(geometryParams)
		{
			switch (geometryParams.type) {
			case CHAIN:
				setGeometryChain();
				break;
			case LADDER:
				setGeometryLadder();
				break;
			case FEAS:
			case FEAS1D:
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

		void bathify(const typename PsimagLite::Vector<RealType>::Type& tb)
		{
			size_t sites = geometryParams_.sites;
			assert(t_.size()!=1);
			if (sites!=t_.n_row() || sites!=t_.n_col())
				throw std::runtime_error("GeometryLibrary::bathify(...)\n");
			size_t nb = tb.size();
			size_t nnew = sites*(1+nb);
			MatrixType tnew(nnew,nnew);
			for (size_t i=0;i<t_.n_row();i++)
			{
				for (size_t j=0;j<t_.n_col();j++)
					tnew(i,j) = t_(i,j);
				for (size_t j=0;j<nb;j++) {
					size_t k = sites + j + nb*i;
					tnew(i,k) = tnew(k,i) = tb[j];
				}
			}
			t_ = tnew;
		}

		template<typename ComplexType>
		void fourierTransform(typename PsimagLite::Vector<ComplexType>::Type& dest,
		                      const MatrixType& src,
		                      size_t leg) const
		{
			size_t n = src.n_row();
			if (n!=geometryParams_.sites) throw std::runtime_error("src must have the same number of sites as lattice\n");
			PsimagLite::Matrix<ComplexType> B(n,n);
			getFourierMatrix(B,leg);
			for (size_t k=0;k<n;k++) {
				ComplexType sum = 0.0;
				for (size_t i=0;i<n;i++) {
					for (size_t j=0;j<n;j++) {
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

		PsimagLite::String name() const
		{

			switch (geometryParams_.type) {
			case CHAIN:
				return "chain";
				break;
			case LADDER:
				return "ladder";
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
			default:
				assert(false);
			}
			return "unknown";
		}

		void addPotential(const RealType& value)
		{
			size_t n = t_.n_row();

			for (size_t i=0;i<n;i++) t_(i,i) += value;
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
				size_t halfSize = size_t(p2.size()/2);
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
			for (size_t i=0;i<p.size();i++) t_(i,i) = p[i];
		}

		void setGeometryFeAs()
		{
			typename PsimagLite::Vector<MatrixType>::Type t;
			size_t edof = geometryParams_.orbitals;
			typename PsimagLite::Vector<RealType>::Type oneSiteHoppings;
			size_t dirs = (geometryParams_.type==FEAS1D) ? 1 : 4;
			readOneSiteHoppings(oneSiteHoppings,geometryParams_.filename,dirs);
			if (oneSiteHoppings.size()!=dirs*edof*edof) {
				throw std::runtime_error("Wrong number of hoppings\n");
			}
			size_t sites = geometryParams_.sites;
			for (size_t i=0;i<edof*edof;i++) { // orbital*orbital cases: aa ab ba bb etc
				MatrixType oneT(sites,sites);
				setGeometryFeAs(oneT,i,oneSiteHoppings);
				reorderLadderX(oneT,geometryParams_.leg);
				assert(isHermitian(oneT,true));
				t.push_back(oneT);
			}
			resizeAndZeroOut(edof*sites,edof*sites);
			for (size_t orbitalPair=0;orbitalPair<edof*edof;orbitalPair++) {
				size_t orb1 = (orbitalPair % geometryParams_.orbitals);
				size_t orb2 = orbitalPair/geometryParams_.orbitals;
				for (size_t i=0;i<sites;i++)
					for (size_t j=0;j<sites;j++)
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
			size_t sites = geometryParams_.sites;
			resizeAndZeroOut(sites,sites);
			for (size_t i=0;i<sites;i++) {
				for (size_t j=0;j<sites;j++) {
					if (i-j==1 || j-i==1) t_(i,j) = geometryParams_.hopping[0];
				}
			}
			if (geometryParams_.option==GeometryParamsType::OPTION_PERIODIC)
				t_(0,sites-1) = t_(sites-1,0) = geometryParams_.hopping[0];
		}

		void setGeometryLadder()
		{
			assert(geometryParams_.hopping.size()==2);
			size_t leg = geometryParams_.leg;
			if (leg<2)
				throw std::runtime_error("GeometryLibrary:: ladder must have leg>1\n");
			assert(!geometryParams_.isPeriodic[DIRECTION_Y] || leg>2);
			size_t sites = geometryParams_.sites;
			resizeAndZeroOut(sites,sites);
			for (size_t i=0;i<sites;i++) {
				PsimagLite::Vector<size_t>::Type v;
				ladderFindNeighbors(v,i,leg,geometryParams_.isPeriodic);
				for (size_t k=0;k<v.size();k++) {
					size_t j = v[k];
					RealType tmp = (ladderSameColumn(i,j,leg))? 
					    geometryParams_.hopping[1] : geometryParams_.hopping[0];
					t_(i,j) = t_(j,i) = tmp; 
				}
			}
		}

		void ladderFindNeighbors(PsimagLite::Vector<size_t>::Type& v,
		                         size_t i,
		                         size_t leg,
		                         const PsimagLite::Vector<bool>::Type& isPeriodic)
		{
			size_t sites = geometryParams_.sites;
			size_t lengthX = static_cast<size_t>(sites/leg);

			v.clear();
			size_t k = i + 1;
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

		bool ladderSameColumn(size_t i,size_t k,size_t leg)
		{
			size_t i1 = i/leg;
			size_t k1 = k/leg;
			if (i1 == k1) return true;
			return false;
		}

		// only 2 orbitals supported
		void setGeometryFeAs(MatrixType& t,
		                     size_t orborb,
		                     const typename PsimagLite::Vector<RealType>::Type& oneSiteHoppings)
		{
			size_t sites = geometryParams_.sites;
			size_t leg = geometryParams_.leg;
			size_t geometryOption = geometryParams_.option;
			size_t orbitalsSquared = geometryParams_.orbitals * geometryParams_.orbitals;
			RealType tx = oneSiteHoppings[orborb+DIRECTION_X*orbitalsSquared];
			size_t lengthx  = sites/leg;
			if (sites%leg!=0) throw std::runtime_error("Leg must divide number of sites.\n");
			for (size_t j=0;j<leg;j++) {
				for (size_t i=0;i<lengthx;i++) {
					if (i+1<lengthx) t(i+1+j*lengthx,i+j*lengthx) = t(i+j*lengthx,i+1+j*lengthx) = tx;
					if (i>0) t(i-1+j*lengthx,i+j*lengthx) = t(i+j*lengthx,i-1+j*lengthx) = tx;
				}
				if (geometryOption==GeometryParamsType::OPTION_PERIODIC) 
					t(j*lengthx,lengthx-1+j*lengthx) = t(lengthx-1+j*lengthx,j*lengthx) = tx;
			}

			if (geometryParams_.type==FEAS1D) return;

			RealType ty = oneSiteHoppings[orborb+DIRECTION_Y*orbitalsSquared];
			for (size_t i=0;i<lengthx;i++) {
				for (size_t j=0;j<leg;j++) {
					if (j>0) t(i+(j-1)*lengthx,i+j*lengthx) = t(i+j*lengthx,i+(j-1)*lengthx) = ty;
					if (j+1<leg) t(i+(j+1)*lengthx,i+j*lengthx) = t(i+j*lengthx,i+(j+1)*lengthx) = ty;
				}
				if (geometryOption==GeometryParamsType::OPTION_PERIODIC)
					t(i,i+(leg-1)*lengthx) = t(i+(leg-1)*lengthx,i) = ty;
			}
			RealType txpy = oneSiteHoppings[orborb+DIRECTION_XPY*orbitalsSquared];
			RealType txmy = oneSiteHoppings[orborb+DIRECTION_XMY*orbitalsSquared];
			for (size_t i=0;i<lengthx;i++) {
				for (size_t j=0;j<leg;j++) {
					if (j+1<leg && i+1<lengthx)
						t(i+1+(j+1)*lengthx,i+j*lengthx) = t(i+j*lengthx,i+1+(j+1)*lengthx) = txpy;
					if (i+1<lengthx && j>0)
						t(i+1+(j-1)*lengthx,i+j*lengthx) = t(i+j*lengthx,i+1+(j-1)*lengthx) = txmy;
					if (geometryOption!=GeometryParamsType::OPTION_PERIODIC || i>0)
						continue;
					if (j+1<leg)
						t((j+1)*lengthx,lengthx-1+j*lengthx) = t(lengthx-1+j*lengthx,(j+1)*lengthx) = txpy;
					if (j>0)
						t((j-1)*lengthx,lengthx-1+j*lengthx) = t(lengthx-1+j*lengthx,(j-1)*lengthx) = txmy;
				}
				if (geometryOption!=GeometryParamsType::OPTION_PERIODIC)
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

		// from 0--1--2--...
		//      N-N+1-N+2--..
		// into:
		//      0--2--4--
		//      1--3--5--
		//
		void reorderLadderX(MatrixType& told,size_t leg)
		{
			size_t sites = geometryParams_.sites;
			MatrixType tnew(told.n_row(),told.n_col());
			for (size_t i=0;i<sites;i++) {
				for (size_t j=0;j<sites;j++) {
					size_t i2 = reorderLadderX(i,leg);
					size_t j2 = reorderLadderX(j,leg);
					tnew(i2,j2) = told(i,j); 
				}
			}
			told = tnew;
		}

		size_t reorderLadderX(size_t i,size_t leg)
		{
			size_t sites = geometryParams_.sites;
			size_t lengthx  = sites/leg;
			// i = x + y*lengthx;
			size_t x = i%lengthx;
			size_t y = i/lengthx;
			return y + x*leg;
		}

		void readOneSiteHoppings(typename PsimagLite::Vector<RealType>::Type& v,
		                         const PsimagLite::String& filename,
		                         size_t dirs)
		{
			typename PsimagLite::IoSimple::In io(filename);
			v.clear();
			for (size_t i=0;i<dirs;i++) {
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

		void appendToVector(typename PsimagLite::Vector<RealType>::Type& v,
		                    const MatrixType& m) const
		{
			for (size_t i=0;i<m.n_row();i++) {
				for (size_t j=0;j<m.n_col();j++) {
					v.push_back(m(i,j));
				}
			}
		}

		template<typename MatrixComplexType>
		void getFourierMatrix(MatrixComplexType& B,size_t leg) const
		{
			typedef typename MatrixComplexType::value_type ComplexType;
			if (geometryParams_.geometryType!=FEAS)
				throw std::runtime_error("getFourierMatrix: unsupported\n");
			size_t n = B.n_row();
			size_t lengthx = n/leg;
			if (n%leg !=0) throw std::runtime_error("Leg must divide number of sites for FEAS\n");
			for (size_t i=0;i<n;i++) {
				size_t rx = i % lengthx;
				size_t ry = i/lengthx;
				for (size_t k=0;k<n;k++) {
					size_t kx = k % lengthx;
					size_t ky = k / lengthx;
					RealType tmpx = rx*kx*2.0*M_PI/lengthx;
					RealType tmpy = ry*ky*2.0*M_PI/leg;
					RealType tmp1 = cos(tmpx+tmpy);
					RealType tmp2 = sin(tmpx+tmpy);
					B(i,k) = ComplexType(tmp1,tmp2);
				}
			}
		}

		void resizeAndZeroOut(size_t nrow,size_t ncol)
		{
			t_.resize(nrow,ncol);
			for (size_t i=0;i<nrow;i++)
				for(size_t j=0;j<ncol;j++)
					t_(i,j)=0;
		}

		void readPotential(typename PsimagLite::Vector<RealType>::Type& v,
		                   const PsimagLite::String& filename)
		{
			typename PsimagLite::Vector<RealType>::Type w;
			PsimagLite::IoSimple::In io(filename);
			try {
				io.read(w,"PotentialT");
			} catch (std::exception& e) {
				std::cerr<<"INFO: No PotentialT in file "<<filename<<"\n";
			}
			io.rewind();

			io.read(v,"potentialV");
			if (w.size()==0) return;
			if (v.size()>w.size()) v.resize(w.size());
			for (size_t i=0;i<w.size();i++) v[i] += w[i];
		}

		const GeometryParamsType& geometryParams_;
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
