// BEGIN LICENSE BLOCK
/*
Copyright  2009 , UT-Battelle, LLC
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
// END LICENSE BLOCK
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

namespace FreeFermions {
	template<typename MatrixType>
	class GeometryLibrary {
			typedef typename MatrixType::value_type FieldType;
		public:
			enum {CHAIN,LADDER,FEAS};
			enum {OPTION_NONE,OPTION_PERIODIC};
			enum {DIRECTION_X=0,DIRECTION_Y=1,DIRECTION_XPY=2,DIRECTION_XMY=3};
			
			GeometryLibrary(size_t n,size_t geometryType,FieldType hopping = 1.0) :
				sites_(n),geometryType_(geometryType),hopping_(hopping)
			{
			}
			void setGeometry(MatrixType& t,size_t geometryOption=OPTION_NONE)
			{
				switch (geometryType_) {
					case CHAIN:
						setGeometryChain(t,geometryOption);
						break;
					case LADDER:
						setGeometryLadder(t,geometryOption);
						break;
					default:
						throw std::runtime_error("GeometryLibrary:: Unknown geometry type\n");
				}
			}
			
			void setGeometry(std::vector<MatrixType>& t,const std::string& filename,size_t leg,size_t geometryOption=OPTION_NONE)
			{
				if (geometryType_!=FEAS) throw std::runtime_error("Unsupported\n");
				size_t edof = 2; // 2 orbitals
				std::vector<FieldType> oneSiteHoppings;
				readOneSiteHoppings(oneSiteHoppings,filename);
	
				t.clear();
				for (size_t i=0;i<edof*edof;i++) { // 4 cases: aa ab ba and bb
					MatrixType oneT(sites_,sites_);
					setGeometryFeAs(oneT,i,oneSiteHoppings,leg,geometryOption);
					reorderLadderX(oneT,leg);
					if (!isHermitian(oneT,true)) throw std::runtime_error("Hopping matrix not hermitian\n");
					t.push_back(oneT);
				}
			}

			void bathify(MatrixType& t,const std::vector<FieldType>& tb)
			{
				if (sites_!=t.n_row() || sites_!=t.n_col())
					throw std::runtime_error("GeometryLibrary::bathify(...)\n");
				size_t nb = tb.size();
				size_t nnew = sites_*(1+nb);
				MatrixType tnew(nnew,nnew);
				for (size_t i=0;i<t.n_row();i++)
				{
					for (size_t j=0;j<t.n_col();j++)
						tnew(i,j) = t(i,j);
					for (size_t j=0;j<nb;j++) {
						size_t k = sites_ + j + nb*i;
						tnew(i,k) = tnew(k,i) = tb[j];
					}
				}
				t = tnew;
				sites_ = nnew;
			}
			
			template<typename SomeRealType>
			void addPotential(MatrixType &t,const std::vector<SomeRealType>& p)
			{
				for (size_t i=0;i<p.size();i++) t(i,i) = p[i];
			}	
			
			template<typename ComplexFieldType>
			void fourierTransform(std::vector<ComplexFieldType>& dest,const MatrixType& src,size_t leg) const
			{
				
				size_t n = src.n_row();
				if (n!=sites_) throw std::runtime_error("src must have the same number of sites as lattice\n");
				PsimagLite::Matrix<ComplexFieldType> B(n,n);
				getFourierMatrix(B,leg);
				for (size_t k=0;k<n;k++) {
					ComplexFieldType sum = 0.0;
					for (size_t i=0;i<n;i++) {
						for (size_t j=0;j<n;j++) {
							sum += std::conj(B(i,k)) * src(i,j) * B(j,k);
						}
					}
					dest[k] = sum;
				}
				 
			}
			
			size_t sites() const { return sites_; }
			
		private:
			void setGeometryChain(MatrixType& t,size_t geometryOption)
			{
				for (size_t i=0;i<sites_;i++) {
					for (size_t j=0;j<sites_;j++) {
						if (i-j==1 || j-i==1) t(i,j) = hopping_;
					}
				}
				if (geometryOption==OPTION_PERIODIC) t(0,sites_-1) = t(sites_-1,0) = hopping_;
			}
			
			void setGeometryLadder(MatrixType& t,size_t leg)
			{
				if (leg<2) throw std::runtime_error("GeometryLibrary:: ladder must have leg>1\n");
				for (size_t i=0;i<sites_;i++) {
					std::vector<size_t> v;
					ladderFindNeighbors(v,i,leg);
					for (size_t k=0;k<v.size();k++) {
						size_t j = v[k];
						t(i,j) = t(j,i) = hopping_;
					}
				}
			}
			
			void ladderFindNeighbors(std::vector<size_t>& v,size_t i,size_t leg)
			{
				v.clear();
				size_t k = i + 1;
				if (k<sites_ && ladderSameColumn(k,i,leg)) v.push_back(k);
				k = i + leg;
				if (k<sites_)  v.push_back(k);
				
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
			void setGeometryFeAs(MatrixType& t,size_t orborb,const std::vector<FieldType>& oneSiteHoppings,
						size_t leg,size_t geometryOption)
			{
				FieldType tx = oneSiteHoppings[orborb+DIRECTION_X*4];
				size_t lengthx  = sites_/leg;
				if (sites_%leg!=0) throw std::runtime_error("Leg must divide number of sites.\n");
				for (size_t j=0;j<leg;j++) {
					for (size_t i=0;i<lengthx;i++) {
						if (i+1<lengthx) t(i+1+j*lengthx,i+j*lengthx) = t(i+j*lengthx,i+1+j*lengthx) = tx;
						if (i>0) t(i-1+j*lengthx,i+j*lengthx) = t(i+j*lengthx,i-1+j*lengthx) = tx;
					}
					if (geometryOption==OPTION_PERIODIC) 
						t(j*lengthx,lengthx-1+j*lengthx) = t(lengthx-1+j*lengthx,j*lengthx) = tx;
				}

				FieldType ty = oneSiteHoppings[orborb+DIRECTION_Y*4];
				for (size_t i=0;i<lengthx;i++) {
					for (size_t j=0;j<leg;j++) {
						if (j>0) t(i+(j-1)*lengthx,i+j*lengthx) = t(i+j*lengthx,i+(j-1)*lengthx) = ty;
						if (j+1<leg) t(i+(j+1)*lengthx,i+j*lengthx) = t(i+j*lengthx,i+(j+1)*lengthx) = ty;
					}
					if (geometryOption==OPTION_PERIODIC)  t(i,i+(leg-1)*lengthx) = t(i+(leg-1)*lengthx,i) = ty;
				}
				FieldType txpy = oneSiteHoppings[orborb+DIRECTION_XPY*4];
				FieldType txmy = oneSiteHoppings[orborb+DIRECTION_XMY*4];
				for (size_t i=0;i<lengthx;i++) {
					for (size_t j=0;j<leg;j++) {
						if (j+1<leg && i+1<lengthx)
							t(i+1+(j+1)*lengthx,i+j*lengthx) = t(i+j*lengthx,i+1+(j+1)*lengthx) = txpy;
						if (i+1<lengthx && j>0)
							t(i+1+(j-1)*lengthx,i+j*lengthx) = t(i+j*lengthx,i+1+(j-1)*lengthx) = txmy;
						if (geometryOption!=OPTION_PERIODIC || i>0) continue;
						if (j+1<leg)
							t((j+1)*lengthx,lengthx-1+j*lengthx) = t(lengthx-1+j*lengthx,(j+1)*lengthx) = txpy;
						if (j>0)
							t((j-1)*lengthx,lengthx-1+j*lengthx) = t(lengthx-1+j*lengthx,(j-1)*lengthx) = txmy;
					}
					if (geometryOption!=OPTION_PERIODIC) continue;
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
				MatrixType tnew(told.n_row(),told.n_col());
				for (size_t i=0;i<sites_;i++) {
					for (size_t j=0;j<sites_;j++) {
						size_t i2 = reorderLadderX(i,leg);
						size_t j2 = reorderLadderX(j,leg);
						tnew(i2,j2) = told(i,j); 
					}
				}
				told = tnew;
			}

			size_t reorderLadderX(size_t i,size_t leg)
			{
				size_t lengthx  = sites_/leg;
				// i = x + y*lengthx;
				size_t x = i%lengthx;
				size_t y = i/lengthx;
				return y + x*leg;
			}
			
			void readOneSiteHoppings(std::vector<FieldType>& v,
							const std::string& filename)
			{
				typename PsimagLite::IoSimple::In io(filename);
				io.read(v,"hoppings");
			}
			
			template<typename MatrixComplexType>
			void getFourierMatrix(MatrixComplexType& B,size_t leg) const
			{
				typedef typename MatrixComplexType::value_type ComplexType;
				if (geometryType_!=FEAS) throw std::runtime_error("getFourierMatrix: unsupported\n");
				size_t n = B.n_row();
				size_t lengthx = n/leg;
				if (n%leg !=0) throw std::runtime_error("Leg must divide number of sites for FEAS\n");
				for (size_t i=0;i<n;i++) {
					size_t rx = i % lengthx;
					size_t ry = i/lengthx;
					for (size_t k=0;k<n;k++) {
						size_t kx = k % lengthx;
						size_t ky = k / lengthx;
						FieldType tmpx = rx*kx*2.0*M_PI/lengthx;
						FieldType tmpy = ry*ky*2.0*M_PI/leg;
						FieldType tmp1 = cos(tmpx+tmpy);
						FieldType tmp2 = sin(tmpx+tmpy);
						B(i,k) = ComplexType(tmp1,tmp2);
					}
				}
			}
			
			size_t sites_;
			size_t geometryType_;
			FieldType hopping_;
	}; // GeometryLibrary
} // namespace FreeFermions

/*@}*/
#endif //GEOMETRY_LIB_H
