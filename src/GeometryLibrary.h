// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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


namespace FreeFermions {
	template<typename MatrixType>
	class GeometryLibrary {
			typedef typename MatrixType::value_type FieldType;
		public:
			enum {CHAIN,LADDER};
			enum {OPTION_NONE,OPTION_PERIODIC};
			
			GeometryLibrary(size_t n,FieldType hopping = 1.0) :
				sites_(n),hopping_(hopping)
			{
			}
			void setGeometry(MatrixType& t,size_t geometryType,size_t geometryOption=OPTION_NONE)
			{
				switch (geometryType) {
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

			size_t sites_;
			FieldType hopping_;
	}; // GeometryLibrary
} // namespace FreeFermions

/*@}*/
#endif //GEOMETRY_LIB_H
