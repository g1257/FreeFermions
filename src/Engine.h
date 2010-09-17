// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
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

/*! \file Engine.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef ENGINE_H
#define ENGINE_H

#include "Utils.h"
#include "HilbertVector.h"
#include "FreeOperator.h"

namespace FreeFermions {
	// All interactions == 0
	template<typename RealType_,typename FieldType_,typename LevelsType>
	class Engine {
			
			typedef psimag::Matrix<FieldType_> MatrixType;
	
		public:
			typedef RealType_ RealType;
			typedef FieldType_ FieldType;
			typedef HilbertVector<RealType,FieldType,LevelsType> HilbertVectorType;
			typedef FreeOperator<HilbertVectorType> FreeOperatorType;
			typedef typename HilbertVectorType::HilbertTermType HilbertTermType;

			template<typename SomeGeometryType>
					Engine(SomeGeometryType& g,size_t geometryOption,size_t dof,bool verbose=false) 
			: sites_(g.sites()),edof_(1),dof_(dof),verbose_(verbose),eigenvectors_(sites_,sites_),
				eigenvalues_(sites_)
			{
				std::vector<MatrixType>* tt = new std::vector<MatrixType>(1);
				(*tt)[0].resize(sites_,sites_);
				g.setGeometry((*tt)[0],geometryOption);
				t_ = tt;
				diagonalize();
				//MatrixType tmp = eigenvectors_;
				//eigenvectors_ = transposeConjugate(tmp);
				if (verbose_) {
					std::cerr<<"#Created core "<<eigenvectors_.n_row();
					std::cerr<<"  times "<<eigenvectors_.n_col()<<"\n";
				}
			}
			
			Engine(const MatrixType& t,size_t dof,bool verbose=false) :
				sites_(t.n_row()),edof_(1),dof_(dof),verbose_(verbose),eigenvectors_(sites_,sites_),
				       eigenvalues_(sites_)
			{
				std::vector<MatrixType>* tt = new std::vector<MatrixType>(1);
				(*tt)[0] = t;
				t_ = tt;
				diagonalize();
				//MatrixType tmp = eigenvectors_;
				//eigenvectors_ = transposeConjugate(tmp);
				if (verbose_) {
					std::cerr<<"#Created core "<<eigenvectors_.n_row();
					std::cerr<<"  times "<<eigenvectors_.n_col()<<"\n";
				}
			}

			Engine(const std::vector<MatrixType>& t,size_t dof,bool verbose=false) :
				t_(&t),sites_(t[0].n_row()),edof_(size_t(sqrt(t.size()))),dof_(dof),verbose_(verbose),
				    eigenvectors_(sites_*edof_,sites_*edof_),eigenvalues_(sites_*edof_)
			{
				if ((edof_) * (edof_) != t.size()) 
					throw std::runtime_error("t.size must be a perfect square\n");
				diagonalize();
				if (verbose_) {
					std::cerr<<"#Created core "<<eigenvectors_.n_row();
					std::cerr<<"  times "<<eigenvectors_.n_col()<<"\n";
				}
			}

			HilbertVectorType newState(bool verbose=false) const
			{
				HilbertVectorType tmp(sites_*edof_,dof_,verbose);
				return tmp;	
			}

			HilbertVectorType newGroundState(const std::vector<size_t>& ne) const
			{
				HilbertVectorType tmp(sites_*edof_,dof_);
				tmp.fill(ne);
				return tmp;
			}

			FreeOperatorType newSimpleOperator(const std::string& label,size_t site,size_t flavor) const
			{
				FreeOperatorType tmp(eigenvectors_,label,site,flavor,dof_);
				return tmp;
			}

			RealType energy(const HilbertTermType& term) const
			{
				return HilbertVectorType::energy(eigenvalues_,term);
			}

			RealType eigenvalue(size_t i) const { return eigenvalues_[i]; }

			size_t dof() const { return dof_; }

		private:
		
			void diagonalize()
			{
				for (size_t i=0;i<edof_*edof_;i++) 
					sumHoppings(i);

				if (!isHermitian(eigenvectors_,true)) throw std::runtime_error("Matrix not hermitian\n");
				
				utils::diag(eigenvectors_,eigenvalues_,'V');
					
				if (verbose_) {
					utils::vectorPrint(eigenvalues_,"eigenvalues",std::cerr);
					std::cerr<<"*************\n";
					std::cerr<<"Eigenvectors:\n";
					std::cerr<<eigenvectors_;
				}
			}
			
			void sumHoppings(size_t orbitalPair)
			{
				size_t orb1 = (orbitalPair & 1);
				size_t orb2 = orbitalPair/2;
				for (size_t i=0;i<sites_;i++)
					for (size_t j=0;j<sites_;j++)
						eigenvectors_(i+orb1*sites_,j+orb2*sites_) += (*t_)[orbitalPair](i,j);
			}

			const std::vector<MatrixType>* t_;
			size_t sites_;
			size_t edof_; // degrees of freedom that are coupled (feas : = 2 , orbitals)
			size_t dof_; // degrees of freedom that are simply repetition (hoppings are diagonal in these)
					// feas =2 spin
			bool verbose_;
			psimag::Matrix<FieldType> eigenvectors_;
			std::vector<RealType> eigenvalues_;
	}; // Engine
} // namespace FreeFermions 

/*@}*/
#endif
