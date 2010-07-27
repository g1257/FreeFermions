// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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
	template<typename RealType_,typename FieldType_,typename UnsignedIntegerType>
	class Engine {
			
			typedef psimag::Matrix<FieldType_> MatrixType;
	
		public:
			typedef RealType_ RealType;
			typedef FieldType_ FieldType;
			typedef HilbertVector<RealType,FieldType,UnsignedIntegerType> HilbertVectorType;
			typedef FreeOperator<RealType,FieldType,UnsignedIntegerType,HilbertVector> FreeOperatorType;
			typedef typename HilbertVectorType::HilbertTermType HilbertTermType;
			
			Engine(const MatrixType& t,size_t dof,bool verbose=false) :
				sites_(t.n_row()),edof_(1),dof_(dof),verbose_(verbose)
			{
				std::vector<MatrixType>* tt = new std::vector<MatrixType>(1);
				(*tt)[0] = t;
				t_ = tt;
				diagonalize();
				if (verbose_) {
					std::cerr<<"#Created core "<<eigenvectors_.n_row();
					std::cerr<<"  times "<<eigenvectors_.n_col()<<"\n";
				}
			}
			
			Engine(const std::vector<MatrixType>& t,size_t dof,bool verbose=false) :
				t_(*t),sites_(t_[0].n_row()),edof_(sqrt(t.size())),dof_(dof),verbose_(verbose)
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
				eigenvectors_.resize(sites_*edof_,sites_*edof_);
				std::vector<MatrixType> eigenvectors(edof_*edof_);
				std::vector<std::vector<RealType> > eigenvalues(edof_*edof_);
				
				for (size_t i=0;i<eigenvectors.size();i++) {
					diagonalize(i,eigenvectors[i],eigenvalues[i]);
				}
				
				size_t counter = 0;
				eigenvalues_.resize(edof_*edof_*sites_);
				for (size_t k=0;k<eigenvalues[0].size();k++) {
					MatrixType ek(edof_,edof_);
					for (size_t gamma=0;gamma<ek.n_row();gamma++) {
						for (size_t gamma2=0;gamma2<ek.n_col();gamma2++) {
							size_t ind = gamma + gamma2*ek.n_row();
							ek(gamma,gamma2) = eigenvalues[ind][k];
						}
					}
					std::vector<RealType> elambda(edof_);
					utils::diag(ek,elambda,'V');
					for (size_t lambda=0;lambda<elambda.size();lambda++) {
						eigenvalues_[counter++] = elambda[lambda];
						for (size_t lambda2=0;lambda2<elambda.size();lambda2++) {
							eigenvecSecondTransform(k,lambda,lambda2,ek,eigenvectors);
						}
					}
				}
				std::vector<size_t> perm(eigenvalues_.size());
				utils::sort(eigenvalues_,perm);
				
				permuteEigenvectors(perm);
				
				if (verbose_) {
					std::cerr<<"Matrix\n";
					std::cerr<<eigenvectors;
					utils::vectorPrint(eigenvalues,"eigenvalues",std::cerr);
					std::cerr<<"*************\n";
					std::cerr<<"Eigenvectors:\n";
					std::cerr<<eigenvectors_;
				}
			}
			
			void diagonalize(size_t orbitalPair,MatrixType& eigenvectors,std::vector<RealType>& eigenvalues)
			{
				eigenvectors = (*t_)[orbitalPair];
				eigenvalues.resize(eigenvectors.n_row());
				utils::diag(eigenvectors,eigenvalues,'V');
			}
			
			void eigenvecSecondTransform(
					size_t k,
					size_t lambda,
					size_t lambda2,
					const MatrixType& ek,
					const std::vector<MatrixType>& eigenvectors)
			{
				for (size_t i=0;i<sites_;i++) {
					eigenvectors_(i+lambda*sites_,k+lambda2*sites_) = 
						eigenvecSecondTransform(lambda,lambda2,i,k,ek,eigenvectors);
				}
			}
			
			FieldType eigenvecSecondTransform(
						size_t lambda,
						size_t lambda2,
						size_t i,
						size_t k,
						const MatrixType& ek,
						const std::vector<MatrixType>& eigenvectors)
			{
				FieldType sum = 0.0;
				for (size_t gamma=0;gamma<ek.n_row();gamma++) {
					for (size_t gamma2=0;gamma2<ek.n_col();gamma2++) {
						sum += std::conj(ek(lambda,gamma))*
							eigenvectors[gamma+gamma2*ek.n_row()](i,k)*
							ek(lambda2,gamma2);
					}
				}
				return sum;
			}
			
			void permuteEigenvectors(const std::vector<size_t>& perm)
			{
				MatrixType a(eigenvectors_.n_row(),eigenvectors_.n_col());
				for (size_t il=0;il<a.n_row();il++) {
					for (size_t il2=0;il2<a.n_col();il2++) {
						a(il,il2)=eigenvectors_(perm[il],perm[il2]);
					}
				}
				eigenvectors_ = a;
			}
			
			const std::vector<MatrixType>* t_;
			size_t sites_;
			size_t edof_;
			size_t dof_;
			bool verbose_;
			psimag::Matrix<FieldType> eigenvectors_;
			std::vector<RealType> eigenvalues_;
	}; // Engine
} // namespace Dmrg 

/*@}*/
#endif
