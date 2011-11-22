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

/*! \file FreeFermions.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef R_DENSITY_MATRIX_H
#define R_DENSITY_MATRIX_H
#include <assert.h>
#include <cstdlib>
#include "Engine.h"
#include "GeometryLibrary.h"
#include "CanonicalStates.h"
#include "CreationOrDestructionOp.h"
#include "Range.h"
#include "RealSpaceState.h"
#include "BLAS.h"

namespace FreeFermions {
	template<typename EngineType>
	class ReducedDensityMatrix {
		
		typedef std::vector<size_t> VectorUintType;
		typedef typename EngineType::FieldType FieldType;
		typedef typename EngineType::RealType RealType;
		typedef std::vector<FieldType> VectorType;
		typedef PsimagLite::Matrix<FieldType> MatrixType;
		typedef typename EngineType::ConcurrencyType ConcurrencyType;
		typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
//		typedef FreeFermions::HilbertState<OperatorType> HilbertStateType;
		typedef FreeFermions::RealSpaceState<OperatorType> HilbertStateType;
		typedef typename EngineType::GeometryLibraryType GeometryLibraryType;
		typedef typename OperatorType::FactoryType OpNormalFactoryType;

		enum {CREATION = OperatorType::CREATION,
				       DESTRUCTION = OperatorType::DESTRUCTION};

		public:
			// note: right and left blocks are assumed equal and of size n
			ReducedDensityMatrix(EngineType& engine,size_t n,size_t ne)
			: engine_(engine),concurrency_(engine.concurrency()),n_(n),ne_(ne)
			{
				assert(engine_.dof()==1);
				calculatePsi(psi_);
				//std::cout<<psi_;
				if (!concurrency_.root()) return;
				calculateRdm(rho_,psi_);
			}
			
			size_t rank() const { return rho_.n_row(); }
			
			void diagonalize(std::vector<RealType>& e)
			{
				if (!concurrency_.root()) return;
				diag(rho_,e,'N');
			}

		private:
			void calculatePsi(MatrixType& psi)
			{
				OpNormalFactoryType opNormalFactory(engine_);
				CanonicalStates aux(n_,ne_);
				size_t states = aux.states();
				psi.resize(states,states);
				
				VectorUintType neV(engine_.dof(),ne_);
				VectorUintType zeroV(engine_.dof(),0);
				HilbertStateType gs(engine_,neV);
				
				std::cout<<"#psi of size"<<states<<"x"<<states<<"\n";
				size_t each = states/10;
				VectorType psiV(states);
				std::vector<VectorType> psiVv(states);
				RealType sum = 0;
				PsimagLite::Range<ConcurrencyType> range(0,states,concurrency_);
				for (;!range.end();range.next()) {
					size_t i = range.index();
					if (i%each ==0 && concurrency_.name()=="serial") {
						std::cerr<<"Done "<<(i*10/each)<<"%\n";
						std::cerr.flush();
					}
					VectorUintType v;
					aux.getSites(v,i);
					HilbertStateType phi(engine_,zeroV);
					psiOneBlock(phi,v,CREATION,opNormalFactory);
					//std::cout<<phi;
					for (size_t j=0;j<states;j++) {
						VectorUintType w;
						aux.getSites(w,j);
						if (v.size()+w.size()!=ne_) {
							psiV[j] = 0;
							continue;
						}

						HilbertStateType phi2 = gs;
						psiOneBlock(phi2,w,DESTRUCTION,opNormalFactory,n_);
						psiV[j] = scalarProduct(phi2,phi);
						//std::cout<<psi(i,j)<<" ";
					}
					psiVv[i]=psiV;
					sum += psiV*psiV;
					//std::cout<<"\n";
				}
				concurrency_.gather(psiVv);
				concurrency_.gather(sum);
				std::cerr<<"sum="<<sum<<"\n";
				for (size_t i=0;i<psiVv.size();i++)
					for (size_t j=0;j<psiVv[i].size();j++)
						psi(i,j) = psiVv[i][j]/sqrt(sum);
			}
			
			void psiOneBlock(HilbertStateType& phi,
			                   const VectorUintType& v,
			                   size_t label,
			                   OpNormalFactoryType& opNormalFactory,
			                   size_t offset = 0)
			{
				size_t sigma = 0; 
				for (size_t i=0;i<v.size();i++) {
					size_t site = v[i] + offset;
					OperatorType& myOp = opNormalFactory(label,site,sigma);
					myOp.applyTo(phi);
				}
			}
			
			void calculateRdm(MatrixType& rho,const MatrixType& psi)
			{
				size_t states=psi.n_row();
				rho.resize(states,states);
				//std::cout<<"#rho of size "<<states<<"x"<<states<<"\n";
				RealType alpha = 1.0;
				RealType beta = 0.0;
				psimag::BLAS::GEMM('N','C',states,states,states,alpha,&(psi(0,0)),states,
				     &(psi(0,0)),states,beta,&(rho(0,0)),states);
				if (!isHermitian(rho))
				  throw std::runtime_error("DensityMatrix not Hermitian\n");
			}
			
			EngineType& engine_;
			ConcurrencyType& concurrency_;
			size_t n_; // number of sites for one block only (both blocks are assumed equal)
			size_t ne_; // number of electrons in the combined lattice (right+left)
			MatrixType psi_; // the overlap of the g.s. with the product state
			MatrixType rho_; // the reduced density matrix (reduced over right block)
	}; // ReducedDensityMatrix
} // FreeFermions namespace
/*@}*/
#endif // R_DENSITY_MATRIX_H
