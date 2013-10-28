/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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
#include "RealSpaceState.h"
#include "BLAS.h"
#include "GeometryParameters.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "InputNg.h"
#include "InputCheck.h"

namespace FreeFermions {
	template<typename EngineType>
	class ReducedDensityMatrix {
		
		typedef typename PsimagLite::Vector<size_t>::Type VectorUintType;
		typedef typename EngineType::FieldType FieldType;
		typedef typename EngineType::RealType RealType;
		typedef typename PsimagLite::Vector<FieldType>::Type VectorType;
		typedef PsimagLite::Matrix<FieldType> MatrixType;
		typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
		typedef FreeFermions::RealSpaceState<OperatorType> HilbertStateType;
		typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
		typedef FreeFermions::GeometryParameters<RealType,InputNgType::Readable> GeometryParamsType;
		typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
		typedef typename OperatorType::FactoryType OpNormalFactoryType;
		typedef PsimagLite::Concurrency ConcurrencyType;

		enum {CREATION = OperatorType::CREATION,
				       DESTRUCTION = OperatorType::DESTRUCTION};

		class MyLoop {

		public:

			MyLoop(EngineType& engine,
			       SizeType n,
			       SizeType ne,
			       CanonicalStates& aux,
			       SizeType states,
			       SizeType nthreads)
			    : engine_(engine),
			      n_(n),
			      ne_(ne),
			      aux_(aux),
			      neV_(engine_.dof(),ne_),
			      zeroV_(engine_.dof(),0),
			      gs_(engine_,neV_),
			      psiV_(states),
			      psiVv_(states),
			      sumV_(ConcurrencyType::storageSize(nthreads),0),
			      opNormalFactory_(engine_)
			{}

			void thread_function_(SizeType threadNum,
			                      SizeType blockSize,
			                      SizeType total,
			                      typename ConcurrencyType::MutexType* myMutex)
			{
				size_t each = total/10;

				SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
				SizeType npthreads = ConcurrencyType::npthreads;
				for (SizeType p=0;p<blockSize;p++) {
					SizeType i = (threadNum+npthreads*mpiRank)*blockSize + p;
					if (i>=total) break;

					if (i%each ==0 && threadNum == 0) {
						std::cerr<<"Done "<<(i*10/each)<<"%\n";
						std::cerr.flush();
					}

					SizeType ind = ConcurrencyType::storageIndex(threadNum);

					VectorUintType v;
					aux_.getSites(v,i);
					HilbertStateType phi(engine_,zeroV_);
					psiOneBlock(phi,v,CREATION,opNormalFactory_);
					//std::cout<<phi;
					for (size_t j=0;j<total;j++) {
						VectorUintType w;
						aux_.getSites(w,j);
						if (v.size()+w.size()!=ne_) {
							psiV_[j] = 0;
							continue;
						}

						HilbertStateType phi2 = gs_;
						psiOneBlock(phi2,w,DESTRUCTION,opNormalFactory_,n_);
						psiV_[j] = scalarProduct(phi2,phi);
						//std::cout<<psi(i,j)<<" ";
					}
					psiVv_[i]=psiV_;
					sumV_[ind] += psiV_*psiV_;
					//std::cout<<"\n";
				}
			}

			FieldType sum() const
			{
				assert(sumV_.size()>0);
				return sumV_[0];
			}

			const typename PsimagLite::Vector<VectorType>::Type& psiVv() const
			{
				return psiVv_;
			}

			void sync()
			{
				SizeType tmp = PsimagLite::sum(sumV_);

				PsimagLite::MPI::pointByPointGather(tmp);
				PsimagLite::MPI::bcast(tmp);

				PsimagLite::MPI::pointByPointGather(psiVv_);
				PsimagLite::MPI::bcast(psiVv_);

				sumV_[0] = tmp;
			}

		private:

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

			EngineType engine_;
			SizeType n_;
			SizeType ne_;
			CanonicalStates aux_;
			VectorUintType neV_;
			VectorUintType zeroV_;
			HilbertStateType gs_;
			VectorType psiV_;
			typename PsimagLite::Vector<VectorType>::Type psiVv_;
			typename PsimagLite::Vector<FieldType>::Type sumV_;
			OpNormalFactoryType opNormalFactory_;
		}; // class MyLoop

		public:

			// note: right and left blocks are assumed equal and of size n
			ReducedDensityMatrix(EngineType& engine,size_t n,size_t ne)
			: engine_(engine),n_(n),ne_(ne)
			{
				assert(engine_.dof()==1);
				calculatePsi(psi_);
				//std::cout<<psi_;
				if (!PsimagLite::Concurrency::root()) return;
				calculateRdm(rho_,psi_);
			}
			
			size_t rank() const { return rho_.n_row(); }
			
			void diagonalize(typename PsimagLite::Vector<RealType>::Type& e)
			{
				if (!PsimagLite::Concurrency::root()) return;
				diag(rho_,e,'N');
			}

		private:
			void calculatePsi(MatrixType& psi)
			{
				CanonicalStates aux(n_,ne_);
				size_t states = aux.states();
				psi.resize(states,states);
				
				std::cout<<"#psi of size"<<states<<"x"<<states<<"\n";

				typedef MyLoop MyLoopType;
				typedef PsimagLite::Parallelizer<MyLoopType> ParallelizerType;
				ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
				                              PsimagLite::MPI::COMM_WORLD);

				MyLoopType myLoop(engine_,n_,ne_,aux,states,threadObject.threads());

				std::cout<<"Using "<<threadObject.name();
				std::cout<<" with "<<threadObject.threads()<<" threads.\n";
				threadObject.loopCreate(states,myLoop);

				myLoop.sync();

				FieldType sum = myLoop.sum();
				const typename PsimagLite::Vector<VectorType>::Type& psiVv = myLoop.psiVv();

//				concurrency_.gather(psiVv);
//				concurrency_.gather(sum);
				std::cerr<<"sum="<<sum<<"\n";
				for (size_t i=0;i<psiVv.size();i++)
					for (size_t j=0;j<psiVv[i].size();j++)
						psi(i,j) = psiVv[i][j]/sqrt(sum);
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
			size_t n_; // number of sites for one block only (both blocks are assumed equal)
			size_t ne_; // number of electrons in the combined lattice (right+left)
			MatrixType psi_; // the overlap of the g.s. with the product state
			MatrixType rho_; // the reduced density matrix (reduced over right block)
	}; // ReducedDensityMatrix
} // FreeFermions namespace
/*@}*/
#endif // R_DENSITY_MATRIX_H
