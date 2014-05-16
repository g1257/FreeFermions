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
/** \ingroup FreeFermions */
/*@{*/

/*! \file ParallelHolonDoublon.h
 *
 * DOC TBW FIXME
 */
#ifndef PARALLEL_HOLON_DOUBLON
#define PARALLEL_HOLON_DOUBLON

#include "Matrix.h"
#include "HilbertState.h"
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"
#include "CreationOrDestructionOp.h"
#include "LibraryOperator.h"
#include "Concurrency.h"

namespace FreeFermions {

template<typename RealType>
struct HolonDoublonParams {
	HolonDoublonParams(const typename PsimagLite::Vector<SizeType>::Type& ne_,
					   const typename PsimagLite::Vector<SizeType>::Type& sites_,
					   SizeType sigma3_,
					   const RealType& offset_,
					   const RealType& step_,
					   bool debug_,
					   bool verbose_)
		: ne(ne_),
		  sites(sites_),
		  sigma3(sigma3_),
		  offset(offset_),
		  step(step_),
		  debug(debug_),
		  verbose(verbose_)
	{}
	typename PsimagLite::Vector<SizeType>::Type ne;
	typename PsimagLite::Vector<SizeType>::Type sites;
	SizeType sigma3;
	RealType offset;
	RealType step;
	bool debug;
	bool verbose;
};

template<typename RealType_,typename FieldType,typename EngineType>
class ParallelHolonDoublon {

	typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
	typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;
	typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
	typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
	typedef FreeFermions::LibraryOperator<OperatorType> LibraryOperatorType;
	typedef typename DiagonalOperatorType::FactoryType OpDiagonalFactoryType;
	typedef typename OperatorType::FactoryType OpNormalFactoryType;
	typedef typename LibraryOperatorType::FactoryType OpLibFactoryType;
	typedef PsimagLite::Concurrency ConcurrencyType;

public:

	typedef RealType_ RealType;
	typedef HolonDoublonParams<RealType> HolonDoublonParamsType;

	ParallelHolonDoublon(const EngineType& engine,const HolonDoublonParamsType& params)
		: engine_(engine),params_(params),gs_(engine,params.ne,params.debug)
	{
		if (params.sites.size()!=3) {
			throw std::runtime_error("ParallelHolonDoublon: expecting 3 sites\n");
		}
	}

	void thread_function_(SizeType threadNum,
	                      SizeType blockSize,
	                      SizeType total,
	                      typename ConcurrencyType::MutexType* myMutex)
	{
		SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
		SizeType npthreads = ConcurrencyType::npthreads;
		for (SizeType p=0;p<blockSize;p++) {
			SizeType it = (threadNum+npthreads*mpiRank)*blockSize + p;
			if (it>=total) continue;

			OpNormalFactoryType opNormalFactory(engine_);
			OpLibFactoryType opLibFactory(engine_);
			OpDiagonalFactoryType opDiagonalFactory(engine_);

			RealType time = it * params_.step + params_.offset;
			EtoTheIhTimeType eih(time,engine_,0);
			DiagonalOperatorType& eihOp = opDiagonalFactory(eih);

			HilbertStateType savedVector = gs_;
			FieldType savedValue = 0;
			FieldType sum = 0;

			for (SizeType sigma = 0;sigma<2;sigma++) {
				HilbertStateType phi = gs_;
				LibraryOperatorType& myOp = opLibFactory(
											  LibraryOperatorType::N,params_.sites[0],1-sigma);
				myOp.applyTo(phi);

				OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,
												 params_.sites[0],sigma);
				myOp2.applyTo(phi);

				for (SizeType sigma2 = 0;sigma2 < 2;sigma2++) {
					HilbertStateType phi3 = phi;


					LibraryOperatorType& myOp3 = opLibFactory(
									 LibraryOperatorType::NBAR,params_.sites[1],1-sigma2);
					myOp3.applyTo(phi3);

					OperatorType& myOp4 = opNormalFactory(
									 OperatorType::DESTRUCTION,params_.sites[1],sigma2);
					myOp4.applyTo(phi3);

					if (params_.verbose) std::cerr<<"Applying exp(iHt)\n";
					eihOp.applyTo(phi3);

					if (params_.verbose) std::cerr<<"Applying c_p\n";
					OperatorType& myOp6 = opNormalFactory(
									  OperatorType::DESTRUCTION,params_.sites[2],params_.sigma3);
					myOp6.applyTo(phi3);

					if (params_.verbose) std::cerr<<"Adding "<<sigma<<" "<<sigma2<<" "<<it<<"\n";

					if (sigma ==0 && sigma2 ==0) savedVector = phi3;
					if (sigma ==1 && sigma2 ==1) {
						savedValue = scalarProduct(phi3,savedVector);
					}
					sum += scalarProduct(phi3,phi3);
					if (params_.verbose) std::cerr<<"Done with scalar product\n";
				}
			}
			sum += 2*real(savedValue);
			std::cout<<time<<" "<<real(sum)<<"\n";
		}
	}

	FieldType calcSuperDensity(SizeType site,
								 SizeType site2)
	{
		HilbertStateType savedVector = gs_;
		FieldType savedValue = 0;
		FieldType sum = 0;
		OpNormalFactoryType opNormalFactory(engine_);
		OpLibFactoryType opLibFactory(engine_);

		for (SizeType sigma = 0;sigma<2;sigma++) {
			HilbertStateType phi = gs_;

			LibraryOperatorType& myOp = opLibFactory(
									 LibraryOperatorType::N,site,1-sigma);

			myOp.applyTo(phi);
			OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,
									 site,sigma);

			myOp2.applyTo(phi);

			for (SizeType sigma2 = 0;sigma2 < 2;sigma2++) {
				HilbertStateType phi3 = phi;
				LibraryOperatorType& myOp3 = opLibFactory(
										 LibraryOperatorType::NBAR,site2,1-sigma2);
				myOp3.applyTo(phi3);

				OperatorType& myOp4 = opNormalFactory(
										 OperatorType::DESTRUCTION,site2,sigma2);

				myOp4.applyTo(phi3);

				if (sigma ==0 && sigma2 ==0) savedVector = phi3;
				sum += scalarProduct(phi3,phi3);

				if (sigma ==1 && sigma2 ==1) {
					savedValue = scalarProduct(phi3,savedVector);
				}
			}
		}
		sum += 2*real(savedValue);
		//std::cerr<<"#sum2="<<scalarProduct(tmpV,tmpV)<<"\n";
		return sum;
	}

private:

	const EngineType& engine_;
	const HolonDoublonParamsType& params_;
	HilbertStateType gs_;
}; // class ParallelHolonDoublon
} // namespace FreeFermions

/*@}*/
#endif // PARALLEL_HOLON_DOUBLON
