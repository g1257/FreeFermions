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

/*! \file ParallelDecay.h
 *
 * DOC TBW FIXME
 */
#ifndef PARALLEL_DECAY_H
#define PARALLEL_DECAY_H

#include "Matrix.h"
#include "HilbertState.h"
#include "EtoTheIhTime.h"
#include "DiagonalOperator.h"
#include "CreationOrDestructionOp.h"
#include "LibraryOperator.h"
#include "Concurrency.h"

namespace FreeFermions {

template<typename RealType>
struct DecayParams {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	DecayParams(SizeType orbitals_,
	            SizeType numberOfSites_,
	            const VectorSizeType& sites_,
	            SizeType sigma3_,
	            const RealType& offset_,
	            const RealType& step_,
	            bool verbose_)
	    : orbitals(orbitals_),
	      numberOfSites(numberOfSites_),
	      sites(sites_),
	      sigma3(sigma3_),
	      offset(offset_),
	      step(step_),
	      verbose(verbose_)
	{}

	SizeType orbitals;
	SizeType numberOfSites;
	VectorSizeType sites;
	SizeType sigma3;
	RealType offset;
	RealType step;
	bool verbose;
};

template<typename RealType_,typename FieldType,typename EngineType>
class ParallelDecay {

	typedef FreeFermions::EToTheIhTime<EngineType> EtoTheIhTimeType;
	typedef FreeFermions::DiagonalOperator<EtoTheIhTimeType> DiagonalOperatorType;
	typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
	typedef FreeFermions::LibraryOperator<OperatorType> LibraryOperatorType;
	typedef typename DiagonalOperatorType::FactoryType OpDiagonalFactoryType;
	typedef typename OperatorType::FactoryType OpNormalFactoryType;
	typedef typename LibraryOperatorType::FactoryType OpLibFactoryType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	typedef FreeFermions::HilbertState<OperatorType,DiagonalOperatorType> HilbertStateType;
	typedef RealType_ RealType;
	typedef DecayParams<RealType> DecayParamsType;

	ParallelDecay(const EngineType& engine,
	              const DecayParamsType& params,
	              const HilbertStateType& gs)
	    : engine_(engine),params_(params),gs_(gs)
	{
		if (params.sites.size()==0) {
			throw std::runtime_error("ParallelDecay\n");
		}
	}

	void thread_function_(SizeType threadNum,
	                      SizeType blockSize,
	                      SizeType total,
	                      typename ConcurrencyType::MutexType* myMutex)
	{
		SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
		SizeType npthreads = ConcurrencyType::npthreads;

		SizeType siteToMeasure = params_.sites[params_.sites.size()-1];

		for (SizeType p=0;p<blockSize;p++) {
			SizeType it = (threadNum+npthreads*mpiRank)*blockSize + p;
			if (it>=total) continue;

			OpNormalFactoryType opNormalFactory(engine_);
			OpDiagonalFactoryType opDiagonalFactory(engine_);

			RealType time = it * params_.step + params_.offset;
			EtoTheIhTimeType eih(time,engine_,0);
			DiagonalOperatorType& eihOp = opDiagonalFactory(eih);

			FieldType sum = 0;

			for (SizeType sigma1 = 0; sigma1 < 2; ++sigma1) {
				for (SizeType sigma2 = 0; sigma2 < 2; ++sigma2) {
					for (SizeType orb1 = 0; orb1 < params_.orbitals; ++orb1) {
						if (orb1 == 1) continue;
						for (SizeType orb2 = 0; orb2 < params_.orbitals; ++orb2) {
							if (orb2 == 1) continue;

							HilbertStateType phi1 = gs_;
							partialDecay(phi1,
							             orb1,
							             sigma1,
							             opNormalFactory);
							eihOp.applyTo(phi1);

							OperatorType& myOp6 = opNormalFactory(OperatorType::DESTRUCTION,
							                                      siteToMeasure,
							                                      params_.sigma3);
							myOp6.applyTo(phi1);

							HilbertStateType phi2 = gs_;
							partialDecay(phi2,
							             orb2,
							             sigma2,
							             opNormalFactory);
							eihOp.applyTo(phi2);

							myOp6.applyTo(phi2);
							sum += scalarProduct(phi1,phi2);
						}
					}
				}
			}
			std::cout<<time<<" "<<real(sum)<<"\n";
		}
	}

	FieldType calcSuperDensity()
	{
		FieldType sum = 0.0;
		OpNormalFactoryType opNormalFactory(engine_);

		for (SizeType sigma1 = 0; sigma1 < 2; ++sigma1) {
			for (SizeType sigma2 = 0; sigma2 < 2; ++sigma2) {
				for (SizeType orb1 = 0; orb1 < params_.orbitals; ++orb1) {
					if (orb1 == 1) continue;
					for (SizeType orb2 = 0; orb2 < params_.orbitals; ++orb2) {
						if (orb2 == 1) continue;
						HilbertStateType phi1 = gs_;
						partialDecay(phi1,orb1,sigma1,opNormalFactory);
						HilbertStateType phi2 = gs_;
						partialDecay(phi2,orb2,sigma2,opNormalFactory);
						sum += scalarProduct(phi1,phi2);
					}
				}
			}
		}

		return sum;
	}

private:

	void partialDecay(HilbertStateType& phi,
	                  SizeType orb,
	                  SizeType sigma,
	                  OpNormalFactoryType& opNormalFactory)
	{
		assert(params_.sites.size()>0);
		for (SizeType i = 0; i < params_.sites.size() - 1; i++) {
			SizeType site1 = params_.sites[i] + 1*params_.numberOfSites;
			OperatorType& myOp1 = opNormalFactory(OperatorType::DESTRUCTION,site1,sigma);
			myOp1.applyTo(phi);
			SizeType site2 = params_.sites[i] + orb*params_.numberOfSites;
			OperatorType& myOp2 = opNormalFactory(OperatorType::CREATION,site2,sigma);
			myOp2.applyTo(phi);
		}
	}

	const EngineType& engine_;
	const DecayParamsType& params_;
	const HilbertStateType& gs_;
}; // class ParallelDecay
} // namespace FreeFermions

/*@}*/
#endif // PARALLEL_DECAY_H

