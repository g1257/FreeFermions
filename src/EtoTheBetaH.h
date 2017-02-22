/*
Copyright  2009-2017, UT-Battelle, LLC
All rights reserved

[FreeFermions, Version 1.]
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

/*! \file EToTheBetaH.h
 *
 *
 *
 */
#ifndef E_TO_THE_BETA_H_H
#define E_TO_THE_BETA_H_H

namespace FreeFermions {
// All interactions == 0
template<typename EngineType_>
class EToTheBetaH {
public:
	typedef EngineType_ EngineType;
	typedef typename EngineType::RealType RealType;
	typedef typename EngineType::FieldType FieldType;

	EToTheBetaH(RealType beta,
	            const EngineType& engine,
	            RealType energyOffset)
	    : beta_(beta),
	      engine_(engine),
	      energyOffset_(energyOffset)
	{}

	template<typename FreeOperatorsType>
	FieldType operator()(const FreeOperatorsType& freeOps,
	                     SizeType loc) const
	{
		RealType sum = 0;
		for (SizeType i=0;i<loc;i++) {
			// 					if (freeOps[i].type != FreeOperatorsType::CREATION &&
			// 						freeOps[i].type != FreeOperatorsType::DESTRUCTION)
			// 						   continue;
			// 					int sign =  (freeOps[i].type ==
			// 							       FreeOperatorsType::CREATION) ? -1 : 1;
			if (freeOps[i].type != FreeOperatorsType::CREATION) continue;
			sum += engine_.eigenvalue(freeOps[i].lambda); //*sign;
		}

		RealType exponent = -beta_*sum;
		return exp(exponent);
	}


	void transpose() { }

private:

	RealType beta_;
	const EngineType& engine_;
	RealType energyOffset_;
}; // EToTheBetaH
} // namespace Dmrg 

/*@}*/
#endif
