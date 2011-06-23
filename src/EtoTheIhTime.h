// BEGIN LICENSE BLOCK
/*
Copyright  2009 , UT-Battelle, LLC
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

/*! \file EToTheIhTime.h
 *
 * 
 *
 */
#ifndef E_TO_THE_I_H_TIME_H
#define E_TO_THE_I_H_TIME_H


namespace FreeFermions {
	// All interactions == 0
	template<typename EngineType>
	class EToTheIhTime {

	public:
			typedef typename EngineType::RealType RealType;
			typedef typename EngineType::FieldType FieldType;

			EToTheIhTime(RealType time,const EngineType& engine,RealType energyOffset) :
				time_(time),engine_(engine),energyOffset_(energyOffset)
			{
			}
			
			template<typename FreeOperatorsType>
			FieldType operator()(const FreeOperatorsType& freeOps,
			                      size_t loc) const
			{
				RealType sum = 0;
				for (size_t i=0;i<loc;i++) {
					if (freeOps[i].type != FreeOperatorsType::CREATION &&
						freeOps[i].type != FreeOperatorsType::DESTRUCTION)
						   continue;
					int sign =  (freeOps[i].type ==
							       FreeOperatorsType::CREATION) ? -1 : 1;
					sum += engine_.eigenvalue(freeOps[i].lambda)*sign;
				}
				//if (fabs(sum)>1e-6) std::cerr<<"ffffffffffffffffff\n";

				RealType exponent = -time_*sum;
				return FieldType(cos(exponent),sin(exponent));
			}

			void transpose() { time_ = -time_; }

		private:
			
			RealType time_;
			const EngineType& engine_;
			RealType energyOffset_;
	}; // EToTheIhTime
} // namespace Dmrg 

/*@}*/
#endif
