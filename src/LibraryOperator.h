// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
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

/*! \file LibraryOperator.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef LIBRARY_OPERATOR_H
#define LIBRARY_OPERATOR_H

#include "Complex.h" // in PsimagLite
#include "OperatorFactory.h"

namespace FreeFermions {
	
	template<typename OperatorType>
	class LibraryOperator {
		typedef LibraryOperator<OperatorType> ThisType;
		typedef OperatorFactory<OperatorType> OpNormalFactoryType;
	public:
		typedef typename OperatorType::EngineType EngineType;
		typedef typename OperatorType::RealType RealType;
		typedef typename OperatorType::FieldType FieldType;
		typedef OperatorFactory<ThisType> FactoryType;

		enum {CREATION = OperatorType::CREATION,
		      DESTRUCTION = OperatorType::DESTRUCTION,
		      N,
		      NBAR};

		friend class OperatorFactory<ThisType>;

		template<typename SomeStateType>
		void applyTo(SomeStateType& state)
		{
			OperatorType& op = opNormalFactory_(DESTRUCTION,ind_,sigma_);
			OperatorType& op2 = opNormalFactory_(CREATION,ind_,sigma_);

			if (type_==N) {
				state.pushInto(op);
				state.pushInto(op2);
			} else if (type_==NBAR){
				state.pushInto(op2);
				state.pushInto(op);
			}

		}

	private:

		//! Use OperatorFactory to create objects of this class
		LibraryOperator(const EngineType& engine,
				size_t type,
				size_t ind,
				size_t sigma)
		: opNormalFactory_(engine),type_(type),ind_(ind),sigma_(sigma)
		{}

		LibraryOperator(const ThisType& x)
		{
			throw std::runtime_error(
			  "LibraryOperator::copyCtor: Don't even think of coming here\n");
		}

		ThisType& operator=(const ThisType& x)
		{
			throw std::runtime_error(
			  "LibraryOperator::assignmentOp: Don't even think of coming here\n");
		}

		OpNormalFactoryType opNormalFactory_;
		size_t type_,ind_,sigma_;
	}; // LibraryOperator
	

} // namespace Dmrg 

/*@}*/
#endif // LIBRARY_OPERATOR_H
