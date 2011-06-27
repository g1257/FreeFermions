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

/*! \file CreationOrDestructionOp.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef CREATION_OR_DEST_H_H
#define CREATION_OR_DEST_H_H

#include "Complex.h" // in PsimagLite
#include "OperatorFactory.h"

namespace FreeFermions {
	
	template<typename EngineType_>
	class CreationOrDestructionOp {
		typedef CreationOrDestructionOp<EngineType_> ThisType;
	public:
		typedef EngineType_ EngineType;
		typedef typename EngineType::RealType RealType;
		typedef typename EngineType::FieldType FieldType;
		typedef OperatorFactory<ThisType> FactoryType;

		enum {CREATION,DESTRUCTION};

		friend class OperatorFactory<ThisType>;

		CreationOrDestructionOp(const EngineType& engine,
		                           size_t type,
		                           size_t ind,
		                           size_t sigma)
		: engine_(engine),type_(type),ind_(ind),sigma_(sigma)
		{}

		size_t type() const { return type_; }

		size_t sigma() const { return sigma_; }

		FieldType const operator()(size_t j) const
		{
			if (type_==CREATION) return engine_.eigenvector(ind_,j);
			return std::conj(engine_.eigenvector(ind_,j));
		}

		template<typename SomeStateType>
		void applyTo(SomeStateType& state)
		{
			state.pushInto(*this);
		}

		void transpose()
		{
			if (type_==CREATION) {
				type_ = DESTRUCTION;
				return;
			}
			type_ = CREATION;
		}

	private:

		//! Use the factory to create objects of this type
		CreationOrDestructionOp(const ThisType* op2)
		: engine_(op2->engine_),
		  type_(op2->type_),
		  ind_(op2->ind_),
		  sigma_(op2->sigma_) {}

		const EngineType& engine_;
		size_t type_,ind_,sigma_;
	}; // CreationOrDestructionOp
} // namespace Dmrg 

/*@}*/
#endif // CREATION_OR_DEST_H_H
