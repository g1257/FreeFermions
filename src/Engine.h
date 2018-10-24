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

/*! \file Engine.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef ENGINE_H
#define ENGINE_H
#include "Matrix.h"
#include "Vector.h"
#include <fstream>

namespace FreeFermions {
// All interactions == 0
template<typename RealType_,typename FieldType_>
class Engine {

public:

	typedef RealType_ RealType;
	typedef FieldType_ FieldType;
	typedef  PsimagLite::Matrix<FieldType> MatrixType;

	enum VerboseEnum {VERBOSE_NO, VERBOSE_YES};

	Engine(const PsimagLite::Matrix<FieldType>& geometry,
	       PsimagLite::String outputFile,
	       SizeType dof,
	       VerboseEnum verbose)
	    : dof_(dof),
	      verbose_(verbose),
	      eigenvectors_(geometry),
	      fout_(outputFile.c_str(), std::ios::app)
	{
		diagonalize();
		if (verbose_) {
			fout_<<"Eigenvalues\n";
			fout_<<eigenvalues_;
			fout_<<"#Created core "<<eigenvectors_.n_row();
			fout_<<"  times "<<eigenvectors_.n_col()<<"\n";
		}
	}

	RealType energy(SizeType ne) const
	{
		RealType sum = 0;
		for (SizeType i=0;i<ne;i++) sum += eigenvalues_[i];
		return sum;
	}

	const RealType& eigenvalue(SizeType i) const { return eigenvalues_[i]; }

	const FieldType& eigenvector(SizeType i,SizeType j) const
	{
		return eigenvectors_(i,j);
	}

	void transform(MatrixType& m) const
	{
		MatrixType tmp = multiplyTransposeConjugate(eigenvectors_, m);
		m = tmp * eigenvectors_;
	}

	SizeType dof() const { return dof_; }

	SizeType size() const { return eigenvalues_.size(); }

private:

	void diagonalize()
	{
		if (!isHermitian(eigenvectors_,true)) throw std::runtime_error("Matrix not hermitian\n");

		diag(eigenvectors_,eigenvalues_,'V');

		if (verbose_) {
			fout_<<"eigenvalues\n";
			fout_<<eigenvalues_;
			fout_<<"*************\n";
			fout_<<"Eigenvectors:\n";
			fout_<<eigenvectors_;
		}
	}

	// degrees of freedom that are simply repetition (hoppings are diagonal in these)
	SizeType dof_;
	bool verbose_;
	MatrixType eigenvectors_;
	typename PsimagLite::Vector<RealType>::Type eigenvalues_;
	std::ofstream fout_;
}; // Engine
} // namespace FreeFermions

/*@}*/
#endif
