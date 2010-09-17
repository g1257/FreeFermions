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
#include "Engine.h"
#include "ObservableLibrary.h"
#include "GeometryLibrary.h"
#include "CanonicalStates.h"

namespace FreeFermions {
	template<typename RealType,typename FieldType>
	class ReducedDensityMatrix {
		
		typedef std::vector<size_t> VectorUintType;
		typedef std::vector<FieldType> VectorType;
		typedef std::vector<bool> LevelsType;
		typedef psimag::Matrix<FieldType> MatrixType;
		typedef FreeFermions::Engine<RealType,FieldType,LevelsType> EngineType;
		typedef typename EngineType::HilbertVectorType HilbertVectorType;
		typedef typename EngineType::FreeOperatorType FreeOperatorType;
		typedef GeometryLibrary<MatrixType> GeometryLibraryType;
		typedef ObservableLibrary<EngineType> ObservableLibraryType;
		
		public:
			// note: right and left blocks are assumed equal and of size n
			ReducedDensityMatrix(size_t n,size_t ne,size_t dof=1)
			: n_(n),ne_(ne),hoppings_(n,n),geometry_(n,GeometryLibraryType::CHAIN),
				engine_(geometry_,GeometryLibraryType::OPTION_NONE,dof,true)
			{
				calculatePsi(psi_);
				calculateRdm(rho_,psi_);
			}
			
			size_t rank() const { return rho_.n_row(); }
			
			void diagonalize(std::vector<RealType>& e)
			{
				utils::diag(rho_,e,'N');
			}

		private:
			void calculatePsi(MatrixType& psi)
			{
				CanonicalStates aux(n_,ne_);
				size_t states = aux.states();
				psi.resize(states,states);
				
				VectorUintType neV(engine_.dof(),ne_);
				HilbertVectorType gs =  engine_.newGroundState(neV);
				
				for (size_t i=0;i<states;i++) {
					VectorUintType v;
					aux.getSites(v,i);
					HilbertVectorType phi =   engine_.newState();
					psiOneBlock(phi,v);
					for (size_t j=0;j<states;j++) {
						VectorUintType w;
						aux.getSites(w,j);
						if (v.size()+w.size()!=ne_) {
							psi(i,j) = 0;
							continue;
						}
						
						HilbertVectorType phi2 = phi;
						psiOneBlock(phi2,w);
						psi(i,j) = scalarProduct(phi2,gs);
					}
				}
			}
			
			void psiOneBlock(HilbertVectorType& phi,const VectorUintType& v)
			{
				size_t sigma = 0; 
				assert(engine_.dof()==1);
				for (size_t i=0;i<v.size();i++) {
					size_t site = v[i];
					FreeOperatorType myOp = engine_.newSimpleOperator("creation",site,sigma);
					HilbertVectorType phi2 = engine_.newState();
					myOp.apply(phi2,phi,FreeOperatorType::SIMPLIFY);
					phi = phi2;
				}
			}
			
			void calculateRdm(MatrixType& rho,const MatrixType& psi)
			{
				size_t states=psi.n_row();
				rho.resize(states,states);
				for (size_t i1=0;i1<states;i1++) {
					for (size_t i2=0;i2<states;i2++) {
						rho(i1,i2) = 0;
						for (size_t j=0;j<states;j++) {
							rho(i1,i2) += psi(i1,j) * std::conj(psi(i2,j));
						}
					}
				}
			}
			
			size_t n_; // number of sites for one block only (both blocks are assumed equal)
			size_t ne_; // number of electrons in the combined lattice (right+left)
			MatrixType hoppings_;
			GeometryLibraryType geometry_;
			EngineType engine_;
			MatrixType psi_; // the overlap of the g.s. with the product state
			MatrixType rho_; // the reduced density matrix (reduced over right block)
	}; // ReducedDensityMatrix
} // FreeFermions namespace
/*@}*/
#endif // R_DENSITY_MATRIX_H
