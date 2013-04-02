/*
Copyright (c) 2011-2012, UT-Battelle, LLC
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

/*! \file DriverHelper.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef DRIVER_HELPER_H
#define DRIVER_HELPER_H

#include <cassert>

namespace FreeFermions {


template<typename GeometryLibraryType>
class DriverHelper {

	typedef typename GeometryLibraryType::GeometryParamsType GeometryParamsType;
	typedef typename GeometryParamsType::RealType RealType;

public:

	static void usage(const std::string& thisFile,const std::string& usage)
	{
		std::cout<<thisFile<<": USAGE IS "<<thisFile<<" ";
		std::cout<<usage<<"\n";
		//std::cout<<"\n";
	}

	static int readLabel(const std::string& filename,const std::string& label)
	{
		PsimagLite::IoSimple::In io(filename);
		int x = 0;
		io.readline(x,label);

		return x;
	}

	static void readPotential(std::vector<RealType>& v,const std::string& filename)
	{
		std::vector<RealType> w;
		PsimagLite::IoSimple::In io(filename);
		try {
			io.read(w,"PotentialT");
		} catch (std::exception& e) {
			std::cerr<<"INFO: No PotentialT in file "<<filename<<"\n";
		}
		io.rewind();

		io.read(v,"potentialV");
		if (w.size()==0) return;
		if (v.size()>w.size()) v.resize(w.size());
		for (size_t i=0;i<w.size();i++) v[i] += w[i];
	}

	static void setMyGeometry(GeometryParamsType& geometryParams,const std::vector<std::string>& vstr)
	{
		// default value
		geometryParams.type = GeometryLibraryType::CHAIN;

		if (vstr.size()<2) {
			// assume chain
			return;
		}

		std::string gName = vstr[0];
		if (gName == "chain") {
			throw std::runtime_error("setMyGeometry: -g chain takes no further arguments\n");
		}

		geometryParams.leg = atoi(vstr[1].c_str());

		if (gName == "ladder") {
			if (vstr.size()!=3) {
				usage("setMyGeometry","-g ladder,leg,isPeriodic");
				throw std::runtime_error("setMyGeometry: usage is: -g ladder,leg,isPeriodic \n");
			}
			geometryParams.type = GeometryLibraryType::LADDER;
			geometryParams.hopping.resize(2);
			geometryParams.hopping[0] =  geometryParams.hopping[1]  = 1.0;
			geometryParams.isPeriodicY = (atoi(vstr[2].c_str())>0);
			return;
		}

		if (vstr.size()<3) {
			usage("setMyGeometry"," -g {feas | ktwoniffour} leg filename");
			throw std::runtime_error("setMyGeometry: usage is: -g {feas | ktwoniffour} leg filename\n");
		}

		geometryParams.filename = vstr[2];

		if (gName == "feas") {
			geometryParams.type = GeometryLibraryType::FEAS;
			if (vstr.size()==4) geometryParams.orbitals=atoi(vstr[3].c_str());
			else geometryParams.orbitals=2;
			return;
		}

		if (gName == "feas1d") {
			geometryParams.type = GeometryLibraryType::FEAS1D;
			if (vstr.size()==4) geometryParams.orbitals=atoi(vstr[3].c_str());
			else geometryParams.orbitals=2;
			return;
		}

		if (gName == "kniffour") {
			geometryParams.type = GeometryLibraryType::KTWONIFFOUR;
			geometryParams.isPeriodicY = geometryParams.leg;
			return;
		}
	}
}; // DriverHelper
} // namespace Dmrg 

/*@}*/
#endif // DRIVER_HELPER_H
