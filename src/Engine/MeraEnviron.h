/*
Copyright (c) 2016, UT-Battelle, LLC

MERA++, Version 0.

This file is part of MERA++.
MERA++ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
MERA++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with MERA++. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MERA_ENVIRON_H
#define MERA_ENVIRON_H
#include <iostream>
#include "ParametersForSolver.h"
#include "TensorSrep.h"

namespace Mera {

class MeraEnviron {

	typedef ParametersForSolver ParametersForSolverType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

public:

	MeraEnviron(PsimagLite::String srep, const ParametersForSolver& params)
	    : params_(params), tensorSrep_(srep)
	{}

	void computeEnvirons()
	{
		std::cout<<"MERA="<<tensorSrep_.sRep()<<"\n";
		for (SizeType i = 0; i < tensorSrep_.size(); ++i)
			environForTensor(i);
	}

	friend std::ostream& operator<<(std::ostream& os, const MeraEnviron& ms)
	{
		os<<ms.tensorSrep_<<"\n";
		return os;
	}

private:

	MeraEnviron(const MeraEnviron&);

	MeraEnviron& operator=(const MeraEnviron&);

	// find Y (environment) for this tensor
	void environForTensor(SizeType ind) const
	{
		if (tensorSrep_(ind).name() == "r") return;
		SizeType id = tensorSrep_(ind).id();
		PsimagLite::String name = tensorSrep_(ind).name();
		SizeType sites = tensorSrep_.maxTag('f');
		VectorStringType vstr(sites,"");
		SizeType terms = 0;
		for (SizeType site = 0; site < sites; ++site) {
			vstr[site] = environForTensor(ind,site);
			if (vstr[site] != "") ++terms;
		}

		std::cout<<"TensorId="<<name<<","<<id<<"\n";
		std::cout<<"Terms="<<terms<<"\n";
		std::cout<<"IgnoreTerm="<<(2*sites+1)<<"\n";
		for (SizeType site = 0; site < sites; ++site)
			if (vstr[site] != "")
				std::cout<<"Environ="<<vstr[site]<<"\n";

		std::cout<<"\n";
	}

	PsimagLite::String environForTensor(SizeType ind, SizeType site) const
	{
		Mera::TensorSrep tensorSrep2(tensorSrep_);
		tensorSrep2.conjugate();
		tensorSrep2.swapFree(0,site);
		tensorSrep2.swapFree(1,site+1);
		PsimagLite::String str3("h0(f");
		str3 += ttos(site+2) + ",f";
		str3 += ttos(site+3) + "|f";
		str3 += ttos(site) + ",f";
		str3 += ttos(site+1) + ")\n";
		Mera::TensorSrep tensorSrep3(str3);
		Mera::TensorSrep::VectorSizeType indicesToContract(2,site);
		indicesToContract[1] = site + 1;
		Mera::TensorSrep tensorSrep4(tensorSrep_);
		tensorSrep4.contract(tensorSrep3,indicesToContract);
		if (!tensorSrep4.isValid(true))
			throw PsimagLite::RuntimeError("Invalid tensor\n");
		tensorSrep4.contract(tensorSrep2);
		std::cerr<<"ENERGY="<<tensorSrep4.sRep()<<"\n";
		if (!tensorSrep4.isValid(true))
			throw PsimagLite::RuntimeError("Invalid tensor\n");
		tensorSrep4.eraseTensor(ind);
		if (!tensorSrep4.isValid(true))
			throw PsimagLite::RuntimeError("Invalid tensor\n");
		if (tensorSrep4.findConjugate(ind) >= tensorSrep4.size()) {
			std::cerr<<"EMPTY_ENVIRON="<<tensorSrep4.sRep()<<"\n";
			return "";
		}

		return tensorSrep4.sRep();
	}

	const ParametersForSolver& params_;
	TensorSrep tensorSrep_;
}; //class

} //namespace

#endif
