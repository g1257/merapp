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
#ifndef MERASOLVER_H
#define MERASOLVER_H
#include <iostream>
#include "ParametersForSolver.h"
#include "TensorSrep.h"

namespace Mera {

class MeraSolver {

	typedef ParametersForSolver ParametersForSolverType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	MeraSolver(PsimagLite::String srep, const ParametersForSolver& params)
	    : params_(params), tensorSrep_(srep)
	{}

	void computeEnvirons()
	{
		for (SizeType i = 0; i < tensorSrep_.size(); ++i)
			environForTensor(i);
	}

	friend std::ostream& operator<<(std::ostream& os, const MeraSolver& ms)
	{
		os<<ms.tensorSrep_<<"\n";
		return os;
	}

private:

	MeraSolver(const MeraSolver&);

	MeraSolver& operator=(const MeraSolver&);

	// find Y (environment) for this tensor
	void environForTensor(SizeType ind) const
	{
		if (tensorSrep_(ind).name() == "r") return;
		SizeType sites = tensorSrep_.maxTag('f') + 1;
		std::cout<<"Y_FOR_TENSOR "<<ind<<"\n";
		for (SizeType site = 0; site < sites; ++site) {
			environForTensor(ind,site);
		}

		std::cout<<"\n";
	}

	void environForTensor(SizeType ind, SizeType site) const
	{
		Mera::TensorSrep tensorSrep2(tensorSrep_);
		tensorSrep2.conjugate();
		PsimagLite::String str3("h0(f");
		str3 += ttos(site+2) + ",f";
		str3 += ttos(site+3) + "|f";
		str3 += ttos(site) + ",f";
		str3 += ttos(site+1) + ")\n";
		Mera::TensorSrep tensorSrep3(str3);
		Mera::TensorSrep::VectorSizeType indicesToContract(2,0);
		indicesToContract[1] = 1;
		Mera::TensorSrep tensorSrep4(tensorSrep_);
		tensorSrep4.contract(tensorSrep3,indicesToContract);
		tensorSrep4.isValid(true);
		tensorSrep4.contract(tensorSrep2);
		tensorSrep4.isValid(true);
		tensorSrep4.eraseTensor(ind);
		std::cout<<"STRING: "<<tensorSrep4.sRep()<<"\n";
	}

	const ParametersForSolver& params_;
	TensorSrep tensorSrep_;
}; //class

} //namespace

#endif
