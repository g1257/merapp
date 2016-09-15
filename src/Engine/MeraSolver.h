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

	void computeGroundState()
	{
		SizeType qiter = 1;

		for (SizeType iter=0; iter<qiter; ++iter)
			for (SizeType i = 0; i < tensorSrep_.size(); ++i)
				optimizeTensor(iter,i);
	}

	friend std::ostream& operator<<(std::ostream& os, const MeraSolver& ms)
	{
		os<<ms.tensorSrep_<<"\n";
		return os;
	}

private:

	MeraSolver(const MeraSolver&);

	MeraSolver& operator=(const MeraSolver&);

	void optimizeTensor(SizeType iter, SizeType ind)
	{
		// find Y (environment) for this tensor
		// Y = USV^+
		// w = -VU^+
	}

	const ParametersForSolver& params_;
	TensorSrep tensorSrep_;
}; //class

} //namespace

#endif
