/*
Copyright (c) 2016-2017, UT-Battelle, LLC

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
#ifndef ParametersForMera_H
#define ParametersForMera_H
#include "Vector.h"
#include "IoSimple.h"

namespace Mera {

template<typename ComplexOrRealType_>
struct ParametersForMera {

	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	ParametersForMera(const VectorType& hTerms,
	                  const VectorSizeType& qOne1,
	                  SizeType m1,
	                  PsimagLite::String eval,
	                  PsimagLite::String model1,
	                  RealType tol)
	    : hamiltonianConnection(hTerms),
	      qOne(qOne1),
	      m(m1),
	      verbose(false),
	      evaluator(eval),
	      model(model1),
	      tolerance(tol)
	{}

	ParametersForMera(PsimagLite::String filename)
	{
		PsimagLite::IoSimple::In io(filename);
		io.readline(options,"MeraOptions=");
		io.read(hamiltonianConnection, "hamiltonianConnection");
		io.read(qOne, "qOne");
		io.readline(m, "m=");
		int x = 0;
		io.readline(x, "verbose=");
		verbose = (x > 0);
		io.readline(evaluator, "evaluator=");
		io.readline(model, "Model=");
		io.readline(tolerance, "Tolerance=");
	}

	PsimagLite::String options;
	VectorType hamiltonianConnection;
	VectorSizeType qOne;
	SizeType m;
	bool verbose;
	PsimagLite::String evaluator;
	PsimagLite::String model;
	RealType tolerance;
}; // struct ParametersForMera

template<typename T>
std::ostream& operator<<(std::ostream& os, const ParametersForMera<T>& p)
{
	os<<"MeraOptions=none\n";
	os<<"hamiltonianConnection ";
	os<<p.hamiltonianConnection;
	os<<"qOne\n";
	os<<p.qOne<<"\n";
	os<<"m="<<p.m<<"\n";
	os<<"verbose="<<((p.verbose) ? 1 : 0)<<"\n";
	os<<"evaluator="<<p.evaluator<<"\n";
	os<<"Tolerance="<<p.tolerance<<"\n";
	return os;
}

} // namespace Mera
#endif // ParametersForMera_H
