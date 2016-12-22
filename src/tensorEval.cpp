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
#include "TensorEval.h"
#include "Vector.h"
#include "SrepEquation.h"

int main()
{
	PsimagLite::String str = "r0(f0) = u0(f0|s0)u1(s0)";

	SizeType dim0 = 5;
	typedef Mera::TensorEval<double> TensorEvalType;
	typedef TensorEvalType::TensorType TensorType;
	TensorEvalType::VectorTensorType vt(3);

	TensorEvalType::VectorSizeType d(2,dim0);
	d[0] = 3;
	vt[0] = new TensorType(d,1);
	vt[1] = new TensorType(dim0,1);
	vt[2] = new TensorType(3,1);
	for (SizeType i = 0; i < vt.size(); ++i) {
		vt[i]->setToIdentity(1.0);
	}

	vt[1]->setToConstant(1.5);

	TensorEvalType::VectorPairStringSizeType idNames;
	idNames.push_back(TensorEvalType::PairStringSizeType("u",0));
	idNames.push_back(TensorEvalType::PairStringSizeType("u",1));
	idNames.push_back(TensorEvalType::PairStringSizeType("r",0));
	TensorEvalType::MapPairStringSizeType nameIdTensor;
	for (SizeType i = 0; i < vt.size(); ++i)
		nameIdTensor[idNames[i]] = i;

	Mera::SrepEquation<double> srepEq(str,vt,idNames,nameIdTensor);
	TensorEvalType tensorEval(srepEq,vt,idNames,nameIdTensor,false);
	TensorEvalType::HandleType handle = tensorEval();

	while (!handle.done());

	tensorEval.printResult(std::cout);

	for (SizeType i = 0; i < vt.size(); ++i) {
		delete vt[i];
		vt[i] = 0;
	}
}
