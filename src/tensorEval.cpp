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
#include "TensorEvalNew.h"
#include "Vector.h"
#include "SrepEquation.h"

int main()
{
	PsimagLite::String str = "r0(f0) = u0(f0|s0)u1(s0)";
	PsimagLite::String evaluator = "slow"; // or you can say "new" here when ready

	typedef Mera::TensorEvalSlow<double> TensorEvalSlowType;
	typedef Mera::TensorEvalNew<double> TensorEvalNewType;

	SizeType dim0 = 5;
	typedef Mera::TensorEvalBase<double> TensorEvalBaseType;
	typedef TensorEvalBaseType::TensorType TensorType;
	TensorEvalBaseType::VectorTensorType vt(3);

	TensorEvalBaseType::VectorSizeType d(2,dim0);
	d[0] = 3;
	vt[0] = new TensorType(d,1);
	vt[1] = new TensorType(dim0,1);
	vt[2] = new TensorType(3,1);
	for (SizeType i = 0; i < vt.size(); ++i) {
		vt[i]->setToIdentity(1.0);
	}

	vt[1]->setToConstant(1.5);

	TensorEvalBaseType::VectorPairStringSizeType idNames;
	idNames.push_back(TensorEvalBaseType::PairStringSizeType("u",0));
	idNames.push_back(TensorEvalBaseType::PairStringSizeType("u",1));
	idNames.push_back(TensorEvalBaseType::PairStringSizeType("r",0));
	TensorEvalBaseType::MapPairStringSizeType nameIdTensor;
	for (SizeType i = 0; i < vt.size(); ++i)
		nameIdTensor[idNames[i]] = i;

	Mera::SrepEquation<double> srepEq(str);
	TensorEvalBaseType* tensorEval = 0;
	if (evaluator == "slow") {
		tensorEval = new TensorEvalSlowType(srepEq,vt,idNames,nameIdTensor,false);
	} else if(evaluator == "new") {
		tensorEval = new TensorEvalNewType(srepEq,vt,idNames,nameIdTensor,false);
	} else {
		throw PsimagLite::RuntimeError("Unknown evaluator " + evaluator + "\n");
	}

	TensorEvalBaseType::HandleType handle = tensorEval->operator()(false);

	while (!handle.done());

	tensorEval->printResult(std::cout);

	for (SizeType i = 0; i < vt.size(); ++i) {
		delete vt[i];
		vt[i] = 0;
	}
}
