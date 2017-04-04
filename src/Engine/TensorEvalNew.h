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
#ifndef TENSOREVALNEW_H
#define TENSOREVALNEW_H
#include "TensorEvalBase.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorEvalNew : public TensorEvalBase<ComplexOrRealType> {

	typedef TensorEvalBase<ComplexOrRealType> TensorEvalBaseType;
	typedef typename TensorEvalBaseType::SrepStatementType SrepStatementType;
	typedef typename TensorEvalBaseType::HandleType HandleType;
	typedef typename TensorEvalBaseType::TensorType TensorType;
	typedef typename TensorEvalBaseType::VectorTensorType VectorTensorType;
	typedef typename TensorEvalBaseType::VectorSizeType VectorSizeType;
	typedef typename TensorEvalBaseType::PairStringSizeType PairStringSizeType;
	typedef typename TensorEvalBaseType::MapPairStringSizeType MapPairStringSizeType;
	typedef typename TensorEvalBaseType::VectorPairStringSizeType VectorPairStringSizeType;

public:

	TensorEvalNew(const SrepStatementType& tSrep,
	              const VectorTensorType& vt,
	              const VectorPairStringSizeType& tensorNameIds,
	              MapPairStringSizeType& nameIdsTensor)
	{}

	HandleType operator()()
	{
		throw PsimagLite::RuntimeError("Not implemented yet\n");
	}

	void printResult(std::ostream& os) const
	{
		throw PsimagLite::RuntimeError("Not implemented yet\n");
	}
};
} // namespace Mera
#endif // TENSOREVALNEW_H
