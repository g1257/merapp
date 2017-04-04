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
#ifndef TENSOREVALBASE_H
#define TENSOREVALBASE_H
#include "Tensor.h"
#include "SrepStatement.h"
#include "TensorEvalHandle.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorEvalBase {

public:

	typedef TensorEvalHandle HandleType;
	typedef Tensor<ComplexOrRealType> TensorType;
	typedef typename PsimagLite::Vector<TensorType*>::Type VectorTensorType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef SrepStatement<ComplexOrRealType> SrepStatementType;
	typedef typename SrepStatementType::PairStringSizeType PairStringSizeType;
	typedef std::map<PairStringSizeType,SizeType> MapPairStringSizeType;
	typedef typename PsimagLite::Vector<PairStringSizeType>::Type VectorPairStringSizeType;

	virtual ~TensorEvalBase() {}

	virtual HandleType operator()() = 0;

	virtual void printResult(std::ostream& os) const = 0;

	static SizeType indexOfOutputTensor(const SrepStatementType& eq,
	                                    const VectorPairStringSizeType& tensorNameIds,
	                                    MapPairStringSizeType& nameIdsTensor)
	{
		SizeType ret = nameIdsTensor[eq.nameIdOfOutput()];
		if (tensorNameIds[ret] != eq.nameIdOfOutput()) {
			PsimagLite::String msg("TensorEvalBase: Could not find ");
			msg += "output tensor " + eq.nameIdOfOutput().first;
			msg += ttos(eq.nameIdOfOutput().second) + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

		return ret;
	}
};
} // namespace Mera
#endif // TENSOREVALBASE_H
