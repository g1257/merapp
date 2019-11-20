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

#include "Tensor.h"
#include "SrepStatement.h"
#include "SymmetryLocal.h"
#include "NameToIndexLut.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorEval {

public:

	typedef bool HandleType;
	typedef Tensor<ComplexOrRealType> TensorType;
	typedef typename PsimagLite::Vector<TensorType*>::Type VectorTensorType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef SrepStatement<ComplexOrRealType> SrepStatementType;
	typedef typename SrepStatementType::PairStringSizeType PairStringSizeType;
	typedef std::map<PairStringSizeType,SizeType> MapPairStringSizeType;
	typedef typename PsimagLite::Vector<PairStringSizeType>::Type VectorPairStringSizeType;
	typedef SymmetryLocal SymmetryLocalType;

	TensorEval(const SrepStatementType& tSrep,
	           const VectorTensorType&,
	           const NameToIndexLut<TensorType>&,
	           SymmetryLocalType*,
	           bool = false) : srepStatement_(tSrep)
	{}

	SizeType nameToIndexLut(PsimagLite::String name)
	{
		throw PsimagLite::RuntimeError("TensorEvalNew::nameToIndexLut Not implemented yet\n");
	}

	HandleType operator()()
	{
		PsimagLite::String copy = srepStatement_.sRep();
		filterForExatn(copy);
		std::cout<<copy<<"\n";
		//Evaluate a tensor network:
		auto evaluated = exatn::evaluateTensorNetwork("srepStatement_.sRep()" , copy);
		TensorType::checkTalshErrorCode(evaluated, "evaluateTensorNetwork");
		PsimagLite::String lhs = srepStatement_.lhs().sRep();
		filterForExatn(lhs);
		auto synced = exatn::sync(lhs);
		return synced;
	}

	void printResult(std::ostream& os) const
	{
		throw PsimagLite::RuntimeError("TensorEvalNew::printResult: Not implemented yet\n");
	}

private:

	static void filterForExatn(PsimagLite::String& str)
	{
		removeAll(str, '*');
		replaceAll(str, '|', ',');
		addMultiplication(str);
	}

	static void removeAll(PsimagLite::String& str, char what)
	{
		PsimagLite::String buffer;
		const SizeType n = str.length();
		for (SizeType i = 0; i < n; ++i)
			if (str[i] != what)
				buffer += str[i];
		str	= buffer;
	}

	static void replaceAll(PsimagLite::String& str, char this1, char that1)
	{
		PsimagLite::String buffer;
		const SizeType n = str.length();
		for (SizeType i = 0; i < n; ++i) {
			if (str[i] == this1)
				buffer += that1;
			else
				buffer += str[i];
		}

		str	= buffer;
	}

	static void addMultiplication(PsimagLite::String& str)
	{
		PsimagLite::String buffer;
		const SizeType n = str.length();
		for (SizeType i = 0; i < n; ++i) {
			if (str[i] == ')' && followsSomething(str, i + 1))
				buffer += ")*";
			else
				buffer += str[i];
		}

		str	= buffer;
	}

	static bool followsSomething(const PsimagLite::String& str, SizeType ind)
	{
		const SizeType n = str.length();
		SizeType i = ind;
		for (; i < n; ++i)
			if (str[i] != ' ' && str[i] != '\t')
				break;

		if (i >= n) return false;

		if (str[i] == '=') return false;

		return true;
	}

	const SrepStatementType& srepStatement_;
};
} // namespace Mera
#endif // TENSOREVALNEW_H
