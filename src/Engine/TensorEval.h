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
#ifndef MERA_TENSOREVAL_H
#define MERA_TENSOREVAL_H

#include "TensorSrep.h"
#include "Tensor.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorEval {

	typedef TensorSrep TensorSrepType;

public:

	typedef Tensor<ComplexOrRealType> TensorType;
	typedef typename TensorType::VectorSizeType VectorSizeType;
	typedef typename PsimagLite::Vector<TensorType*>::Type VectorTensorType;
	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef PsimagLite::Vector<PairSizeType>::Type VectorPairSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef std::pair<PsimagLite::String,SizeType> PairStringSizeType;
	typedef PsimagLite::Vector<PairStringSizeType>::Type VectorPairStringSizeType;

	TensorEval(PsimagLite::String srep,
	           const VectorTensorType& vt,
	           const VectorPairStringSizeType& tensorNameIds)
	    : tensorSrep_(srep), data_(vt), tensorNameIds_(tensorNameIds)
	{}

	ComplexOrRealType operator()(const VectorSizeType& free) const
	{
		SizeType total = tensorSrep_.maxTag('s') + 1;
		VectorSizeType summed(total,0);
		VectorSizeType dimensions(total,0);

		prepare(dimensions);

		ComplexOrRealType sum = 0.0;
		do {
			sum += evalInternal(summed,free);
			//std::cerr<<summed;
			//std::cerr<<"---------\n";
		} while (nextIndex(summed,dimensions));

		return sum;
	}

	static bool nextIndex(VectorSizeType& summed,
	                      const VectorSizeType& dimensions)
	{
		for (SizeType i = 0; i < summed.size(); ++i) {
			summed[i]++;
			if (summed[i] < dimensions[i]) break;
			if (i +1 == summed.size()) return false;
			summed[i] = 0;
		}

		return true;
	}

private:

	void prepare(VectorSizeType& dimensions) const
	{
		SizeType ntensors = tensorSrep_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType id = tensorSrep_(i).id();
			SizeType mid = idNameToIndex(tensorSrep_(i).name(),id);
			assert(mid < data_.size());
			SizeType ins = tensorSrep_(i).ins();
			for (SizeType j = 0; j < ins; ++j) {
				if (tensorSrep_(i).legType(j,TensorStanza::INDEX_DIR_IN) !=
				        TensorStanza::INDEX_TYPE_SUMMED) continue;
				SizeType sIndex = tensorSrep_(i).legTag(j,TensorStanza::INDEX_DIR_IN);

				assert(j < data_[mid]->args());
				assert(sIndex < dimensions.size());
				dimensions[sIndex] = data_[mid]->argSize(j);
			}

			SizeType outs = tensorSrep_(i).outs();
			for (SizeType j = 0; j < outs; ++j) {
				if (tensorSrep_(i).legType(j,TensorStanza::INDEX_DIR_OUT) !=
				        TensorStanza::INDEX_TYPE_SUMMED) continue;
				SizeType sIndex = tensorSrep_(i).legTag(j,TensorStanza::INDEX_DIR_OUT);

				assert(sIndex < dimensions.size());
				assert(j + ins < data_[mid]->args());
				dimensions[sIndex] = data_[mid]->argSize(j+ins);
			}
		}
	}

	ComplexOrRealType evalInternal(const VectorSizeType& summed,
	                               const VectorSizeType& free) const
	{
		ComplexOrRealType prod = 1.0;
		SizeType ntensors = tensorSrep_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			prod *= evalThisTensor(tensorSrep_(i),summed,free);
			if (prod == 0) break;
		}

		return prod;
	}

	ComplexOrRealType evalThisTensor(const TensorStanza& ts,
	                                 const VectorSizeType& summed,
	                                 const VectorSizeType& free) const
	{
		SizeType id = ts.id();
		SizeType mid = idNameToIndex(ts.name(),id);
		SizeType ins = ts.ins();
		VectorSizeType args(data_[mid]->args(),0);

		for (SizeType j = 0; j < ins; ++j) {
			SizeType index = ts.legTag(j,TensorStanza::INDEX_DIR_IN);

			switch (ts.legType(j,TensorStanza::INDEX_DIR_IN)) {

			case TensorStanza::INDEX_TYPE_SUMMED:
				assert(index < summed.size());
				assert(j < args.size());
				args[j] = summed[index];
				break;

			case TensorStanza::INDEX_TYPE_FREE:
				assert(index < free.size());
				assert(j < args.size());
				args[j] = free[index];
				break;

			case  TensorStanza::INDEX_TYPE_DUMMY:
				assert(j < args.size());
				args[j] = 0;
				break;
			}
		}

		SizeType outs = ts.outs();
		for (SizeType j = 0; j < outs; ++j) {
			SizeType index = ts.legTag(j,TensorStanza::INDEX_DIR_OUT);

			switch (ts.legType(j,TensorStanza::INDEX_DIR_OUT)) {

			case TensorStanza::INDEX_TYPE_SUMMED:
				assert(index < summed.size());
				assert(j+ins < args.size());
				args[j+ins] = summed[index];
				break;

			case TensorStanza::INDEX_TYPE_FREE:
				assert(index < free.size());
				assert(j+ins < args.size());
				args[j+ins] = free[index];
				break;

			case  TensorStanza::INDEX_TYPE_DUMMY:
				assert(j+ins < args.size());
				args[j+ins] = 0;
				break;
			}
		}

		return data_[mid]->operator()(args);
	}

	SizeType idNameToIndex(PsimagLite::String name, SizeType id) const
	{
		PairStringSizeType p(name,id);
		VectorPairStringSizeType::const_iterator it = std::find(tensorNameIds_.begin(),
		                                                        tensorNameIds_.end(),
		                                                        p);
		if (it == tensorNameIds_.end()) {
			PsimagLite::String str("TensorEval: ");
			str += "could not find index for " + name + " " + ttos(id) + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		return it - tensorNameIds_.begin();
	}

	TensorSrepType tensorSrep_;
	const VectorTensorType& data_;
	const VectorPairStringSizeType& tensorNameIds_;
};
}
#endif // MERA_TENSOREVAL_H
