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
#include <map>
#include "Tokenizer.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorEval {

	typedef TensorSrep TensorSrepType;

	class SrepEquation {

	public:

		typedef Tensor<ComplexOrRealType> TensorType;
		typedef typename TensorType::VectorSizeType VectorSizeType;
		typedef typename PsimagLite::Vector<TensorType*>::Type VectorTensorType;
		typedef std::pair<SizeType, SizeType> PairSizeType;
		typedef PsimagLite::Vector<PairSizeType>::Type VectorPairSizeType;
		typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
		typedef std::pair<PsimagLite::String,SizeType> PairStringSizeType;
		typedef PsimagLite::Vector<PairStringSizeType>::Type VectorPairStringSizeType;
		typedef std::map<PairStringSizeType,SizeType> MapPairStringSizeType;
		typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

		SrepEquation(PsimagLite::String str,
		             const VectorTensorType& vt,
		             const VectorPairStringSizeType& tensorNameIds,
		             MapPairStringSizeType& nameIdsTensor)
		    : lhs_(0),rhs_(0),outputTensor_(0)
		{
			VectorStringType vstr;
			PsimagLite::tokenizer(str,vstr,"=");
			if (vstr.size() != 2)
				PsimagLite::RuntimeError("SrepEquation:: syntax error " + str + "\n");
			lhs_ = new TensorSrepType(vstr[0]);
			rhs_ = new TensorSrepType(vstr[1]);

			if (lhs_->size() != 1)
				PsimagLite::RuntimeError("SrepEquation:: LHS should have exactly 1 tensor\n");

			PairStringSizeType nameIdOfOutput(lhs_->operator ()(0).name(),
			                                  lhs_->operator ()(0).id());
			SizeType outputTensorIndex = nameIdsTensor[nameIdOfOutput];
			assert(outputTensorIndex < vt.size());
			outputTensor_ = vt[outputTensorIndex];
		}

		~SrepEquation()
		{
			delete lhs_;
			delete rhs_;
		}

		void fillOutput(const VectorSizeType& free,
		                ComplexOrRealType value)
		{
			outputTensor_->operator()(free) = value;
		}

		const TensorSrepType& lhs() const
		{
			assert(lhs_);
			return *lhs_;
		}

		const TensorSrepType& rhs() const
		{
			assert(rhs_);
			return *rhs_;
		}

		const TensorType& outputTensor() const
		{
			assert(outputTensor_);
			return *outputTensor_;
		}

	private:

		SrepEquation(const SrepEquation&);

		SrepEquation& operator=(const SrepEquation&);

		TensorSrepType* lhs_;
		TensorSrepType* rhs_;
		TensorType* outputTensor_;
	};

	class TensorEvalHandle {


	public:

		enum Status {STATUS_IDLE, STATUS_IN_PROGRESS, STATUS_DONE};

		TensorEvalHandle(Status status = STATUS_IDLE)
		    : status_(status)
		{}

		bool done() const
		{
			return (status_ == STATUS_DONE);
		}

	private:

		Status status_;
	};

public:

	typedef SrepEquation SrepEquationType;
	typedef TensorEvalHandle HandleType;
	typedef typename SrepEquationType::VectorTensorType VectorTensorType;
	typedef typename SrepEquationType::VectorPairStringSizeType VectorPairStringSizeType;
	typedef typename SrepEquationType::MapPairStringSizeType MapPairStringSizeType;
	typedef typename SrepEquationType::VectorSizeType VectorSizeType;
	typedef typename SrepEquationType::PairStringSizeType PairStringSizeType;
	typedef typename SrepEquationType::TensorType TensorType;

	TensorEval(PsimagLite::String srep,
	           const VectorTensorType& vt,
	           const VectorPairStringSizeType& tensorNameIds,
	           MapPairStringSizeType& nameIdsTensor)
	    : srepEq_(new SrepEquationType(srep,vt,tensorNameIds,nameIdsTensor)),
	      ownsSrep_(true),
	      data_(vt),
	      tensorNameIds_(tensorNameIds),
	      nameIdsTensor_(nameIdsTensor)
	{}

	TensorEval(SrepEquationType& tSrep,
	           const VectorTensorType& vt,
	           const VectorPairStringSizeType& tensorNameIds,
	           MapPairStringSizeType& nameIdsTensor)
	    : srepEq_(&tSrep),
	      ownsSrep_(false),
	      data_(vt),
	      tensorNameIds_(tensorNameIds),
	      nameIdsTensor_(nameIdsTensor)
	{}

	~TensorEval()
	{
		if (ownsSrep_) delete srepEq_;
	}

	HandleType operator()()
	{
		SizeType total = srepEq_->outputTensor().args();
		VectorSizeType dimensions(total,0);
		for (SizeType i = 0; i < total; ++i)
			dimensions[i] = srepEq_->outputTensor().argSize(i);

		VectorSizeType free(total,0);

		do {
			srepEq_->fillOutput(free,slowEvaluator(free,srepEq_->rhs()));
		} while (nextIndex(free,dimensions));

		HandleType handle(HandleType::STATUS_DONE);
		return handle;
	}

	void printResult(std::ostream& os) const
	{
		SizeType total = srepEq_->outputTensor().args();
		VectorSizeType dimensions(total,0);
		for (SizeType i = 0; i < total; ++i)
			dimensions[i] = srepEq_->outputTensor().argSize(i);

		VectorSizeType free(total,0);

		do {
			SizeType index = srepEq_->outputTensor().index(free);
			std::cout<<index<<" "<<srepEq_->outputTensor()(free)<<"\n";
		} while (nextIndex(free,dimensions));
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

	ComplexOrRealType slowEvaluator(const VectorSizeType& free,
	                                const TensorSrepType& srep)
	{
		SizeType total = srep.maxTag('s') + 1;
		VectorSizeType summed(total,0);
		VectorSizeType dimensions(total,0);

		prepare(dimensions,srep);

		ComplexOrRealType sum = 0.0;
		do {
			sum += evalInternal(summed,free,srep);
		} while (nextIndex(summed,dimensions));

		return sum;
	}

	void prepare(VectorSizeType& dimensions, const TensorSrepType& tensorSrep)
	{
		SizeType ntensors = tensorSrep.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType id = tensorSrep(i).id();
			SizeType mid = idNameToIndex(tensorSrep(i).name(),id);
			assert(mid < data_.size());
			SizeType ins = tensorSrep(i).ins();
			for (SizeType j = 0; j < ins; ++j) {
				if (tensorSrep(i).legType(j,TensorStanza::INDEX_DIR_IN) !=
				        TensorStanza::INDEX_TYPE_SUMMED) continue;
				SizeType sIndex = tensorSrep(i).legTag(j,TensorStanza::INDEX_DIR_IN);

				assert(j < data_[mid]->args());
				assert(sIndex < dimensions.size());
				dimensions[sIndex] = data_[mid]->argSize(j);
			}

			SizeType outs = tensorSrep(i).outs();
			for (SizeType j = 0; j < outs; ++j) {
				if (tensorSrep(i).legType(j,TensorStanza::INDEX_DIR_OUT) !=
				        TensorStanza::INDEX_TYPE_SUMMED) continue;
				SizeType sIndex = tensorSrep(i).legTag(j,TensorStanza::INDEX_DIR_OUT);

				assert(sIndex < dimensions.size());
				assert(j + ins < data_[mid]->args());
				dimensions[sIndex] = data_[mid]->argSize(j+ins);
			}
		}
	}

	ComplexOrRealType evalInternal(const VectorSizeType& summed,
	                               const VectorSizeType& free,
	                               const TensorSrepType& tensorSrep)
	{
		ComplexOrRealType prod = 1.0;
		SizeType ntensors = tensorSrep.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			prod *= evalThisTensor(tensorSrep(i),summed,free);
			if (prod == 0) break;
		}

		return prod;
	}

	ComplexOrRealType evalThisTensor(const TensorStanza& ts,
	                                 const VectorSizeType& summed,
	                                 const VectorSizeType& free)
	{
		SizeType id = ts.id();
		SizeType mid = idNameToIndex(ts.name(),id);
		SizeType ins = ts.ins();
		SizeType outs = ts.outs();
		assert(data_[mid]->args() == ins + outs);
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
			default:
				PsimagLite::RuntimeError("evalThisTensor: Wrong index type\n");
			}
		}

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
			default:
				PsimagLite::RuntimeError("evalThisTensor: Wrong index type\n");
			}
		}

		return data_[mid]->operator()(args);
	}

	SizeType idNameToIndex(PsimagLite::String name, SizeType id)
	{
		return nameIdsTensor_[PairStringSizeType(name,id)];
	}

	TensorEval(const TensorEval& other);

	TensorEval& operator=(const TensorEval& other);

	SrepEquationType* srepEq_;
	bool ownsSrep_;
	const VectorTensorType& data_;
	const VectorPairStringSizeType& tensorNameIds_;
	MapPairStringSizeType& nameIdsTensor_;
};
}
#endif // MERA_TENSOREVAL_H
