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
#include "SrepEquation.h"
#include "TensorBreakup.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorEval {

	typedef TensorSrep TensorSrepType;

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

	typedef SrepEquation<ComplexOrRealType> SrepEquationType;
	typedef typename PsimagLite::Vector<SrepEquationType*>::Type VectorSrepEquationType;
	typedef TensorBreakup::VectorStringType VectorStringType;
	typedef TensorEvalHandle HandleType;
	typedef typename SrepEquationType::VectorTensorType VectorTensorType;
	typedef typename SrepEquationType::VectorPairStringSizeType VectorPairStringSizeType;
	typedef typename SrepEquationType::MapPairStringSizeType MapPairStringSizeType;
	typedef typename SrepEquationType::VectorSizeType VectorSizeType;
	typedef typename SrepEquationType::PairStringSizeType PairStringSizeType;
	typedef typename SrepEquationType::TensorType TensorType;

	TensorEval(SrepEquationType& tSrep,
	           const VectorTensorType& vt,
	           const VectorPairStringSizeType& tensorNameIds,
	           MapPairStringSizeType& nameIdsTensor,
	           bool modify)
	    : srepEq_(tSrep),
	      data_(vt), // deep copy
	      tensorNameIds_(tensorNameIds), // deep copy
	      nameIdsTensor_(nameIdsTensor) // deep copy
	{
		if (!modify) return;
		TensorBreakup tensorBreakup(srepEq_.lhs(), srepEq_.rhs());
		// get t0, t1, etc definitions and result
		VectorStringType vstr;
		tensorBreakup(vstr);

		// loop over temporaries definitions
		assert(!(vstr.size() & 1));
		for (SizeType i = 0; i < vstr.size(); i += 2) {
			// add them to tensorNameIds nameIdsTensor
			PsimagLite::String temporaryName = vstr[i];
			if (temporaryName == tSrep.lhs().sRep()) {
				std::cout<<"Definition of "<<srepEq_.rhs().sRep()<<" is ";
				std::cout<<vstr[i + 1]<<"\n";
				srepEq_.rhs() = TensorSrep(vstr[i + 1]);
			}

			if (temporaryName[0] != 't') continue;
			PsimagLite::String str = temporaryName.substr(1,temporaryName.length());
			SizeType temporaryId = atoi(str.c_str());
			temporaryName = "t";
			tensorNameIds_.push_back(PairStringSizeType(temporaryName,temporaryId));
			// warning: nameIdsTensor_ is out of sync with tensorNameIds_

			// add storage for this temporary
			VectorSizeType args;
			SizeType ins = findArgsAndIns(args,vstr[i],vstr[i+1]);
			TensorType* t = new TensorType(args, ins);
			garbage_.push_back(t);
			data_.push_back(t);
		}

		// sync nameIdsTensor_ with tensorNameIds_
		for (SizeType i = 0; i < tensorNameIds_.size(); ++i)
			nameIdsTensor_[tensorNameIds_[i]] = i;

		VectorSrepEquationType veqs;

		for (SizeType i = 0; i < vstr.size(); i += 2) {
			veqs.push_back(new SrepEquationType(vstr[i] + "=" + vstr[i+1],
			               data_,
			               tensorNameIds_,
			               nameIdsTensor_));
			SizeType j = veqs.size() - 1;
			veqs[j]->canonicalize();
			TensorEval tEval(*(veqs[j]),
			                 data_,
			                 tensorNameIds_,
			                 nameIdsTensor_,
			                 false);
			std::cerr<<"Evaluation of "<<veqs[j]->sRep()<<"\n";
			tEval(false); //handle the handle here
		}

		std::cout.flush();
		for (SizeType i = 0; i < veqs.size(); ++i) {
			delete veqs[i];
			veqs[i] = 0;
		}
	}

	~TensorEval()
	{
		for (SizeType i = 0; i < garbage_.size(); ++i) {
			delete garbage_[i];
			garbage_[i] = 0;
		}
	}

	HandleType operator()(bool cached)
	{
		HandleType handle(HandleType::STATUS_DONE);
		if (cached) return handle;

		SizeType total = srepEq_.outputTensor().args();
		static VectorSizeType dimensions;
		if (total > dimensions.size()) dimensions.resize(total,0);
		// fillFreeDimensions(dimensions);
		for (SizeType i = 0; i < total; ++i)
			dimensions[i] = srepEq_.outputTensor().argSize(i);

		static VectorSizeType free;
		if (total > free.size()) free.resize(total,0);
		else std::fill(free.begin(), free.end(), 0);

		do {
			srepEq_.fillOutput(free,slowEvaluator(free,srepEq_.rhs()));
		} while (nextIndex(free,dimensions,total));


		return handle;
	}

	void printResult(std::ostream& os) const
	{
		SizeType total = srepEq_.outputTensor().args();
		static VectorSizeType dimensions;
		if (total > dimensions.size()) dimensions.resize(total,0);

		for (SizeType i = 0; i < total; ++i)
			dimensions[i] = srepEq_.outputTensor().argSize(i);

		static VectorSizeType free;
		if (total > free.size()) free.resize(total,0);
		else std::fill(free.begin(), free.end(), 0);

		do {
			SizeType index = srepEq_.outputTensor().index(free);
			std::cout<<index<<" "<<srepEq_.outputTensor()(free)<<"\n";
		} while (nextIndex(free,dimensions,total));
	}

	static bool nextIndex(VectorSizeType& summed,
	                      const VectorSizeType& dimensions,
	                      SizeType total)
	{
		assert(total <= summed.size());
		for (SizeType i = 0; i < total; ++i) {
			summed[i]++;
			if (summed[i] < dimensions[i]) break;
			if (i +1 == total) return false;
			summed[i] = 0;
		}

		return true;
	}

private:

	ComplexOrRealType slowEvaluator(const VectorSizeType& free,
	                                const TensorSrepType& srep)
	{
		SizeType total = srep.maxTag('s') + 1;
		static VectorSizeType summed;
		if (summed.size() < total) summed.resize(total,0);
		else std::fill(summed.begin(), summed.end(), 0);

		static VectorSizeType dimensions;
		if (dimensions.size() < total) dimensions.resize(total,0);
		else std::fill(dimensions.begin(), dimensions.end(), 0);

		prepare(dimensions,srep);

		ComplexOrRealType sum = 0.0;
		do {
			sum += evalInternal(summed,free,srep);
		} while (nextIndex(summed,dimensions,total));

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

	SizeType findArgsAndIns(VectorSizeType& args,
	                        PsimagLite::String strLeft,
	                        PsimagLite::String strRight) const
	{
		TensorStanza lhs(strLeft);
		TensorSrep rhs(strRight);

		TensorStanza::IndexDirectionEnum in = TensorStanza::INDEX_DIR_IN;
		TensorStanza::IndexDirectionEnum out = TensorStanza::INDEX_DIR_OUT;

		SizeType ins = lhs.ins();
		for (SizeType i = 0; i < ins; ++i) {
			if (lhs.legType(i,in) != TensorStanza::INDEX_TYPE_FREE) continue;
			args.push_back(findDimensionOfFreeLeg(rhs,lhs.legTag(i,in)));
		}

		SizeType outs = lhs.outs();
		for (SizeType i = 0; i < outs; ++i) {
			if (lhs.legType(i,out) != TensorStanza::INDEX_TYPE_FREE) continue;
			args.push_back(findDimensionOfFreeLeg(rhs,lhs.legTag(i,out)));
		}

		return ins;
	}

	SizeType findDimensionOfFreeLeg(const TensorSrep& srep, SizeType leg) const
	{
		TensorStanza::IndexDirectionEnum in = TensorStanza::INDEX_DIR_IN;
		TensorStanza::IndexDirectionEnum out = TensorStanza::INDEX_DIR_OUT;

		SizeType ntensors = srep.size();
		for (SizeType j = 0; j < ntensors; ++j) {
			TensorStanza stanza = srep(j);
			PsimagLite::String name = stanza.name();
			SizeType id = stanza.id();
			SizeType ins = stanza.ins();
			for (SizeType i = 0; i < ins; ++i) {
				if (stanza.legType(i,in) != TensorStanza::INDEX_TYPE_FREE) continue;
				if (stanza.legTag(i,in) != leg) continue;
				SizeType tensorIndex = nameIdsTensor_[PairStringSizeType(name,id)];
				return data_[tensorIndex]->argSize(i);
			}

			SizeType outs = stanza.outs();
			for (SizeType i = 0; i < outs; ++i) {
				if (stanza.legType(i,out) != TensorStanza::INDEX_TYPE_FREE) continue;
				if (stanza.legTag(i,out) != leg) continue;
				SizeType tensorIndex = nameIdsTensor_[PairStringSizeType(name,id)];
				return data_[tensorIndex]->argSize(i+ins);
			}
		}

		throw PsimagLite::RuntimeError("findDimensionOfFreeLeg\n");
	}

	void fillFreeDimensions(VectorSizeType& dimensions) const
	{
		TensorStanza::IndexDirectionEnum in = TensorStanza::INDEX_DIR_IN;
		TensorStanza::IndexDirectionEnum out = TensorStanza::INDEX_DIR_OUT;

		SizeType ins = srepEq_.lhs().ins();
		for (SizeType j = 0; j < ins; ++j) {
			if (srepEq_.lhs().legTag(j,in) != TensorStanza::INDEX_TYPE_FREE)
				continue;
			SizeType ind = srepEq_.lhs().legTag(j,in);
			assert(ind < dimensions.size());
			dimensions[ind] = srepEq_.outputTensor().argSize(j);
		}

		SizeType outs = srepEq_.lhs().outs();
		for (SizeType j = 0; j < outs; ++j) {
			if (srepEq_.lhs().legTag(j,out) != TensorStanza::INDEX_TYPE_FREE)
				continue;
			SizeType ind = srepEq_.lhs().legTag(j,out);
			assert(ind < dimensions.size());
			dimensions[ind] = srepEq_.outputTensor().argSize(j);
		}
	}

	TensorEval(const TensorEval& other);

	TensorEval& operator=(const TensorEval& other);

	SrepEquationType& srepEq_;
	VectorTensorType data_;
	VectorPairStringSizeType tensorNameIds_;
	mutable MapPairStringSizeType nameIdsTensor_;
	VectorTensorType garbage_;
};
}
#endif // MERA_TENSOREVAL_H
