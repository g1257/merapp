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
#ifndef MERA_TensorEvalSlow_H
#define MERA_TensorEvalSlow_H

#include "TensorSrep.h"
#include <map>
#include "Tokenizer.h"
#include "SrepEquation.h"
#include "TensorBreakup.h"
#include "TensorEvalBase.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorEvalSlow : public TensorEvalBase<ComplexOrRealType> {

	typedef TensorSrep TensorSrepType;

public:

	typedef TensorEvalBase<ComplexOrRealType> TensorEvalBaseType;
	typedef typename TensorEvalBaseType::SrepEquationType SrepEquationType;
	typedef typename TensorEvalBaseType::HandleType HandleType;
	typedef typename TensorEvalBaseType::TensorType TensorType;
	typedef typename TensorEvalBaseType::VectorTensorType VectorTensorType;
	typedef typename TensorEvalBaseType::VectorSizeType VectorSizeType;
	typedef typename TensorEvalBaseType::PairStringSizeType PairStringSizeType;
	typedef typename TensorEvalBaseType::MapPairStringSizeType MapPairStringSizeType;
	typedef typename PsimagLite::Vector<SrepEquationType*>::Type VectorSrepEquationType;
	typedef TensorBreakup::VectorStringType VectorStringType;
	typedef typename PsimagLite::Vector<PairStringSizeType>::Type VectorPairStringSizeType;

	static const SizeType EVAL_BREAKUP = TensorBreakup::EVAL_BREAKUP;

	TensorEvalSlow(SrepEquationType& tSrep,
	           const VectorTensorType& vt,
	           const VectorPairStringSizeType& tensorNameIds,
	           MapPairStringSizeType& nameIdsTensor,
	           bool modify)
	    : srepEq_(tSrep),
	      data_(vt), // deep copy
	      tensorNameIds_(tensorNameIds), // deep copy
	      nameIdsTensor_(nameIdsTensor) // deep copy

	{
		indexOfOutputTensor_ = indexOfOutputTensor(tSrep, tensorNameIds, nameIdsTensor);

		if (!modify) return;

		TensorBreakup tensorBreakup(srepEq_.lhs(), srepEq_.rhs());
		// get t0, t1, etc definitions and result
		VectorStringType vstr;
		tensorBreakup(vstr);
		//		PsimagLite::String brokenResult = tensorBreakup.brokenResult();
		// loop over temporaries definitions
		assert(!(vstr.size() & 1));
		SizeType outputLocation = 1 + vstr.size();
		for (SizeType i = 0; i < vstr.size(); i += 2) {
			// add them to tensorNameIds nameIdsTensor
			PsimagLite::String temporaryName = vstr[i];
			if (temporaryName == tSrep.lhs().sRep()) {
				//				vstr[i + 1] = brokenResult;
				std::cout<<"Definition of "<<srepEq_.rhs().sRep()<<" is ";
				std::cout<<vstr[i + 1]<<"\n";
				srepEq_.rhs() = TensorSrep(vstr[i + 1]);
				outputLocation = i;
			}

			if (temporaryName[0] != 't') continue;
			PsimagLite::String str = temporaryName.substr(1,temporaryName.length());
			SizeType temporaryId = atoi(str.c_str());
			temporaryName = "t";
			PairStringSizeType tmpPair = PairStringSizeType(temporaryName,temporaryId);
			tensorNameIds_.push_back(tmpPair);
			nameIdsTensor_[tmpPair] = tensorNameIds_.size() - 1;

			// add this temporary, call setDimensions for output tensor later
			TensorStanza tmpStanza(vstr[i]);
			VectorSizeType args(1,1); // bogus
			TensorType* t = new TensorType(args, tmpStanza.ins());
			garbage_.push_back(t);
			data_.push_back(t);
		}

		VectorSrepEquationType veqs;
		TensorSrepType::VectorPairSizeType empty;
		for (SizeType i = 0; i < vstr.size(); i += 2) {
			veqs.push_back(new SrepEquationType(vstr[i] + "=" + vstr[i+1]));
			SizeType j = veqs.size() - 1;
			if (i != outputLocation)
				veqs[j]->canonicalize();
			veqs[j]->rhs().simplify(empty);

			TensorEvalSlow tEval(*(veqs[j]),data_,tensorNameIds_,nameIdsTensor_,false);

			std::cerr<<"Evaluation of "<<veqs[j]->sRep()<<"\n";
			tEval(false); //handle the handle here
		}

		std::cout.flush();
		for (SizeType i = 0; i < veqs.size(); ++i) {
			delete veqs[i];
			veqs[i] = 0;
		}
	}

	~TensorEvalSlow()
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

		SizeType total = srepEq_.lhs().maxTag('f') + 1;

		static VectorSizeType dimensions;
		if (total != dimensions.size()) dimensions.resize(total,0);
		else std::fill(dimensions.begin(), dimensions.end(), 0);

		bool hasFree = srepEq_.lhs().hasLegType('f');
		if (hasFree) {
			prepare(dimensions,srepEq_.rhs(),TensorStanza::INDEX_TYPE_FREE);
		} else {
			assert(dimensions.size() == 1);
			dimensions[0] = 1;
		}

		static VectorSizeType free;
		if (total != free.size()) free.resize(total,0);
		else std::fill(free.begin(), free.end(), 0);

		if (dimensions.size() == 1 && dimensions[0] == 0)
			dimensions[0] = 1;
		outputTensor().setSizes(dimensions);

		do {
			outputTensor()(free) = slowEvaluator(free,srepEq_.rhs());
		} while (nextIndex(free,dimensions,total));

		return handle;
	}

	void printResult(std::ostream& os) const
	{
		SizeType total = outputTensor().args();
		static VectorSizeType dimensions;
		if (total > dimensions.size()) dimensions.resize(total,0);

		for (SizeType i = 0; i < total; ++i)
			dimensions[i] = outputTensor().argSize(i);

		static VectorSizeType free;
		if (total > free.size()) free.resize(total,0);
		else std::fill(free.begin(), free.end(), 0);

		do {
			SizeType index = outputTensor().index(free);
			std::cout<<index<<" "<<outputTensor()(free)<<"\n";
		} while (nextIndex(free,dimensions,total));
	}

	static bool nextIndex(VectorSizeType& summed,
	                      const VectorSizeType& dimensions,
	                      SizeType total)
	{
		assert(total <= summed.size());
		for (SizeType i = 0; i < total; ++i)
			assert(dimensions[i] == 0 || summed[i] < dimensions[i]);

		for (SizeType i = 0; i < total; ++i) {
			summed[i]++;
			if (summed[i] < dimensions[i]) break;
			summed[i] = 0;
			if (i + 1 == total) return false;
		}

		for (SizeType i = 0; i < total; ++i)
			assert(dimensions[i] == 0 || summed[i] < dimensions[i]);

		return true;
	}

	static SizeType indexOfOutputTensor(const SrepEquationType& eq,
	                                    const VectorPairStringSizeType& tensorNameIds,
	                                    MapPairStringSizeType& nameIdsTensor)
	{
		SizeType ret = nameIdsTensor[eq.nameIdOfOutput()];
		if (tensorNameIds[ret] != eq.nameIdOfOutput()) {
			PsimagLite::String msg("SrepEquation: Could not find ");
			msg += "output tensor " + eq.nameIdOfOutput().first;
			msg += ttos(eq.nameIdOfOutput().second) + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

		return ret;
	}

private:

	ComplexOrRealType slowEvaluator(const VectorSizeType& free,
	                                const TensorSrepType& srep)
	{
		SizeType total = srep.maxTag('s') + 1;
		static VectorSizeType summed;
		if (summed.size() != total) summed.resize(total,0);
		else std::fill(summed.begin(), summed.end(), 0);

		static VectorSizeType dimensions;
		if (dimensions.size() != total) dimensions.resize(total,0);
		else std::fill(dimensions.begin(), dimensions.end(), 0);

		bool hasSummed = srep.hasLegType('s');
		if (hasSummed) {
			prepare(dimensions, srep, TensorStanza::INDEX_TYPE_SUMMED);
		} else {
			assert(dimensions.size() == 1);
			dimensions[0] = 1;
		}

		ComplexOrRealType sum = 0.0;
		do {
			sum += evalInternal(summed,free,srep);
		} while (nextIndex(summed,dimensions,total));

		return sum;
	}

	void prepare(VectorSizeType& dimensions,
	             const TensorSrepType& tensorSrep,
	             TensorStanza::IndexTypeEnum type) const
	{
		SizeType ntensors = tensorSrep.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			prepareStanza(dimensions, tensorSrep(i), type);
		}
	}

	void prepareStanza(VectorSizeType& dimensions,
	                   const TensorStanza& stanza,
	                   TensorStanza::IndexTypeEnum type) const
	{
		SizeType id = stanza.id();
		SizeType mid = idNameToIndex(stanza.name(),id);
		assert(mid < data_.size());
		SizeType legs = stanza.legs();
		for (SizeType j = 0; j < legs; ++j) {
			if (stanza.legType(j) != type)
				continue;
			SizeType sIndex = stanza.legTag(j);

			assert(j < data_[mid]->args());
			assert(sIndex < dimensions.size());
			dimensions[sIndex] = data_[mid]->argSize(j);
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
		SizeType legs = ts.legs();
		assert(legs == 0 || data_[mid]->args() == legs);

		VectorSizeType args(data_[mid]->args(),0);

		for (SizeType j = 0; j < legs; ++j) {
			SizeType index = ts.legTag(j);

			switch (ts.legType(j)) {

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

		return data_[mid]->operator()(args);
	}

	SizeType idNameToIndex(PsimagLite::String name, SizeType id) const
	{
		typename MapPairStringSizeType::iterator it =
		        nameIdsTensor_.find(PairStringSizeType(name,id));
		if (it == nameIdsTensor_.end())
			throw PsimagLite::RuntimeError("idNameToIndex: key not found\n");
		return it->second;
	}

	TensorType& outputTensor()
	{
		assert(indexOfOutputTensor_ < data_.size());
		return *(data_[indexOfOutputTensor_]);
	}

	const TensorType& outputTensor() const
	{
		assert(indexOfOutputTensor_ < data_.size());
		return *(data_[indexOfOutputTensor_]);
	}

	TensorEvalSlow(const TensorEvalSlow& other);

	TensorEvalSlow& operator=(const TensorEvalSlow& other);

	SrepEquationType srepEq_;
	VectorTensorType data_;
	VectorPairStringSizeType tensorNameIds_;
	mutable MapPairStringSizeType nameIdsTensor_;
	SizeType indexOfOutputTensor_;
	VectorTensorType garbage_;
};
}
#endif // MERA_TensorEvalSlow_H
