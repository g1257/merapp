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
#include "TensorBreakup.h"
#include "TensorEvalBase.h"
#include "SymmetryLocal.h"

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
	typedef typename TensorEvalBaseType::VectorPairStringSizeType VectorPairStringSizeType;
	typedef typename TensorEvalBaseType::MapPairStringSizeType MapPairStringSizeType;
	typedef typename PsimagLite::Vector<SrepEquationType*>::Type VectorSrepEquationType;
	typedef TensorBreakup::VectorStringType VectorStringType;
	typedef SymmetryLocal SymmetryLocalType;
	typedef SymmetryLocalType::VectorVectorSizeType VectorVectorSizeType;

	static const SizeType EVAL_BREAKUP = TensorBreakup::EVAL_BREAKUP;

	TensorEvalSlow(const SrepEquationType& tSrep,
	               const VectorTensorType& vt,
	               const VectorPairStringSizeType& tensorNameIds,
	               MapPairStringSizeType& nameIdsTensor,
	               SymmetryLocalType* symmLocal,
	               bool modify = EVAL_BREAKUP)
	    : srepEq_(tSrep),
	      data_(vt), // deep copy
	      tensorNameIds_(tensorNameIds), // deep copy
	      nameIdsTensor_(nameIdsTensor), // deep copy
	      symmLocal_(symmLocal),
	      modify_(modify)
	{
		indexOfOutputTensor_ = TensorEvalBaseType::indexOfOutputTensor(tSrep,
		                                                               tensorNameIds,
		                                                               nameIdsTensor);

		if (!modify_) return;

		TensorBreakup tensorBreakup(srepEq_.lhs(), srepEq_.rhs());
		// get t0, t1, etc definitions and result
		VectorStringType vstr;
		tensorBreakup(vstr);

		// loop over temporaries definitions
		assert(!(vstr.size() & 1));
		SizeType outputLocation = 1 + vstr.size();
		for (SizeType i = 0; i < vstr.size(); i += 2) {
			// add them to tensorNameIds nameIdsTensor
			PsimagLite::String temporaryName = vstr[i];
			if (temporaryName == tSrep.lhs().sRep()) {
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

			TensorEvalSlow tEval(*(veqs[j]),
			                     data_,
			                     tensorNameIds_,
			                     nameIdsTensor_,
			                     symmLocal_,
			                     false);

			std::cerr<<"Evaluation of "<<veqs[j]->sRep()<<"\n";
			tEval(); //handle the handle here
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

		for (SizeType i = 0; i < garbage2_.size(); ++i) {
			delete garbage2_[i];
			garbage2_[i] = 0;
		}
	}

	HandleType operator()()
	{
		HandleType handle(HandleType::STATUS_DONE);
		if (modify_ && EVAL_BREAKUP) return handle;

		SizeType total = srepEq_.lhs().maxTag('f') + 1;

		static VectorSizeType dimensions;
		static VectorVectorSizeType q;
		if (total != dimensions.size()) {
			dimensions.resize(total, 0);
			q.resize(total, 0);
		} else {
			std::fill(dimensions.begin(), dimensions.end(), 0);
			std::fill(q.begin(), q.end(), static_cast<VectorSizeType*>(0));
		}

		bool hasFree = srepEq_.lhs().hasLegType('f');
		if (hasFree) {
			prepare(dimensions,q,srepEq_.rhs(),TensorStanza::INDEX_TYPE_FREE);
			setQnsForOutput(q);
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
		} while (ProgramGlobals::nextIndex(free,dimensions,total));

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
		} while (ProgramGlobals::nextIndex(free,dimensions,total));
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

		VectorVectorSizeType q;
		bool hasSummed = srep.hasLegType('s');
		if (hasSummed) {
			prepare(dimensions, q, srep, TensorStanza::INDEX_TYPE_SUMMED);
		} else {
			assert(dimensions.size() == 1);
			dimensions[0] = 1;
		}

		ComplexOrRealType sum = 0.0;
		do {
			sum += evalInternal(summed,free,srep);
		} while (ProgramGlobals::nextIndex(summed,dimensions,total));

		return sum;
	}

	void prepare(VectorSizeType& dimensions,
	             VectorVectorSizeType& q,
	             const TensorSrepType& tensorSrep,
	             TensorStanza::IndexTypeEnum type) const
	{
		SizeType ntensors = tensorSrep.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			prepareStanza(dimensions, q, tensorSrep(i), type);
		}
	}

	void prepareStanza(VectorSizeType& dimensions,
	                   VectorVectorSizeType& q,
	                   const TensorStanza& stanza,
	                   TensorStanza::IndexTypeEnum type) const
	{
		SizeType id = stanza.id();
		SizeType mid = idNameToIndex(stanza.name(),id);
		SizeType tensorIndex = symmLocal_->nameIdToIndex(stanza.name() + ttos(id));
		if (tensorIndex >= symmLocal_->size())
			assert(false);

		assert(mid < data_.size());
		SizeType legs = stanza.legs();
		for (SizeType j = 0; j < legs; ++j) {
			if (stanza.legType(j) != type)
				continue;
			SizeType sIndex = stanza.legTag(j);

			assert(j < data_[mid]->args());
			assert(sIndex < dimensions.size());
			dimensions[sIndex] = data_[mid]->argSize(j);
			if (type == TensorStanza::INDEX_TYPE_FREE) {
				const VectorSizeType* qSrc = symmLocal_->q(tensorIndex, j);
				assert(qSrc);
				assert(sIndex < q.size());
				VectorSizeType* ptr = new VectorSizeType(qSrc->size());
				//garbage2_.push_back(ptr);
				*(ptr) = *qSrc;
				q[sIndex] = ptr;
			}
		}
	}

	ComplexOrRealType evalInternal(const VectorSizeType& summed,
	                               const VectorSizeType& free,
	                               const TensorSrepType& tensorSrep)
	{
		if (symmLocal_ && !symmetriesPass(summed, free, tensorSrep))
			return 0.0;

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

	bool symmetriesPass(const VectorSizeType& summed,
	                    const VectorSizeType& free,
	                    const TensorSrepType& tensorSrep) const
	{
		assert(symmLocal_);
		SizeType ntensors = tensorSrep.size();
		for (SizeType i = 0; i < ntensors; ++i)
			if (!symmetriesPass(tensorSrep(i),summed,free))
				return false;

		return true;
	}

	bool symmetriesPass(const TensorStanza& ts,
	                    const VectorSizeType& summed,
	                    const VectorSizeType& free) const
	{
		assert(symmLocal_);
		PsimagLite::String tensorNameId = ts.name() + ttos(ts.id());
		SizeType tensorIndex = symmLocal_->nameIdToIndex(tensorNameId);
		if (ts.maxTag('f') == 0) return true;
		if (tensorIndex >= symmLocal_->size())
			assert(false);

		SizeType legs = ts.legs();
		SizeType ins = ts.ins();
		SizeType qin = 0;
		SizeType qout = 0;

		for (SizeType j = 0; j < legs; ++j) {
			SizeType index = ts.legTag(j);
			SizeType tmp = 0;

			switch (ts.legType(j)) {

			case TensorStanza::INDEX_TYPE_SUMMED:
				assert(index < summed.size());
				tmp = symmLocal_->q(tensorIndex,j)->operator[](summed[index]);
				break;

			case TensorStanza::INDEX_TYPE_FREE:
				assert(index < free.size());
				tmp = symmLocal_->q(tensorIndex,j)->operator[](free[index]);
				break;

			case  TensorStanza::INDEX_TYPE_DUMMY:
				break;

			default:
				PsimagLite::RuntimeError("symmetriesPass: Wrong index type\n");
			}

			if (j < ins)
				qin += tmp;
			else
				qout += tmp;
		}

		return (qin == qout);
	}

	SizeType idNameToIndex(PsimagLite::String name, SizeType id) const
	{
		typename MapPairStringSizeType::iterator it =
		        nameIdsTensor_.find(PairStringSizeType(name,id));
		if (it == nameIdsTensor_.end())
			throw PsimagLite::RuntimeError("idNameToIndex: key not found\n");
		return it->second;
	}

	void setQnsForOutput(VectorVectorSizeType& q)
	{
		// remap tensor indexing into symm local indexing
		PairStringSizeType p = tensorNameIds_[indexOfOutputTensor_];
		PsimagLite::String str = p.first + ttos(p.second);

		SizeType legs = srepEq_.lhs().legs();
		VectorSizeType v(legs, 0);
		for (SizeType j = 0; j < legs; ++j) {
			assert(srepEq_.lhs().legType(j) == TensorStanza::INDEX_TYPE_FREE);
			v[j] = srepEq_.lhs().legTag(j);
		}

		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType iperm(v.size(), 0);
		sort.sort(v, iperm);

		// set qs
		symmLocal_->addTensor(str, q, iperm);
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
	SymmetryLocalType* symmLocal_;
	bool modify_;
	SizeType indexOfOutputTensor_;
	VectorTensorType garbage_;
	mutable VectorVectorSizeType garbage2_;
};
}
#endif // MERA_TensorEvalSlow_H
