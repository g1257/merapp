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
#include "TensorBreakup.h"
#include "TensorEvalBase.h"
#include "SymmetryLocal.h"
#include "BLAS.h"
#include "PsimagLite.h"
#include "NameToIndexLut.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorEvalSlow : public TensorEvalBase<ComplexOrRealType> {

	typedef TensorSrep TensorSrepType;

public:

	typedef TensorEvalBase<ComplexOrRealType> TensorEvalBaseType;
	typedef typename TensorEvalBaseType::SrepStatementType SrepStatementType;
	typedef typename TensorEvalBaseType::HandleType HandleType;
	typedef typename TensorEvalBaseType::TensorType TensorType;
	typedef typename TensorEvalBaseType::VectorTensorType VectorTensorType;
	typedef typename TensorEvalBaseType::VectorSizeType VectorSizeType;
	typedef typename TensorEvalBaseType::PairStringSizeType PairStringSizeType;
	typedef typename TensorEvalBaseType::VectorPairStringSizeType VectorPairStringSizeType;
	typedef typename TensorEvalBaseType::MapPairStringSizeType MapPairStringSizeType;
	typedef typename PsimagLite::Vector<SrepStatementType*>::Type VectorSrepStatementType;
	typedef TensorBreakup::VectorStringType VectorStringType;
	typedef typename TensorType::MatrixType MatrixType;
	typedef SymmetryLocal SymmetryLocalType;
	typedef SymmetryLocalType::VectorVectorSizeType VectorVectorSizeType;

	static const SizeType EVAL_BREAKUP = TensorBreakup::EVAL_BREAKUP;

	TensorEvalSlow(const SrepStatementType& tSrep,
	               const VectorTensorType& vt,
	               NameToIndexLut<TensorType>& nameToIndexLUT,
	               SymmetryLocalType* symmLocal,
	               bool modify = EVAL_BREAKUP)
	    : srepStatement_(tSrep),
	      data_(vt), // deep copy
	      nameToIndexLUT_(nameToIndexLUT),
	      symmLocal_(symmLocal),
	      modify_(modify)
	{
		if (!modify_) return;

		TensorBreakup tensorBreakup(srepStatement_.lhs(), srepStatement_.rhs());
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
				//std::cout<<"Definition of "<<srepStatement_.rhs().sRep()<<" is ";
				//std::cout<<vstr[i + 1]<<"\n";
				srepStatement_.rhs() = TensorSrep(vstr[i + 1]);
				outputLocation = i;
			}

			if (temporaryName[0] != 't') continue;

			// add this temporary, call setDimensions for output tensor later
			TensorStanza tmpStanza(vstr[i]);
			VectorSizeType args(1,1); // bogus
			const PsimagLite::String nameAndId = upToParens(temporaryName);
			TensorType* t = new TensorType(nameAndId, args, tmpStanza.ins());
			garbage_.push_back(t);
			data_.push_back(t);
			nameToIndexLUT.push(nameAndId);
		}

		VectorSrepStatementType veqs;
		TensorSrepType::VectorPairSizeType empty;
		for (SizeType i = 0; i < vstr.size(); i += 2) {
			veqs.push_back(new SrepStatementType(vstr[i] + "=" + vstr[i+1]));
			SizeType j = veqs.size() - 1;
			if (i != outputLocation)
				veqs[j]->canonicalize();
			veqs[j]->rhs().simplify(empty);

			TensorEvalSlow tEval(*(veqs[j]),
			                     data_,
			                     nameToIndexLUT_,
			                     symmLocal_,
			                     false);

			//std::cerr<<"Evaluation of "<<veqs[j]->sRep()<<"\n";
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
	}

	SizeType nameToIndexLut(PsimagLite::String name)
	{
		return nameToIndexLUT_(name);
	}

	HandleType operator()()
	{
		HandleType handle(HandleType::STATUS_DONE);
		if (modify_ && EVAL_BREAKUP) return handle;

		SizeType total = srepStatement_.lhs().maxTag('f') + 1;

		//if (srepStatement_.rhs().size() == 2 && total > 0 && !symmLocal_)
		//	return operatorParensFast();

		VectorSizeType dimensions(total, 0);
		VectorVectorSizeType q(total, 0);

		SizeType indexOfOutputTensor = nameToIndexLUT_(nameOfOutputTensor_);
		bool hasFree = srepStatement_.lhs().hasLegType('f');
		if (hasFree) {
			prepare(dimensions,q,srepStatement_.rhs(),TensorStanza::INDEX_TYPE_FREE);
			setQnsForOutput(q, indexOfOutputTensor);
		} else {
			assert(dimensions.size() == 1);
			dimensions[0] = 1;
		}

		VectorSizeType free(total, 0);

		if (dimensions.size() == 1 && dimensions[0] == 0)
			dimensions[0] = 1;
		outputTensor(indexOfOutputTensor).setSizes(dimensions);

		do {
			outputTensor(indexOfOutputTensor)(free) = slowEvaluator(free,srepStatement_.rhs());
		} while (ProgramGlobals::nextIndex(free,dimensions,total));

		return handle;
	}

	void printResult(std::ostream& os) const
	{
		SizeType indexOfOutputTensor = nameToIndexLUT_(nameOfOutputTensor_);
		SizeType total = outputTensor(indexOfOutputTensor).args();
		static VectorSizeType dimensions;
		if (total > dimensions.size()) dimensions.resize(total,0);

		for (SizeType i = 0; i < total; ++i)
			dimensions[i] = outputTensor(indexOfOutputTensor).argSize(i);

		static VectorSizeType free;
		if (total > free.size()) free.resize(total,0);
		else std::fill(free.begin(), free.end(), 0);

		do {
			SizeType index = outputTensor(indexOfOutputTensor).index(free);
			std::cout<<index<<" "<<outputTensor(indexOfOutputTensor)(free)<<"\n";
		} while (ProgramGlobals::nextIndex(free,dimensions,total));
	}

private:

	ComplexOrRealType slowEvaluator(const VectorSizeType& free,
	                                const TensorSrepType& srep)
	{
		SizeType total = srep.maxTag('s') + 1;
		VectorSizeType summed(total, 0);

		VectorSizeType dimensions(total, 0);

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
		SizeType mid = nameToIndexLUT_(stanza.name() + ttos(id));
		SizeType tensorIndex = (symmLocal_) ?
		            symmLocal_->nameIdToIndex(stanza.name() + ttos(id)) : 0;
		if (symmLocal_ && tensorIndex >= symmLocal_->size())
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
			if (symmLocal_ && type == TensorStanza::INDEX_TYPE_FREE) {
				const VectorSizeType* qSrc = symmLocal_->q(tensorIndex, j);
				assert(qSrc);
				assert(sIndex < q.size());
				VectorSizeType* ptr = new VectorSizeType(qSrc->size());
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

	ComplexOrRealType& evalThisTensor(const TensorStanza& ts,
	                                  const VectorSizeType& summed,
	                                  const VectorSizeType& free) const
	{
		SizeType id = ts.id();
		SizeType mid = nameToIndexLUT_(ts.name() + ttos(id));
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
		for (SizeType i = 0; i < ntensors; ++i) {
			// tensor r (root tensor) has no out legs, so different symmetry
			// other tensors might have different symmetry also
			// Therefore, symmetry as implemented only applies to u and w and h
			PsimagLite::String name = tensorSrep(i).name();
			if (name != "u" && name != "w" && name != "h")
				continue;
			if (!symmetriesPass(tensorSrep(i),summed,free))
				return false;
		}

		return true;
	}

	bool symmetriesPass(const TensorStanza& ts,
	                    const VectorSizeType& summed,
	                    const VectorSizeType& free) const
	{
		assert(symmLocal_);
		PsimagLite::String tensorNameId = ts.name() + ttos(ts.id());
		SizeType tensorIndex = symmLocal_->nameIdToIndex(tensorNameId);
		//if (ts.maxTag('f') == 0) return true;
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

	HandleType operatorParensFast()
	{
		HandleType handle(HandleType::STATUS_DONE);
		assert(!modify_ || !EVAL_BREAKUP);
		assert(!symmLocal_);

		SizeType total = srepStatement_.lhs().maxTag('f') + 1;

		assert(srepStatement_.rhs().size() == 2 && total > 0);

		VectorSizeType dimensions(total, 0);
		VectorVectorSizeType q(total, 0);

		bool hasFree = srepStatement_.lhs().hasLegType('f');
		if (hasFree) {
			prepare(dimensions,q,srepStatement_.rhs(),TensorStanza::INDEX_TYPE_FREE);
			setQnsForOutput(q);
		} else {
			assert(dimensions.size() == 1);
			dimensions[0] = 1;
		}

		if (dimensions.size() == 1 && dimensions[0] == 0)
			dimensions[0] = 1;
		outputTensor().setSizes(dimensions);

		MatrixType m1;
		reshapeIntoMatrix(m1, srepStatement_.rhs()(0), dimensions);
		SizeType frees1 = m1.rows();
		SizeType summed = m1.cols();
		MatrixType m2;
		reshapeIntoMatrix(m2, srepStatement_.rhs()(1), dimensions);
		assert(summed == m2.rows());
		SizeType frees2 = m2.cols();
		assert(frees1*frees2 == volumeOf(dimensions));
		MatrixType m3(frees1, frees2);
		const ComplexOrRealType alpha = 1.0;
		const ComplexOrRealType beta = 0.0;
		psimag::BLAS::GEMM('N',
		                   'N',
		                   frees1,
		                   frees2,
		                   summed,
		                   alpha,
		                   &(m1(0,0)),
		                   frees1,
		                   &(m2(0,0)),
		                   summed,
		                   beta,
		                   &(m3(0,0)),
		                   m3.rows());

		reshapeIntoTensor(outputTensor(), m3, srepStatement_.rhs());

		return handle;
	}

	void reshapeIntoMatrix(MatrixType& m,
	                       const TensorStanza& ts,
	                       const VectorSizeType& dimensions) const
	{
		SizeType totalSummed = ts.maxTag('s') + 1;
		VectorSizeType dimensionsSummed(totalSummed, 0);
		VectorVectorSizeType q;
		bool hasSummed = ts.hasLegType('s');
		if (hasSummed) {
			prepareStanza(dimensionsSummed, q, ts, TensorStanza::INDEX_TYPE_SUMMED);
		} else {
			assert(dimensionsSummed.size() == 1);
			dimensionsSummed[0] = 1;
		}

		VectorSizeType reducedDimensions;
		reduceDimensions(reducedDimensions, dimensions, ts);
		m.resize(volumeOf(reducedDimensions), volumeOf(dimensionsSummed));
		SizeType total = reducedDimensions.size();
		VectorSizeType free(total, 0);

		do {
			SizeType row = vectorToIndex(free, reducedDimensions);
			reshapeIntoMatrix(m, ts, free, row, dimensionsSummed);
		} while (ProgramGlobals::nextIndex(free,reducedDimensions,total));
	}

	void reshapeIntoMatrix(MatrixType& m,
	                       const TensorStanza& ts1,
	                       const VectorSizeType& free,
	                       SizeType row,
	                       const VectorSizeType& dimensions) const
	{
		TensorStanza ts = ts1;
		ts.canonicalize();
		SizeType total = ts.maxTag('s') + 1;
		VectorSizeType summed(total, 0);

		do {
			SizeType col = vectorToIndex(summed, dimensions);
			m(row, col) = evalThisTensor(ts, summed, free);
		} while (ProgramGlobals::nextIndex(summed,dimensions,total));
	}

	void reshapeIntoTensor(TensorType& tensor,
	                       const MatrixType& src,
	                       const TensorSrep& srep) const
	{
		SizeType rows = src.rows();
		SizeType cols = src.cols();
		SizeType total1 = srep(0).maxTag('f') + 1;
		SizeType total2 = srep(1).maxTag('f') + 1;
		VectorVectorSizeType q;
		VectorSizeType dimensions1;
		prepareStanza(dimensions1,q,srep(0),TensorStanza::INDEX_TYPE_FREE);
		VectorSizeType dimensions2;
		prepareStanza(dimensions2,q,srep(1),TensorStanza::INDEX_TYPE_FREE);
		MatrixType frees1(rows, total1);
		MatrixType frees2(cols, total2);
		for (SizeType i = 0; i < rows; ++i)
			getFrees(frees1, i, dimensions1);

		for (SizeType j = 0; j < cols; ++j)
			getFrees(frees2, j, dimensions2);

		VectorSizeType frees(total1 + total2, 0);

		for (SizeType i = 0; i < rows; ++i) {
			for (SizeType j = 0; j < cols; ++j) {
				combineFrees(frees, frees1, frees2, i, j);
				// update output tensor
				err("need code to update output tensor\n");
			}
		}
	}

	SizeType volumeOf(const VectorSizeType& v) const
	{
		SizeType n = v.size();
		if (n == 0) return 1;
		SizeType prod = v[0];
		for (SizeType i = 1; i < n; ++i) {
			prod *= v[i];
		}

		return prod;
	}

	void reduceDimensions(VectorSizeType& rd,
	                      const VectorSizeType& v,
	                      const TensorStanza& ts) const
	{
		SizeType n = v.size();
		if (n == 0) return;
		PsimagLite::Vector<bool>::Type free(n, false);
		SizeType legs = ts.legs();
		SizeType count = 0;
		for (SizeType i = 0; i < legs; ++i) {
			if (ts.legType(i) != TensorStanza::INDEX_TYPE_FREE)
				continue;
			free[ts.legTag(i)] = true;
			++count;
		}

		rd.resize(count);
		SizeType j = 0;
		for (SizeType i = 0; i < n; ++i) {
			if (!free[i]) continue;
			rd[j++] = v[i];
		}
	}

	void getFrees(MatrixType& frees, SizeType row, const VectorSizeType& d) const
	{
		SizeType n = d.size();
		assert(n > 0);
		assert(frees.rows() == n);
		SizeType prod = d[0];
		for (SizeType i = 1; i < n - 1; ++i) {
			prod *= d[i-1];
		}

		SizeType temp = row;
		SizeType j = 0;
		SizeType k = (n >= 2) ? n - 2 : 0;
		while (j < n && k < n) {
			div_t q = div(temp, prod);
			frees(row, j++) = q.quot;
			temp = q.rem;
			prod /= d[k--];
		}

		frees(row, j++) = temp;
	}

	void combineFrees(VectorSizeType& frees,
	                  const MatrixType& frees1,
	                  const MatrixType& frees2,
	                  SizeType ind,
	                  SizeType jnd) const
	{
		SizeType total1 = frees1.cols();
		SizeType total2 = frees2.cols();
		assert(frees.size() == total1 + total2);
		for (SizeType k = 0; k < total1; ++k)
			frees[k] = frees1(ind, k);
		for (SizeType k = 0; k < total2; ++k)
			frees[k] = frees2(jnd, k);
	}

	SizeType vectorToIndex(const VectorSizeType& v,
	                       const VectorSizeType& d) const
	{
		assert(v.size() > 0);
		assert(v.size() == d.size());
		SizeType sum = v[0];
		for (SizeType i = 1; i < d.size(); ++i) {
			sum += v[i]*d[i-1];
		}

		return sum;
	}

	void setQnsForOutput(VectorVectorSizeType& q,
	                     SizeType indexOfOutputTensor)
	{
		// remap tensor indexing into symm local indexing
		PsimagLite::String str = data_[indexOfOutputTensor]->name();

		SizeType legs = srepStatement_.lhs().legs();
		VectorSizeType v(legs, 0);
		for (SizeType j = 0; j < legs; ++j) {
			assert(srepStatement_.lhs().legType(j) == TensorStanza::INDEX_TYPE_FREE);
			v[j] = srepStatement_.lhs().legTag(j);
		}

		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType iperm(v.size(), 0);
		sort.sort(v, iperm);

		// set qs
		if (symmLocal_)
			symmLocal_->addTensor(str, q, iperm);
	}

	TensorType& outputTensor(SizeType indexOfOutputTensor)
	{
		assert(indexOfOutputTensor < data_.size());
		return *(data_[indexOfOutputTensor]);
	}

	const TensorType& outputTensor(SizeType indexOfOutputTensor) const
	{
		assert(indexOfOutputTensor < data_.size());
		return *(data_[indexOfOutputTensor]);
	}

	static PsimagLite::String upToParens(PsimagLite::String str)
	{
		const SizeType l = str.length();
		PsimagLite::String buffer = "";
		for (SizeType i = 0; i < l; ++i) {
			if (str[i] == '(') return buffer;
			buffer += str[i];
		}

		return buffer;
	}

	TensorEvalSlow(const TensorEvalSlow& other);

	TensorEvalSlow& operator=(const TensorEvalSlow& other);

	SrepStatementType srepStatement_;
	VectorTensorType data_;
	NameToIndexLut<TensorType>& nameToIndexLUT_;
	SymmetryLocalType* symmLocal_;
	bool modify_;
	VectorTensorType garbage_;
};
}
#endif // MERA_TensorEvalSlow_H
