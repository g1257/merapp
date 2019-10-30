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
#ifndef PARALLELENVIRONHELPER_H
#define PARALLELENVIRONHELPER_H
#include "Matrix.h"
#include "TensorEvalSlow.h"
#include "TensorEvalBase.h"
#include "TensorEvalNew.h"
#include "Vector.h"
#include "TensorStanza.h"

namespace  Mera {

template<typename ComplexOrRealType>
class ParallelEnvironHelper {

public:

	typedef TensorEvalBase<ComplexOrRealType> TensorEvalBaseType;
	typedef TensorEvalSlow<ComplexOrRealType> TensorEvalSlowType;
	typedef TensorEvalNew<ComplexOrRealType> TensorEvalNewType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename TensorEvalBaseType::PairStringSizeType PairStringSizeType;
	typedef typename TensorEvalBaseType::TensorType TensorType;
	typedef typename TensorEvalBaseType::VectorTensorType VectorTensorType;
	typedef typename TensorEvalBaseType::SrepStatementType SrepStatementType;
	typedef typename PsimagLite::Vector<SrepStatementType*>::Type VectorSrepStatementType;
	typedef PsimagLite::Vector<TensorStanza::IndexDirectionEnum>::Type VectorDirType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<MatrixType*>::Type VectorMatrixType;
	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef typename TensorEvalBaseType::MapPairStringSizeType MapPairStringSizeType;
	typedef typename TensorEvalBaseType::VectorPairStringSizeType VectorPairStringSizeType;
	typedef typename TensorEvalSlowType::SymmetryLocalType SymmetryLocalType;

	ParallelEnvironHelper(VectorSrepStatementType& tensorSrep,
	                      PsimagLite::String evaluator,
	                      SizeType ignore,
	                      VectorTensorType& tensors,
	                      NameToIndexLut<TensorType>& nameToIndexLut,
	                      SymmetryLocalType* symmLocal)
	    : tensorSrep_(tensorSrep),
	      evaluator_(evaluator),
	      ignore_(ignore),
	      tensors_(tensors),
	      nameToIndexLut_(nameToIndexLut),
	      symmLocal_(symmLocal),
	      m_(PsimagLite::Concurrency::codeSectionParams.npthreads, 0)
	{
		for (SizeType i = 0; i < m_.size(); ++i)
			m_[i] = new MatrixType;
	}

	~ParallelEnvironHelper()
	{
		for (SizeType i = 0; i < m_.size(); ++i) {
			delete m_[i];
			m_[i] = 0;
		}
	}

	void doTask(SizeType taskNumber, SizeType threadNum)
	{
		if (taskNumber == ignore_) return;
		appendToMatrix(*(m_[threadNum]), *(tensorSrep_[taskNumber]), evaluator_);
	}

	SizeType tasks() const { return tensorSrep_.size(); }

	const MatrixType& matrix() const
	{
		if (m_.size() == 0 || !(m_[0]))
			throw PsimagLite::RuntimeError("ParallelEnvironHelper::matrix() failed\n");
		return *(m_[0]);
	}

	void sync()
	{
		if (m_.size() == 0) return;
		if (!(m_[0])) return;
		for (SizeType i = 1; i < m_.size(); ++i)
			checkAndAccumulate(*(m_[i]));
	}

	void appendToMatrix(MatrixType& m,
	                    SrepStatementType& eq,
	                    PsimagLite::String evaluator)
	{
		SizeType total = eq.rhs().maxTag('f') + 1;
		VectorSizeType freeIndices(total,0);
		VectorDirType directions(total,TensorStanza::TensorLegType::INDEX_DIR_IN);
		VectorSizeType dimensions(total,0);
		VectorBoolType conjugate(total,false);
		prepareFreeIndices(directions,conjugate,dimensions,eq.rhs());
		modifyDirections(directions,conjugate);
		PairSizeType rc = getRowsAndCols(dimensions,directions);
		if (m.n_row() == 0) {
			m.resize(rc.first, rc.second);
			m.setTo(0.0);
		} else if (m.n_row() != rc.first || m.n_col() != rc.second) {
			PsimagLite::String str("Hamiltonian terms environ \n");
			throw PsimagLite::RuntimeError(str);
		}

		assert(m.n_row() > 0 && m.n_col() > 0);

		// prepare output tensor for evaluator
		outputTensor(eq).setSizes(dimensions);

		// evaluate environment
		TensorEvalBaseType* tensorEval =  getTensorEvalPtr(evaluator,
		                                                   eq,
		                                                   tensors_,
		                                                   nameToIndexLut_,
		                                                   symmLocal_);

		typename TensorEvalBaseType::HandleType handle = tensorEval->operator()();
		while (!handle.done());

		delete tensorEval;
		tensorEval = 0;

		// copy result into m
		SizeType count = 0;
		do {
			PairSizeType rc = getRowAndColFromFree(freeIndices,dimensions,directions);
			ComplexOrRealType tmp = outputTensor(eq)(freeIndices);
			m(rc.first,rc.second) += tmp;
			count++;
		} while (ProgramGlobals::nextIndex(freeIndices,dimensions,total));
	}

	static TensorEvalBaseType* getTensorEvalPtr(PsimagLite::String evaluator,
	                                            const SrepStatementType& srep,
	                                            VectorTensorType& tensors,
	                                            NameToIndexLut<TensorType>& nameToIndexLut,
	                                            SymmetryLocalType* symmLocal)
	{
		TensorEvalBaseType* tensorEval = 0;
		if (evaluator == "slow") {
			tensorEval = new TensorEvalSlowType(srep, tensors, nameToIndexLut, symmLocal);
		} else if (evaluator == "new") {
			tensorEval = new TensorEvalNewType(srep, tensors);
		} else {
			throw PsimagLite::RuntimeError("Unknown evaluator " + evaluator + "\n");
		}

		return tensorEval;
	}

private:

	void prepareFreeIndices(VectorDirType& directions,
	                        VectorBoolType& conjugate,
	                        VectorSizeType& dimensions,
	                        const TensorSrep& t) const
	{
		assert(dimensions.size() == directions.size());
		assert(dimensions.size() == conjugate.size());

		SizeType ntensors = t.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType ind = nameToIndexLut_(t(i).fullName());
			assert(ind < tensors_.size());

			SizeType ins = t(i).ins();
			SizeType legs = t(i).legs();
			bool conjugate1 = t(i).isConjugate();
			for (SizeType j = 0; j < legs; ++j) {
				TensorStanza::IndexTypeEnum legType = t(i).legType(j);
				if (legType != TensorStanza::INDEX_TYPE_FREE) continue;
				SizeType index = t(i).legTag(j);
				assert(index < dimensions.size());
				dimensions[index] = tensors_[ind]->dimension(j);
				assert(index < directions.size());
				directions[index] = (j < ins) ? TensorStanza::INDEX_DIR_IN :
				                                TensorStanza::INDEX_DIR_OUT;
				conjugate[index] = conjugate1;
			}
		}
	}

	PairSizeType getRowAndColFromFree(VectorSizeType& freeIndices,
	                                  const VectorSizeType& dimensions,
	                                  const VectorDirType& dirs) const
	{
		SizeType n = freeIndices.size();
		assert(n == dirs.size());
		assert(n == dimensions.size());
		SizeType row = 0;
		SizeType col = 0;
		SizeType prodRow = 1;
		SizeType prodCol = 1;
		for (SizeType i = 0; i < n; ++i) {
			if (dirs[i] == TensorStanza::INDEX_DIR_IN) {
				row += freeIndices[i]*prodRow;
				prodRow *= dimensions[i];
			} else {
				col += freeIndices[i]*prodCol;
				prodCol *= dimensions[i];
			}
		}

		return PairSizeType(row,col);
	}

	PairSizeType getRowsAndCols(const VectorSizeType& dimensions,
	                            const VectorDirType& dirs) const
	{
		SizeType n = dimensions.size();
		assert(n == dirs.size());

		VectorSizeType freeIndices(n,0);
		for (SizeType i = 0; i < n; ++i) {
			if (dimensions[i] == 0) continue;
			freeIndices[i] = dimensions[i] - 1;
		}

		PairSizeType p = getRowAndColFromFree(freeIndices,dimensions,dirs);
		return PairSizeType(p.first + 1, p.second + 1);
	}

	void modifyDirections(VectorDirType& dirs,
	                      const VectorBoolType& conjugate) const
	{
		SizeType n = dirs.size();
		assert(n == conjugate.size());

		SizeType counter = 0;
		for (SizeType i = 0; i < n; ++i) {
			if (conjugate[i]) dirs[i] = oppositeDir(dirs[i]);
			if (dirs[i] == dirs[0]) counter++;
		}

		if (counter < n) return;
		counter /= 2;
		for (SizeType i = counter; i < n; ++i)
			dirs[i] = oppositeDir(dirs[0]);
	}

	TensorStanza::IndexDirectionEnum oppositeDir(TensorStanza::IndexDirectionEnum dir) const
	{
		TensorStanza::IndexDirectionEnum in = TensorStanza::INDEX_DIR_IN;
		TensorStanza::IndexDirectionEnum out = TensorStanza::INDEX_DIR_OUT;
		return (dir == in) ? out : in;
	}

	TensorType& outputTensor(const SrepStatementType& eq)
	{
		SizeType indexOfOutputTensor = nameToIndexLut_(eq.nameIdOfOutput());
		assert(indexOfOutputTensor < tensors_.size());
		return *(tensors_[indexOfOutputTensor]);
	}

	void checkAndAccumulate(const MatrixType& m) const
	{
		if (m.n_row() == 0 && m.n_col() == 0) return;
		assert(m_.size() > 0);
		SizeType nrow = m_[0]->n_row();
		SizeType ncol = m_[0]->n_col();

		if (m.n_row() != nrow || m.n_col() != ncol)
			throw PsimagLite::RuntimeError("checkAndAccumulate matrix m failed\n");

		for (SizeType i = 0; i < nrow; ++i)
			for (SizeType j = 0; j < ncol; ++j)
				m_[0]->operator()(i,j) += m(i,j);
	}

	VectorSrepStatementType& tensorSrep_;
	PsimagLite::String evaluator_;
	SizeType ignore_;
	VectorTensorType& tensors_;
	NameToIndexLut<TensorType>& nameToIndexLut_;
	SymmetryLocalType* symmLocal_;
	VectorMatrixType m_;
}; // class ParallelEnvironHelper
}

#endif // PARALLELENVIRONHELPER_H
