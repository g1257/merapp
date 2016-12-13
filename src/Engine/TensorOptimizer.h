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
#ifndef TENSOROPTIMIZER_H
#define TENSOROPTIMIZER_H
#include "Vector.h"
#include "TensorSrep.h"
#include "TensorEval.h"
#include <algorithm>
#include "Sort.h"
#include "Matrix.h"
#include "LanczosSolver.h"
#include "CrsMatrix.h"

namespace Mera {

template<typename ComplexOrRealType, typename IoInType>
class TensorOptimizer {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<TensorStanza::IndexDirectionEnum>::Type VectorDirType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;

public:

	typedef Mera::TensorEval<ComplexOrRealType> TensorEvalType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename TensorEvalType::PairStringSizeType PairStringSizeType;
	typedef typename TensorEvalType::VectorPairStringSizeType VectorPairStringSizeType;
	typedef typename TensorEvalType::TensorType TensorType;
	typedef typename TensorEvalType::VectorTensorType VectorTensorType;
	typedef typename TensorEvalType::SrepEquationType SrepEquationType;
	typedef typename PsimagLite::Vector<SrepEquationType*>::Type VectorSrepEquationType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef typename TensorEvalType::MapPairStringSizeType MapPairStringSizeType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef PsimagLite::LanczosSolver<ParametersForSolverType,
	SparseMatrixType,
	VectorType> LanczosSolverType;
	TensorOptimizer(IoInType& io,
	                PsimagLite::String nameToOptimize,
	                SizeType idToOptimize,
	                const VectorPairStringSizeType& tensorNameAndIds,
	                MapPairStringSizeType& nameIdsTensor,
	                VectorTensorType& tensors,
	                const ParametersForSolverType& params)
	    : tensorToOptimize_(nameToOptimize,idToOptimize),
	      tensorNameIds_(tensorNameAndIds),
	      nameIdsTensor_(nameIdsTensor),
	      tensors_(tensors),
	      indToOptimize_(nameIdsTensor_[tensorToOptimize_]),
	      layer_(0),
	      indexOfRootTensor_(0),
	      params_(params)
	{
		io.readline(layer_,"Layer=");
		io.readline(ignore_,"IgnoreTerm=");
		SizeType terms = 0;
		io.readline(terms,"Terms=");
		tensorSrep_.resize(terms,0);

		PsimagLite::String findStr = "Environ=";
		for (SizeType i = 0; i < terms; ++i) {
			PsimagLite::String srep;
			io.readline(srep,findStr);
			size_t index = srep.find("equal");
			if (index != PsimagLite::String::npos)
				srep.replace(index,5,"=");
			tensorSrep_[i] = new SrepEquationType(srep,
			                                      tensors,
			                                      tensorNameAndIds,
			                                      nameIdsTensor);
		}

		bool flag = false;
		for (SizeType i = 0; i < tensorNameIds_.size(); ++i) {
			if (tensorNameIds_[i].first == "r") {
				indexOfRootTensor_ = i;
				flag = true;
				break;
			}
		}

		if (!flag)
			throw PsimagLite::RuntimeError("TensorOptimizer: Root tensor not found\n");
	}

	~TensorOptimizer()
	{
		SizeType terms = tensorSrep_.size();
		for (SizeType i = 0; i < terms; ++i) {
			delete tensorSrep_[i];
			tensorSrep_[i] = 0;
		}
	}

	void optimize(SizeType iters, SizeType upIter)
	{
		SizeType ins = tensors_[indToOptimize_]->ins();
		SizeType outs = tensors_[indToOptimize_]->args() - ins;
		PsimagLite::String cond = conditionToSrep(tensorToOptimize_,ins,outs);
		std::cout<<"cond="<<cond<<"\n";
		TensorSrep condSrep(cond);

		RealType eprev = 0.0;
		for (SizeType iter = 0; iter < iters; ++iter) {
			RealType e = optimizeInternal(iter, upIter);
			if (tensorToOptimize_.first == "r") {
				std::cout<<"energy="<<e<<"\n";
			}

			if (iter > 0 && fabs(eprev-e)<1e-4)
				break;
			eprev = e;

			if (condSrep.maxTag('f') == 0) continue;

			//			MatrixType condMatrix;
			//			appendToMatrix(condMatrix,condSrep);
			//			assert(isTheIdentity(condMatrix));
		}
	}

	const PairSizeType& nameId() const { return tensorToOptimize_; }

	SizeType layer() const { return layer_; }

private:

	PsimagLite::String conditionToSrep(PairStringSizeType nameId,
	                                   SizeType ins,
	                                   SizeType outs) const
	{
		PsimagLite::String name = nameId.first;
		SizeType id = nameId.second;
		PsimagLite::String srep = name + ttos(id) + "(";
		for (SizeType j = 0; j < ins; ++j) {
			srep += "s" + ttos(j);
			if (j < ins - 1) srep += ",";
		}

		if (outs > 0) srep += "|";
		for (SizeType j = 0; j < outs; ++j) {
			srep += "f" + ttos(j);
			if (j < outs - 1) srep += ",";
		}

		srep += ")" + name + "*" + ttos(id) + "(";
		for (SizeType j = 0; j < ins; ++j) {
			srep += "s" + ttos(j);
			if (j < ins - 1) srep += ",";
		}

		if (outs > 0) srep += "|";
		for (SizeType j = 0; j < outs; ++j) {
			srep += "f" + ttos(j + outs);
			if (j < outs - 1) srep += ",";
		}

		srep += ")";
		return srep;
	}

	RealType optimizeInternal(SizeType iter, SizeType upIter)
	{
		SizeType terms = tensorSrep_.size();
		MatrixType m;
		std::cerr<<"ignore="<<ignore_<<"\n";
		for (SizeType i = 0; i < terms; ++i) {
			if (i == ignore_) continue;
			appendToMatrix(m,*(tensorSrep_[i]));
		}

		MatrixType mSrc = m;
		VectorRealType s(m.n_row(),0);
		if (tensorToOptimize_.first == "r") { // diagonalize
			if (!isHermitian(m,true))
				throw PsimagLite::RuntimeError("Not Hermitian H\n");

			SizeType rows = tensors_[indToOptimize_]->argSize(0);
			SizeType cols = tensors_[indToOptimize_]->argSize(1);
			assert(rows*cols == m.n_row());
			MatrixType t(rows,cols);
			bool fullDiag = (params_.options.find("fulldiag") != PsimagLite::String::npos);
			if (fullDiag) {
				diag(m,s,'V');
				topTensorFoldVector(t,m);
			} else {
				lanczosDiag(t,s,m);
			}

			tensors_[indToOptimize_]->setToMatrix(t);
			assert(0 < s.size());

			RealType e = computeRyR(mSrc);
			std::cerr<<"r*Y(r)r="<<e<<"\n";

			return s[0];
		}

		RealType tmp = PsimagLite::norm2(m);
		std::cerr<<"About to do svd matrix with norm2= "<<tmp<<"\n";
		MatrixType vt;
		svd('S',m,s,vt);
		page14StepL3(m,vt);
		RealType result = 0.0;
		for (SizeType i = 0; i < s.size(); ++i)
			result += s[i];
		std::cerr<<"TensorOptimizer["<<indToOptimize_<<"] svdSumOfS= "<<result<<"\n";
		return result;
	}

	void page14StepL3(const MatrixType& m,
	                  const MatrixType& vt)
	{
		SizeType r1 = vt.n_col();
		SizeType r2 = m.n_row();
		SizeType cols = m.n_col();
		assert(cols == vt.n_row());
		MatrixType t(r1,r2);
		for (SizeType i = 0; i < r1; ++i) {
			for (SizeType j = 0; j < r2; ++j) {
				ComplexOrRealType sum = 0.0;
				for (SizeType k = 0; k < cols; ++k)
					sum -= vt(k,i)*m(j,k);
				t(i,j) = PsimagLite::conj(sum);
			}
		}

		tensors_[indToOptimize_]->setToMatrix(t);
	}

	void topTensorFoldVector(MatrixType& t,
	                         const MatrixType& eigenvector) const
	{
		assert(tensorToOptimize_.first == "r");
		SizeType rows = t.n_row();
		SizeType cols = t.n_col();
		for (SizeType i = 0; i < rows; ++i)
			for (SizeType j = 0; j < cols; ++j)
				t(i,j) = eigenvector(i + j*rows,0);
	}

	void lanczosDiag(MatrixType& t, VectorRealType& s, const MatrixType& src) const
	{
		SizeType n = src.n_row();
		assert(tensorToOptimize_.first == "r");
		SizeType rows = t.n_row();
		SizeType cols = t.n_col();
		assert(n == rows*cols);
		t.resize(rows,cols);
		SparseMatrixType srcSparse(src);
		LanczosSolverType lanczosSolver(srcSparse,params_);
		VectorType gsVector(n,0.0);
		if (s.size() == 0) s.resize(1,0.0);
		lanczosSolver.computeGroundState(s[0],gsVector);
		for (SizeType i = 0; i < rows; ++i)
			for (SizeType j = 0; j < cols; ++j)
				t(i,j) = gsVector[i + j*rows];
	}

	void appendToMatrix(MatrixType& m, SrepEquationType& eq) const
	{
		SizeType total = eq.rhs().maxTag('f') + 1;
		VectorSizeType freeIndices(total,0);
		VectorDirType directions(total,TensorStanza::INDEX_DIR_IN);
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

		// prepare output tensor for evaluator
		eq.outputTensor().setSizes(dimensions);

		// evaluate environment
		TensorEvalType tensorEval(eq,tensors_,tensorNameIds_,nameIdsTensor_);
		typename TensorEvalType::HandleType handle = tensorEval();
		while (!handle.done());

		// copy result into m
		SizeType count = 0;
		do {
			PairSizeType rc = getRowAndColFromFree(freeIndices,dimensions,directions);
			ComplexOrRealType tmp = eq.outputTensor()(freeIndices);
			//std::cerr<<"MATRIX i=" <<rc.first<<" j="<<rc.second<<"\n";
			m(rc.first,rc.second) += tmp;
			count++;
		} while (TensorEvalType::nextIndex(freeIndices,dimensions));
		//std::cerr<<count<<"\n";
	}

	void prepareFreeIndices(VectorDirType& directions,
	                        VectorBoolType& conjugate,
	                        VectorSizeType& dimensions,
	                        const TensorSrep& t) const
	{
		assert(dimensions.size() == directions.size());
		assert(dimensions.size() == conjugate.size());

		SizeType ntensors = t.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			PsimagLite::String name = t(i).name();
			SizeType id = t(i).id();
			PairStringSizeType p(name,id);

			SizeType ind = nameIdsTensor_[p];
			assert(ind < tensors_.size());

			SizeType ins = t(i).ins();
			SizeType outs = t(i).outs();
			bool conjugate1 = t(i).isConjugate();
			for (SizeType j = 0; j < ins; ++j) {
				TensorStanza::IndexTypeEnum legType = t(i).legType(j,
				                                                   TensorStanza::INDEX_DIR_IN);
				if (legType != TensorStanza::INDEX_TYPE_FREE) continue;
				SizeType index = t(i).legTag(j,
				                             TensorStanza::INDEX_DIR_IN);
				assert(index < dimensions.size());
				dimensions[index] = tensors_[ind]->dimension(j);
				assert(index < directions.size());
				directions[index] = TensorStanza::INDEX_DIR_IN;
				conjugate[index] = conjugate1;
			}

			for (SizeType j = 0; j < outs; ++j) {
				TensorStanza::IndexTypeEnum legType = t(i).legType(j,
				                                                   TensorStanza::INDEX_DIR_OUT);
				if (legType != TensorStanza::INDEX_TYPE_FREE) continue;
				SizeType index = t(i).legTag(j,
				                             TensorStanza::INDEX_DIR_OUT);
				assert(index < dimensions.size());
				dimensions[index] = tensors_[ind]->dimension(j+ins);
				assert(index < directions.size());
				directions[index] = TensorStanza::INDEX_DIR_OUT;
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

	RealType computeRyR(const MatrixType& y) const
	{
		if (indToOptimize_ != indexOfRootTensor_)
			throw PsimagLite::RuntimeError("Don't call computeRyR unless optimizing root\n");

		RealType sum = 0.0;
		const TensorType& r = *tensors_[indToOptimize_];
		assert(r.args() == 2);
		SizeType rows = r.argSize(0);
		SizeType cols = r.argSize(1);
		assert(rows*cols == y.n_col());
		assert(y.n_row() == y.n_col());
		VectorSizeType args1(2,0);
		VectorSizeType args2(2,0);
		for (SizeType i = 0; i < rows; ++i) {
			args1[0] = i;
			for (SizeType j = 0; j < cols; ++j) {
				args1[1] = j;
				for (SizeType k = 0; k < rows; ++k) {
					args2[0] = k;
					for (SizeType l = 0; l < cols; ++l) {
						args2[1] = l;
						sum += PsimagLite::conj(r(args1))*y(i+j*rows,k+l*rows)*r(args2);
					}
				}
			}
		}

		return sum;
	}

	TensorOptimizer(const TensorOptimizer&);

	TensorOptimizer& operator=(const TensorOptimizer&);

	PairStringSizeType tensorToOptimize_;
	VectorSrepEquationType tensorSrep_;
	const VectorPairStringSizeType& tensorNameIds_;
	MapPairStringSizeType& nameIdsTensor_;
	VectorTensorType& tensors_;
	SizeType indToOptimize_;
	SizeType ignore_;
	SizeType layer_;
	SizeType indexOfRootTensor_;
	const ParametersForSolverType& params_;
}; // class TensorOptimizer
} // namespace Mera
#endif // TENSOROPTIMIZER_H
