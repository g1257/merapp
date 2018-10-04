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
#include "TensorEvalSlow.h"
#include "TensorEvalNew.h"
#include <algorithm>
#include "Sort.h"
#include "Matrix.h"
#include "LanczosSolver.h"
#include "CrsMatrix.h"
#include "BLAS.h"
#include "SymmetryLocal.h"
#include "ParallelEnvironHelper.h"
#include "Parallelizer.h"
#include "ParametersForMera.h"
#include "Random48.h"

namespace Mera {

template<typename ComplexOrRealType, typename IoInType>
class TensorOptimizer {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<TensorStanza::IndexDirectionEnum>::Type VectorDirType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef ParallelEnvironHelper<ComplexOrRealType> ParallelEnvironHelperType;

public:

	typedef ParametersForMera<ComplexOrRealType> ParametersForMeraType;
	typedef TensorEvalBase<ComplexOrRealType> TensorEvalBaseType;
	typedef TensorEvalSlow<ComplexOrRealType> TensorEvalSlowType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename TensorEvalBaseType::PairStringSizeType PairStringSizeType;
	typedef typename TensorEvalBaseType::VectorPairStringSizeType VectorPairStringSizeType;
	typedef typename TensorEvalBaseType::TensorType TensorType;
	typedef typename TensorEvalBaseType::VectorTensorType VectorTensorType;
	typedef typename TensorEvalBaseType::SrepStatementType SrepStatementType;
	typedef typename PsimagLite::Vector<SrepStatementType*>::Type VectorSrepStatementType;
	typedef typename ParallelEnvironHelperType::MatrixType MatrixType;
	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef typename TensorEvalBaseType::MapPairStringSizeType MapPairStringSizeType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef PsimagLite::LanczosSolver<ParametersForSolverType,SparseMatrixType,VectorType>
	LanczosSolverType;
	typedef typename TensorEvalSlowType::SymmetryLocalType SymmetryLocalType;
	typedef typename PsimagLite::Stack<VectorType>::Type StackVectorType;

	TensorOptimizer(IoInType& io,
	                PsimagLite::String nameToOptimize,
	                SizeType idToOptimize,
	                const VectorPairStringSizeType& tensorNameAndIds,
	                MapPairStringSizeType& nameIdsTensor,
	                VectorTensorType& tensors,
	                const ParametersForSolverType& params,
	                const ParametersForMeraType& paramsForMera,
	                SymmetryLocalType* symmLocal)
	    : tensorToOptimize_(nameToOptimize,idToOptimize),
	      tensorNameIds_(tensorNameAndIds),
	      nameIdsTensor_(nameIdsTensor),
	      tensors_(tensors),
	      indToOptimize_(nameIdsTensor_[tensorToOptimize_]),
	      layer_(0),
	      firstOfLayer_(0),
	      indexOfRootTensor_(0),
	      params_(params),
	      paramsForMera_(paramsForMera),
	      symmLocal_(symmLocal),
	      verbose_(false),
	      rng_(1234)
	{
		io.readline(layer_,"Layer=");
		io.readline(firstOfLayer_,"FirstOfLayer=");

		io.readline(ignore_,"IgnoreTerm=");
		SizeType terms = 0;
		io.readline(terms,"Terms=");
		try {
			int tmp = 0;
			io.readline(tmp,"Verbose=");
			verbose_ = (tmp > 0);
		} catch (std::exception&) {}

		tensorSrep_.resize(terms,0);

		PsimagLite::String findStr = "Environ=";
		for (SizeType i = 0; i < terms; ++i) {
			PsimagLite::String srep;
			io.readline(srep,findStr);
			size_t index = srep.find("equal");
			if (index != PsimagLite::String::npos)
				srep.replace(index,5,"=");
			tensorSrep_[i] = new SrepStatementType(srep);
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

	void optimize(SizeType iters,
	              SizeType upIter,
	              PsimagLite::String evaluator)
	{
		if (tensorSrep_.size() == 0) return;

		saveTensor();

		SizeType ins = tensors_[indToOptimize_]->ins();
		SizeType outs = tensors_[indToOptimize_]->args() - ins;
		PsimagLite::String cond = conditionToSrep(tensorToOptimize_,ins,outs);

		SrepStatementType* condSrep = 0;
		if (cond != "") {
			condSrep = new SrepStatementType(cond);
			std::cout<<"cond="<<cond<<"\n";
		}

		RealType tolerance = paramsForMera_.tolerance;
		RealType eprev = 0.0;
		for (SizeType iter = 0; iter < iters; ++iter) {

			RealType e = optimizeInternal(iter, upIter, evaluator);
			if (tensorToOptimize_.first == "r") {
				std::cout<<"energy="<<e<<"\n";
				break;
			}

			if (iter > 0 && tolerance > 0 && fabs(eprev-e) < tolerance)
				break;
			std::cerr<<"TensorOptimizer: e="<<e;
			std::cerr<<" eprev="<<eprev<<" tolerance="<<tolerance<<"\n";
			eprev = e;

			if (!condSrep) continue;

			if (condSrep->lhs().maxTag('f') == 0) continue;

			ParallelEnvironHelperType parallelEnvironHelper(tensorSrep_,
			                                                evaluator,
			                                                ignore_,
			                                                tensorNameIds_,
			                                                nameIdsTensor_,
			                                                tensors_,
			                                                symmLocal_);

			MatrixType condMatrix;
			parallelEnvironHelper.appendToMatrix(condMatrix, *condSrep, evaluator);
			if (!isTheIdentity(condMatrix))
				std::cerr<<"not a isometry or unitary\n";
		}

		delete condSrep;
	}

	const PairStringSizeType& nameId() const { return tensorToOptimize_; }

	SizeType layer() const { return layer_; }

	static TensorEvalBaseType* getTensorEvalPtr(PsimagLite::String evaluator,
	                                            const SrepStatementType& srep,
	                                            VectorTensorType& tensors,
	                                            const VectorPairStringSizeType& tensorNameIds,
	                                            MapPairStringSizeType& nameIdsTensor,
	                                            SymmetryLocalType* symmLocal)
	{
		return ParallelEnvironHelperType::getTensorEvalPtr(evaluator,
		                                                   srep,
		                                                   tensors,
		                                                   tensorNameIds,
		                                                   nameIdsTensor,
		                                                   symmLocal);
	}

	void restoreTensor()
	{
		if (stack_.size() == 0)
			throw PsimagLite::RuntimeError("restoreTensor: stack is empty\n");
		tensors_[indToOptimize_]->data() = stack_.top();
		stack_.pop();
	}

	const SizeType& firstOfLayer() const { return firstOfLayer_; }

	void copyFirstOfLayer(PsimagLite::String name, SizeType id)
	{
		SizeType ind = nameIdsTensor_[PairStringSizeType(name, id)];
		tensors_[indToOptimize_]->data() = tensors_[ind]->data();
	}

private:

	void saveTensor()
	{
		stack_.push(tensors_[indToOptimize_]->data());
	}

	PsimagLite::String conditionToSrep(PairStringSizeType nameId,
	                                   SizeType ins,
	                                   SizeType outs) const
	{
		PsimagLite::String name = nameId.first;
		SizeType id = nameId.second;
		PsimagLite::String srep = name + ttos(id) + "(";
		PsimagLite::String args("");

		for (SizeType j = 0; j < ins; ++j) {
			srep += "s" + ttos(j);
			if (j < ins - 1) srep += ",";
		}

		if (outs > 0) srep += "|";

		for (SizeType j = 0; j < outs; ++j) {
			srep += "f" + ttos(j);
			args += "f" + ttos(j);
			if (j < outs - 1) {
				srep += ",";
				args += ",";
			}
		}

		srep += ")" + name + "*" + ttos(id) + "(";
		for (SizeType j = 0; j < ins; ++j) {
			srep += "s" + ttos(j);
			if (j < ins - 1) srep += ",";
		}

		if (outs > 0) {
			srep += "|";
			args += "|";
		}

		for (SizeType j = 0; j < outs; ++j) {
			srep += "f" + ttos(j + outs);
			args += "f" + ttos(j + outs);
			if (j < outs - 1) {
				srep += ",";
				args += ",";
			}
		}

		srep += ")";
		return (args == "") ? "" : "u1000(" + args + ")= " + srep;
	}

	RealType optimizeInternal(SizeType iter, SizeType upIter, PsimagLite::String evaluator)
	{
		if (verbose_)
			std::cerr<<"ignore="<<ignore_<<"\n";
		typedef PsimagLite::Parallelizer<ParallelEnvironHelperType> ParallelizerType;
		ParallelizerType threadedEnviron(PsimagLite::Concurrency::codeSectionParams);

		ParallelEnvironHelperType parallelEnvironHelper(tensorSrep_,
		                                                evaluator,
		                                                ignore_,
		                                                tensorNameIds_,
		                                                nameIdsTensor_,
		                                                tensors_,
		                                                symmLocal_);

		threadedEnviron.loopCreate(parallelEnvironHelper);
		parallelEnvironHelper.sync();

		MatrixType m = parallelEnvironHelper.matrix();
		MatrixType mSrc = m;
		VectorRealType s(m.n_row(),0);
		if (tensorToOptimize_.first == "r") { // diagonalize
			std::cout<<"MATRIX_MAY_FOLLOW\n";
			if (!isHermitian(m,true)) {
				if (m.n_row() < 512) std::cout<<m;
				throw PsimagLite::RuntimeError("Not Hermitian H\n");
			}

			bool printmatrix = (params_.options.find("printmatrix") != PsimagLite::String::npos);
			if (params_.options.find("printMatrix") != PsimagLite::String::npos)
				printmatrix = true;
			if (printmatrix)
				if (m.n_row() < 512) std::cout<<m;

			SizeType args = tensors_[indToOptimize_]->args();
			if ((args & 1) || args < 2)
				throw PsimagLite::RuntimeError("r tensor should have even lengs\n");
			SizeType argsOver2 = args/2;
			SizeType rows = tensors_[indToOptimize_]->argSize(0);
			for (SizeType i = 1; i < argsOver2; ++i)
				rows *= tensors_[indToOptimize_]->argSize(i);

			SizeType cols = tensors_[indToOptimize_]->argSize(argsOver2);
			for (SizeType i = argsOver2 + 1; i < args; ++i)
				cols *= tensors_[indToOptimize_]->argSize(i);

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
			if (verbose_) std::cerr<<"r*Y(r)r="<<e<<"\n";

			bool stopEarly = (params_.options.find("stopEarly") != PsimagLite::String::npos);
			if (params_.options.find("stopearly") != PsimagLite::String::npos)
				stopEarly = true;
			if (stopEarly)
				throw PsimagLite::RuntimeError("stopEarly requested by user\n");

			return s[0];
		}

		RealType tmp = PsimagLite::norm2(m);
		if (verbose_)
			std::cerr<<"About to do svd matrix with norm2= "<<tmp<<"\n";
		MatrixType vt;
		svd('S',m,s,vt);
		page14StepL3(m,vt);
		RealType result = 0.0;
		for (SizeType i = 0; i < s.size(); ++i)
			result += s[i];
		std::cerr<<"ITER="<<iter<<" TensorOptimizer[";
		std::cerr<<indToOptimize_<<"] svdSumOfS= "<<result<<"\n";
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
		psimag::BLAS::GEMM('C','C',r1,cols,cols,-1.0,&(vt(0,0)),cols,&(m(0,0)),r2,0.0,&(t(0,0)),r1);

#ifndef DNDEBUG
		if (r1 == r2) {
			MatrixType tmp(r1,r2);
			psimag::BLAS::GEMM('C','N',r1,r1,r1,1.0,&(t(0,0)),r1,&(t(0,0)),r1,0.0,&(tmp(0,0)),r1);
			if (!isTheIdentity(tmp))
				throw PsimagLite::RuntimeError("internal error\n");
		}
#endif

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

		VectorType initialV(n);
		randomizeVector(initialV, 1, 0, rng_);
		lanczosSolver.computeOneState(s[0], gsVector, initialV, 0);
		for (SizeType i = 0; i < rows; ++i)
			for (SizeType j = 0; j < cols; ++j)
				t(i,j) = gsVector[i + j*rows];
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
	VectorSrepStatementType tensorSrep_;
	const VectorPairStringSizeType& tensorNameIds_;
	MapPairStringSizeType& nameIdsTensor_;
	VectorTensorType& tensors_;
	SizeType indToOptimize_;
	SizeType ignore_;
	SizeType layer_;
	SizeType firstOfLayer_;
	SizeType indexOfRootTensor_;
	const ParametersForSolverType& params_;
	const ParametersForMeraType& paramsForMera_;
	SymmetryLocalType* symmLocal_;
	bool verbose_;
	PsimagLite::Random48<double> rng_;
	StackVectorType stack_;
}; // class TensorOptimizer
} // namespace Mera
#endif // TENSOROPTIMIZER_H
