#ifndef TENSOROPTIMIZER_H
#define TENSOROPTIMIZER_H
#include "Vector.h"
#include "IoSimple.h"
#include "TensorSrep.h"
#include "TensorEval.h"
#include <algorithm>
#include "Sort.h"
#include "Matrix.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorOptimizer {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<TensorSrep*>::Type VectorTensorSrepType;
	typedef std::pair<PsimagLite::String,SizeType> PairStringSizeType;
	typedef PsimagLite::Vector<TensorStanza::IndexDirectionEnum>::Type VectorDirType;
	typedef PsimagLite::Vector<PairStringSizeType>::Type VectorPairStringSizeType;

public:

	typedef PsimagLite::IoSimple::In IoInType;
	typedef Mera::TensorEval<double> TensorEvalType;
	typedef TensorEvalType::TensorType TensorType;
	typedef TensorEvalType::VectorTensorType VectorTensorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef std::pair<SizeType,SizeType> PairSizeType;

	TensorOptimizer(IoInType& io,
	                PsimagLite::String dstr,
	                PsimagLite::String nameToOptimize,
	                SizeType idToOptimize)
	{
		initTensorSreps(io,nameToOptimize,idToOptimize);

		initTensorNameIds();

		initTensors(dstr);

		std::cerr<<"TensorOptimizer::ctor() done\n";
	}

	~TensorOptimizer()
	{
		SizeType terms = tensorSrep_.size();
		for (SizeType i = 0; i < terms; ++i) {
			delete tensorSrep_[i];
			tensorSrep_[i] = 0;
		}

		terms = tensors_.size();
		for (SizeType i = 0; i < terms; ++i) {
			delete tensors_[i];
			tensors_[i] = 0;
		}
	}

	void optimize()
	{
		SizeType terms = tensorSrep_.size();
		for (SizeType i = 0; i < terms; ++i) {
			appendToMatrix(*(tensorSrep_[i]));
			std::cerr<<"WARNING: DOING ONE TERM ONLY!!\n";
			break;
		}
	}

private:

	void initTensorSreps(IoInType& io,
	                     PsimagLite::String nameToOptimize,
	                     SizeType idToOptimize)
	{
		SizeType terms = 0;
		io.readline(terms,"TERMS=");
		std::cerr<<"Read "<<terms<<" for tensor id "<<idToOptimize<<"\n";
		tensorSrep_.resize(terms,0);
		for (SizeType i = 0; i < terms; ++i) {
			PsimagLite::String srep;
			io.readline(srep,"STRING=");
			std::cerr<<"Read string "<<srep<<"\n";
			tensorSrep_[i] = new TensorSrep(srep);
			findTensors(*(tensorSrep_[i]),nameToOptimize,idToOptimize);
			std::cerr<<"Free indices= "<<(1+tensorSrep_[i]->maxTag('f'))<<"\n";
		}
	}

	void initTensorNameIds()
	{
		PsimagLite::Sort<VectorPairStringSizeType> sort;
		VectorSizeType perm(tensorNameIds_.size(),0);
		sort.sort(tensorNameIds_,perm);
		SizeType end = (std::unique(tensorNameIds_.begin(),
		                            tensorNameIds_.end()) -
		                tensorNameIds_.begin());
		tensorNameIds_.resize(end);
	}

	void initTensors(PsimagLite::String dstr)
	{
		tensors_.resize(tensorNameIds_.size());
		SizeType ntensors = tensors_.size();

		TensorSrep td(dstr);
		if (td.size() != ntensors) {
			PsimagLite::String str("TensorOptimizer dimension string " + ttos(td.size()));
			str += ", was expecting " + ttos(ntensors) + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		for (SizeType i = 0; i < ntensors; ++i) {
			PsimagLite::String name = td(i).name();
			SizeType id = td(i).id();
			PairStringSizeType p(name,id);
			VectorPairStringSizeType::iterator x = std::find(tensorNameIds_.begin(),
			                                                 tensorNameIds_.end(),
			                                                 p);
			if (x == tensorNameIds_.end()) {
				std::cerr<<"WARNING: Unused tensor name= "<<name<<" id= "<<id<<"\n";
				continue;
			}

			SizeType ind = x - tensorNameIds_.begin();
			assert(ind < tensors_.size());

			SizeType ins = td(i).ins();
			SizeType outs = td(i).outs();
			VectorSizeType dimensions(ins + outs);
			for (SizeType j = 0; j < ins; ++j) {
				SizeType legTag = td(i).legTag(j,TensorStanza::INDEX_DIR_IN);
				dimensions[j] = legTag;
			}

			for (SizeType j = 0; j < outs; ++j) {
				SizeType legTag = td(i).legTag(j,TensorStanza::INDEX_DIR_OUT);
				dimensions[j+ins] = legTag;
			}

			assert(ind < tensors_.size());
			tensors_[ind] = new TensorType(dimensions);
		}
	}

	void findTensors(const TensorSrep& t,
	                 PsimagLite::String nameToOptimize,
	                 SizeType idToOptimize)
	{
		SizeType ntensors = t.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			PsimagLite::String name = t(i).name();
			SizeType id = t(i).id();
			bool conjugate = t(i).isConjugate();
			bool b = (name == nameToOptimize && id == idToOptimize && conjugate);
			if (!b && conjugate) continue;
			PairStringSizeType p(name,id);
			tensorNameIds_.push_back(p);
		}
	}

	void appendToMatrix(const TensorSrep& t)
	{
		SizeType total = t.maxTag('f') + 1;
		VectorSizeType freeIndices(total,0);
		VectorDirType directions(total,TensorStanza::INDEX_DIR_IN);
		VectorSizeType dimensions(total,0);
		prepareFreeIndices(dimensions,directions,t);
		PairSizeType rc = getRowsAndCols(dimensions,directions);
		if (m_.n_row() == 0) {
			m_.resize(rc.first, rc.second);
		} else if (m_.n_row() != rc.first || m_.n_col() != rc.second) {
			PsimagLite::String str("Hamiltonian terms environ \n");
			throw PsimagLite::RuntimeError(str);
		}

		do {
			PairSizeType rc = getRowAndColFromFree(freeIndices,dimensions,directions);
			m_(rc.first,rc.second) += 0.0;
		} while (TensorEvalType::nextIndex(freeIndices,dimensions));
	}

	void prepareFreeIndices(VectorSizeType& dimensions,
	                        VectorDirType& directions,
	                        const TensorSrep& t) const
	{
		SizeType n = dimensions.size();
		assert(n == directions.size());
		SizeType ntensors = t.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			PsimagLite::String name = t(i).name();
			SizeType id = t(i).id();
			PairStringSizeType p(name,id);
			VectorPairStringSizeType::const_iterator x = std::find(tensorNameIds_.begin(),
			                                                       tensorNameIds_.end(),
			                                                       p);
			if (x == tensorNameIds_.end()) {
				std::cerr<<"WARNING: Unused tensor name= "<<name<<" id= "<<id<<"\n";
				continue;
			}

			SizeType ind = x - tensorNameIds_.begin();
			assert(ind < tensors_.size());

			SizeType ins = t(i).ins();
			SizeType outs = t(i).outs();
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

	TensorOptimizer(const TensorOptimizer&);

	TensorOptimizer& operator=(const TensorOptimizer&);

	VectorTensorSrepType tensorSrep_;
	VectorPairStringSizeType tensorNameIds_;
	VectorTensorType tensors_;
	MatrixType m_;

}; // class TensorOptimizer
} // namespace Mera
#endif // TENSOROPTIMIZER_H
