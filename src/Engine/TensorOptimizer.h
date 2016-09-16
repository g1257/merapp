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
	typedef PsimagLite::Vector<TensorStanza::IndexDirectionEnum>::Type VectorDirType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;

public:

	typedef PsimagLite::IoSimple::In IoInType;
	typedef Mera::TensorEval<ComplexOrRealType> TensorEvalType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename TensorEvalType::PairStringSizeType PairStringSizeType;
	typedef typename TensorEvalType::VectorPairStringSizeType VectorPairStringSizeType;
	typedef typename TensorEvalType::TensorType TensorType;
	typedef typename TensorEvalType::VectorTensorType VectorTensorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef std::pair<SizeType,SizeType> PairSizeType;

	TensorOptimizer(IoInType& io,
	                PsimagLite::String nameToOptimize,
	                SizeType idToOptimize)
	    : tensorToOptimize_(nameToOptimize,idToOptimize),twoSiteHam_(4,4)
	{
		bool makeHamTheIdentity = false;
		setTwoSiteHam(makeHamTheIdentity);

		initTensorSreps(io,nameToOptimize,idToOptimize);

		initTensorNameIds();

		PsimagLite::String dstr("");
		io.rewind();
		io.readline(dstr,"DIMENSION_SREP=");
		initTensors(dstr);

		std::cerr<<"TensorOptimizer::ctor() done\n";
	}

	~TensorOptimizer()
	{
		SizeType terms = tensorSrep_.size();
		for (SizeType i = 0; i < terms; ++i) {
			delete tensorSrep_[i];
			delete energySrep_[i];
			tensorSrep_[i] = energySrep_[i] = 0;
		}

		terms = tensors_.size();
		for (SizeType i = 0; i < terms; ++i) {
			delete tensors_[i];
			tensors_[i] = 0;
		}
	}

	void optimize(SizeType iters)
	{
		for (SizeType iter = 0; iter < iters; ++iter) {
			optimizeInternal(iter);
			std::cout<<"energy="<<calcEnergy()<<"\n";
		}
	}

private:

	void optimizeInternal(SizeType)
	{
		SizeType terms = tensorSrep_.size();
		MatrixType m;
		for (SizeType i = 0; i < terms; ++i) {
			appendToMatrix(m,*(tensorSrep_[i]));
			std::cerr<<m;
			std::cerr<<"------------------\n";
		}

		std::cerr<<"About to do svd...\n";
		VectorRealType s(m.n_row(),0);
		MatrixType vt(m.n_col(),m.n_col());
		svd('A',m,s,vt);
		page14StepL3(m,vt);
	}

	// FIXME: pick up model dependency here
	void setTwoSiteHam(bool testWithIdentity)
	{
		for (SizeType i = 0; i < 4; ++i)
			twoSiteHam_(i,i) = 1.0;
		if (testWithIdentity) return;

		for (SizeType i = 0; i < 4; ++i) {
			// Sz Sz
			twoSiteHam_(i,i) = (i == 0 || i ==3) ? 0.25 : -0.25;
			if (i == 3) continue;
			for (SizeType j = 0; j < 3; ++j) {
				if (i == j) continue;
				twoSiteHam_(i,j) = 0.5; // S+S- + S-S+
			}
		}
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
					sum += vt(k,i)*m(j,k);
				t(i,j) = -std::conj(sum);
			}
		}

		typename VectorPairStringSizeType::const_iterator it = std::find(tensorNameIds_.begin(),
		                                                                 tensorNameIds_.end(),
		                                                                 tensorToOptimize_);
		if (it == tensorNameIds_.end()) {
			PsimagLite::String str("Not found tensor to optimize in tensors\n");
			throw PsimagLite::RuntimeError(str);
		}

		SizeType ind = it - tensorNameIds_.begin();
		SizeType ins = 2; // FIXME
		tensors_[ind]->setToMatrix(ins,t);
	}

	ComplexOrRealType calcEnergy() const
	{
		VectorSizeType freeIndices;
		ComplexOrRealType sum = 0.0;
		for (SizeType i = 0; i < energySrep_.size(); ++i) {
			TensorEvalType eval(energySrep_[i]->sRep(),tensors_,tensorNameIds_);
			sum += eval(freeIndices);
		}

		return sum;
	}

	void initTensorSreps(IoInType& io,
	                     PsimagLite::String nameToOptimize,
	                     SizeType idToOptimize)
	{
		SizeType terms = 0;
		io.readline(terms,"TERMS=");
		std::cerr<<"Read "<<terms<<" for tensor id "<<idToOptimize<<"\n";
		tensorSrep_.resize(terms,0);
		energySrep_.resize(terms,0);

		for (SizeType i = 0; i < terms; ++i) {
			PsimagLite::String srep;
			io.readline(srep,"ENERGY=");
			energySrep_[i] = new TensorSrep(srep);
			std::cerr<<"Free indices= "<<(1+energySrep_[i]->maxTag('f'))<<"\n";

			io.readline(srep,"STRING=");
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
			typename VectorPairStringSizeType::iterator x = std::find(tensorNameIds_.begin(),
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
			if (name == "h") {
				tensors_[ind]->setToMatrix(ins,twoSiteHam_);
			} else if (name == "r") {
				assert(0 < dimensions.size());
				tensors_[ind]->setToIdentity(1,1.0/sqrt(2.0));
			} else {
				tensors_[ind]->setToIdentity(ins,1.0);
			}
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

	void appendToMatrix(MatrixType& m, const TensorSrep& t)
	{
		//std::cerr<<"SREP= "<<t.sRep()<<"\n";
		SizeType total = t.maxTag('f') + 1;
		VectorSizeType freeIndices(total,0);
		VectorDirType directions(total,TensorStanza::INDEX_DIR_IN);
		VectorSizeType dimensions(total,0);
		VectorBoolType conjugate(total,false);
		prepareFreeIndices(dimensions,directions,conjugate,t);
		modifyDirections(directions,conjugate);
		PairSizeType rc = getRowsAndCols(dimensions,directions);
		bool invert = false;
		if (m.n_row() == 0) {
			m.resize(rc.first, rc.second);
		} else if (m.n_row() == rc.second && m.n_col() == rc.first) {
			invert = true;
		} else if (m.n_row() != rc.first || m.n_col() != rc.second) {
			PsimagLite::String str("Hamiltonian terms environ \n");
			throw PsimagLite::RuntimeError(str);
		}

		SizeType count = 0;
		do {
			PairSizeType rc = getRowAndColFromFree(freeIndices,dimensions,directions);
			TensorEvalType eval(t.sRep(),tensors_,tensorNameIds_);
			ComplexOrRealType tmp = eval(freeIndices);
			if (invert)
				m(rc.second,rc.first) += tmp;
			else
				m(rc.first,rc.second) += tmp;
			count++;
		} while (TensorEvalType::nextIndex(freeIndices,dimensions));
		std::cerr<<count<<"\n";
	}

	void prepareFreeIndices(VectorSizeType& dimensions,
	                        VectorDirType& directions,
	                        VectorBoolType& conjugate,
	                        const TensorSrep& t) const
	{
		SizeType n = dimensions.size();
		assert(n == directions.size());
		assert(n = conjugate.size());

		SizeType ntensors = t.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			PsimagLite::String name = t(i).name();
			SizeType id = t(i).id();
			PairStringSizeType p(name,id);
			typename VectorPairStringSizeType::const_iterator x = std::find(tensorNameIds_.begin(),
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

		for (SizeType i = 0; i < n; ++i)
			if (conjugate[i]) dirs[i] = oppositeDir(dirs[i]);
	}

	TensorStanza::IndexDirectionEnum oppositeDir(TensorStanza::IndexDirectionEnum dir) const
	{
		TensorStanza::IndexDirectionEnum in = TensorStanza::INDEX_DIR_IN;
		TensorStanza::IndexDirectionEnum out = TensorStanza::INDEX_DIR_OUT;
		return (dir == in) ? out : in;
	}

	TensorOptimizer(const TensorOptimizer&);

	TensorOptimizer& operator=(const TensorOptimizer&);

	PairStringSizeType tensorToOptimize_;
	VectorTensorSrepType tensorSrep_;
	VectorTensorSrepType energySrep_;
	VectorPairStringSizeType tensorNameIds_;
	VectorTensorType tensors_;
	MatrixType twoSiteHam_;
}; // class TensorOptimizer
} // namespace Mera
#endif // TENSOROPTIMIZER_H
