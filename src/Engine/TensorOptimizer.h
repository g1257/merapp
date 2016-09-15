#ifndef TENSOROPTIMIZER_H
#define TENSOROPTIMIZER_H
#include "Vector.h"
#include "IoSimple.h"
#include "TensorSrep.h"
#include "TensorEval.h"
#include <algorithm>
#include "Sort.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorOptimizer {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<TensorSrep*>::Type VectorTensorSrepType;
	typedef std::pair<PsimagLite::String,SizeType> PairStringSizeType;
	typedef PsimagLite::Vector<PairStringSizeType>::Type VectorPairStringSizeType;

public:

	typedef PsimagLite::IoSimple::In IoInType;
	typedef Mera::TensorEval<double> TensorEvalType;
	typedef TensorEvalType::TensorType TensorType;
	typedef TensorEvalType::VectorTensorType VectorTensorType;

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

	TensorOptimizer(const TensorOptimizer&);

	TensorOptimizer& operator=(const TensorOptimizer&);

	VectorTensorSrepType tensorSrep_;
	VectorPairStringSizeType tensorNameIds_;
	VectorTensorType tensors_;

}; // class TensorOptimizer
} // namespace Mera
#endif // TENSOROPTIMIZER_H
