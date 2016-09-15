#ifndef TENSOROPTIMIZER_H
#define TENSOROPTIMIZER_H
#include "Vector.h"
#include "IoSimple.h"
#include "TensorSrep.h"
#include "TensorEval.h"
#include <algorithm>
#include "Sort.h"

namespace Mera {

bool operator==(const std::pair<PsimagLite::String,SizeType>& p1,
                const std::pair<PsimagLite::String,SizeType>& p2)
{
	if (p1.first != p2.first) return false;
	if (p1.second != p2.second) return false;
	return true;
}

bool operator!=(const std::pair<PsimagLite::String,SizeType>& p1,
                const std::pair<PsimagLite::String,SizeType>& p2)
{
	if (p1.first != p2.first) return true;
	if (p1.second != p2.second) return true;
	return false;
}

template<typename ComplexOrRealType>
class TensorOptimizer {

	//typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
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
		}

		PsimagLite::Sort<VectorPairStringSizeType> sort;
		VectorSizeType perm(tensorNameIds_.size(),0);
		sort.sort(tensorNameIds_,perm);
		SizeType end = (std::unique(tensorNameIds_.begin(),
		                            tensorNameIds_.end()) -
		                tensorNameIds_.begin());
		tensorNameIds_.resize(end);
		std::cerr<<"TensorOptimizer::ctor() done\n";
	}

	~TensorOptimizer()
	{
		SizeType terms = tensorSrep_.size();
		for (SizeType i = 0; i < terms; ++i) {
			delete tensorSrep_[i];
			tensorSrep_[i] = 0;
		}
	}

	void optimize()
	{}

private:

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
			if (name == "r") continue;
			PairStringSizeType p(name,id);
			tensorNameIds_.push_back(p);
		}
	}

	TensorOptimizer(const TensorOptimizer&);

	TensorOptimizer& operator=(const TensorOptimizer&);

	VectorTensorSrepType tensorSrep_;
	VectorTensorType tensors_;
	VectorPairStringSizeType tensorNameIds_;
}; // class TensorOptimizer
} // namespace Mera
#endif // TENSOROPTIMIZER_H
