#ifndef NAMETOINDEXLUT_H
#define NAMETOINDEXLUT_H
#include "Vector.h"
#include "Sort.h"

namespace Mera {

template<typename TensorType>
class NameToIndexLut {

public:

	typedef typename PsimagLite::Vector<TensorType*>::Type VectorTensorType;
	typedef std::pair<PsimagLite::String, SizeType> PairStringSizeType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

	NameToIndexLut(const VectorTensorType& v)
	{
		redo(v);
	}

	void push(PsimagLite::String name)
	{
		const SizeType last = map_.size();
		map_[name] = last;
	}

	void redo(const VectorTensorType& v)
	{
		map_.clear();
		const SizeType n = v.size();
		for (SizeType i = 0; i < n; ++i)
			map_[v[i]->name()] = i;
	}

	void redo2(const VectorTensorType& v)
	{
		map_.clear();
		const SizeType n = v.size();
		VectorStringType names(n);
		VectorSizeType indices(n);
		for (SizeType i = 0; i < n; ++i) {
			names[i] = v[i]->name();
			indices[i] = i;
		}

		PsimagLite::Sort<VectorStringType> sort;
		sort.sort(names, indices);

		VectorSizeType inverse(n);
		for (SizeType i = 0; i < n; ++i)
			inverse[indices[i]] = i;

		for (SizeType i = 0; i < n; ++i)
			map_[names[i]]  = inverse[i];
	}

	SizeType operator()(PsimagLite::String str)
	{
		return map_.at(str);
	}

	SizeType operator()(const PairStringSizeType& p)
	{
		return map_.at(p.first + ttos(p.second));
	}

private:

	std::map<PsimagLite::String, SizeType> map_;
};
}
#endif // NAMETOINDEXLUT_H
