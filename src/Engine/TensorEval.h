#ifndef MERA_TENSOREVAL_H
#define MERA_TENSOREVAL_H

#include "TensorSrep.h"
#include "Tensor.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorEval {

	typedef TensorSrep TensorSrepType;

public:

	typedef Tensor<ComplexOrRealType> TensorType;
	typedef typename TensorType::VectorSizeType VectorSizeType;
	typedef typename PsimagLite::Vector<TensorType*>::Type VectorTensorType;
	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef PsimagLite::Vector<PairSizeType> VectorPairSizeType;

	TensorEval(PsimagLite::String srep,const VectorTensorType& vt)
	    : tensorSrep_(srep), data_(vt)
	{}

	ComplexOrRealType eval() const
	{
		SizeType total = tensorSrep_.maxTag('s');
		VectorSizeType summed(total,0);
		VectorPairSizeType whereSummed;
		VectorSizeType free;
		ComplexOrRealType sum = 0.0;
		do {
			sum += evalInternal(summed,whereSummed,free);
		} while (nextSummed(summed));

		return sum;
	}

private:

	bool nextSummed(VectorSizeType& summed) const
	{
		return false;
	}

	ComplexOrRealType evalInternal(const VectorSizeType& summed,
	                               const VectorPairSizeType& whereSummed,
	                               const VectorSizeType& free) const
	{
		return 0.0;
	}


	TensorSrepType tensorSrep_;
	const VectorTensorType& data_;
};
}
#endif // MERA_TENSOREVAL_H
