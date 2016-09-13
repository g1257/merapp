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
	typedef typename PsimagLite::Vector<TensorType*>::Type VectorTensorType;

	TensorEval(PsimagLite::String srep,const VectorTensorType& vt)
	    : tensorSrep_(srep), data_(vt)
	{}

	ComplexOrRealType eval() const
	{
		return 0.0;
	}

private:

	TensorSrepType tensorSrep_;
	const VectorTensorType& data_;
};
}
#endif // MERA_TENSOREVAL_H
