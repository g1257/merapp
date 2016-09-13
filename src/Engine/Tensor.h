#ifndef TENSOR_H
#define TENSOR_H
#include "Vector.h"
#include "ProgramGlobals.h"

namespace Mera {

template<typename ComplexOrRealType>
class Tensor {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorComplexOrRealType;

	Tensor(SizeType dim0)
	    : dimensions_(1,dim0),data_(dim0)
	{}

	// FIXME: Take into account Hermitian and unitary properties
	void setToRandom()
	{
		for (SizeType i = 0; i < data_.size(); ++i)
			data_[i] = 10.0*ProgramGlobals::rng();
	}

private:

	VectorSizeType dimensions_;
	VectorComplexOrRealType data_;
};
}
#endif // TENSOR_H
