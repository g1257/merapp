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

	SizeType args() const { return dimensions_.size(); }

	SizeType argSize(SizeType ind) const
	{
		assert(ind < dimensions_.size());
		return dimensions_[ind];
	}

	const ComplexOrRealType& operator()(const VectorSizeType& args) const
	{
		SizeType index = pack(args);
		assert(index < data_.size());
		return data_[index];
	}

private:

	SizeType pack(const VectorSizeType& args) const
	{
		SizeType index = args[0];
		SizeType prod = dimensions_[0];

		for (SizeType i = 1; i < args.size(); ++i) {
			index += args[i]*prod;
			prod *= dimensions_[i];
		}

		return index;
	}

	VectorSizeType dimensions_;
	VectorComplexOrRealType data_;
};
}
#endif // TENSOR_H
