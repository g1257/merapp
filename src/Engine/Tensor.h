#ifndef TENSOR_H
#define TENSOR_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "Matrix.h"

namespace Mera {

template<typename ComplexOrRealType>
class Tensor {

public:

	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorComplexOrRealType;

	Tensor(SizeType dim0)
	    : dimensions_(1,dim0),data_(dim0)
	{}

	Tensor(const VectorSizeType& d) : dimensions_(d)
	{
		SizeType n = dimensions_.size();
		if (n == 0) return;
		assert(0 < n);
		SizeType prod = dimensions_[0];
		for (SizeType i = 1; i < dimensions_.size(); ++i)
			prod *= dimensions_[i];
		data_.resize(prod,0);
	}

	// FIXME: Take into account Hermitian and unitary properties
	void setToRandom()
	{
		for (SizeType i = 0; i < data_.size(); ++i)
			data_[i] = 10.0*ProgramGlobals::rng();
	}

	void setTo(ComplexOrRealType value)
	{
		for (SizeType i = 0; i < data_.size(); ++i)
			data_[i] = value;
	}

	void setToIdentity(SizeType ins)
	{
		if (ins == 0) return;
		if (dimensions_.size() <= ins) return;

		SizeType dins = 1;
		for (SizeType i = 0; i < ins; ++i)
			dins *= dimensions_[i];

		SizeType douts = 1;
		for (SizeType i = ins; i < dimensions_.size(); ++i)
			douts *= dimensions_[i];
		for (SizeType x = 0; x < dins; ++x)
			if (x < douts) data_[x + x*dins] = 1.0;
	}

	void setToMatrix(SizeType ins, const MatrixType& m)
	{
		if (ins == 0) return;
		if (dimensions_.size() <= ins) return;

		SizeType dins = 1;
		for (SizeType i = 0; i < ins; ++i)
			dins *= dimensions_[i];

		SizeType douts = 1;
		for (SizeType i = ins; i < dimensions_.size(); ++i)
			douts *= dimensions_[i];
		for (SizeType x = 0; x < dins; ++x)
			for (SizeType y = 0; y < douts; ++y)
				data_[x + y*dins] = m(x,y);
	}

	const SizeType& dimension(SizeType j) const
	{
		assert(j < dimensions_.size());
		return dimensions_[j];
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
