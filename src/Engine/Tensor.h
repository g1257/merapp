/*
Copyright (c) 2016, UT-Battelle, LLC

MERA++, Version 0.

This file is part of MERA++.
MERA++ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
MERA++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with MERA++. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef TENSOR_H
#define TENSOR_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "Matrix.h"
#include "RandomForTests.h"

namespace Mera {

template<typename ComplexOrRealType>
class Tensor {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorComplexOrRealType;

	Tensor(SizeType dim0, SizeType ins)
	    : dimensions_(1,dim0),data_(dim0),ins_(ins),rng_(1234)
	{}

	Tensor(const VectorSizeType& d, SizeType ins)
	    : dimensions_(d),ins_(ins),rng_(1234)
	{
		SizeType n = dimensions_.size();
		if (n == 0) return;
		assert(0 < n);
		SizeType prod = dimensions_[0];
		for (SizeType i = 1; i < dimensions_.size(); ++i)
			prod *= dimensions_[i];
		data_.resize(prod,0);
	}

	void setToIdentity(ComplexOrRealType value)
	{
		setToIdentityEx(value,ins_);
	}

	void setToIdentityEx(ComplexOrRealType value, SizeType ins)
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
			if (x < douts)
				data_[x + x*dins] = value;
	}

	void setRandom()
	{
		SizeType n = data_.size();
		ComplexOrRealType sum = 0.0;
		for (SizeType i = 0; i < n; ++i) {
			ComplexOrRealType value = rng_();
			data_[i] = value;
			sum += value*PsimagLite::conj(value);
		}

		RealType tmp = PsimagLite::real(sum);
		assert(tmp > 1e-6);
		tmp = 2.0/sqrt(tmp);
		for (SizeType i = 0; i < n; ++i)
			data_[i] *= tmp;
	}

	void setToMatrix(const MatrixType& m)
	{
		if (ins_ == 0) return;
		if (dimensions_.size() < ins_) return;
		if (dimensions_.size() == ins_) {
			data_ = m();
			return;
		}

		SizeType dins = 1;
		for (SizeType i = 0; i < ins_; ++i)
			dins *= dimensions_[i];

		SizeType douts = 1;
		for (SizeType i = ins_; i < dimensions_.size(); ++i)
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

	SizeType ins() const { return ins_; }

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
	SizeType ins_;
	PsimagLite::RandomForTests<RealType> rng_;
};
}
#endif // TENSOR_H
