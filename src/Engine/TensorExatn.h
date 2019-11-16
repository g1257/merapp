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
#ifndef TENSOR_EXATN_H
#define TENSOR_EXATN_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "Matrix.h"
#include "RandomForTests.h"
#include "exatn.hpp"
#include "../../exatn/tpls/ExaTensor/include/talshxx.hpp"

namespace Mera {

template<typename ComplexOrRealType>
class Tensor {

	Tensor(const Tensor&) = delete;

	Tensor& operator=(const Tensor&) = delete;

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorComplexOrRealType;
	typedef std::pair<PsimagLite::String, SizeType> PairStringSizeType;

	class BlobAndSizeConst {

	public:

		BlobAndSizeConst() : size_(0), ptr_(0) {}

		BlobAndSizeConst(SizeType size, const ComplexOrRealType* ptr)
		    : size_(size), ptr_(ptr)
		{}

		const ComplexOrRealType& operator[](SizeType ind) const
		{
			assert(ptr_);
			assert(ind < size_);
			return ptr_[ind];
		}

		SizeType size_;
		const ComplexOrRealType* ptr_;
	};

	typedef BlobAndSizeConst TensorBlobType;

	// Tensor with only one dimension
	Tensor(PsimagLite::String name, SizeType dim0, SizeType ins)
	    : name_(name),
	      dimensions_(1, dim0),
	      ins_(ins)
	{
		exatn::createTensor(name_, exatn::TensorElementType::REAL64, exatn::numerics::TensorShape(dimensions_));
	}

	Tensor(PsimagLite::String name, const VectorSizeType& d, SizeType ins)
	    : name_(name),
	      dimensions_(d),
	      ins_(ins)
	{
		exatn::createTensor(name_, exatn::TensorElementType::REAL64, exatn::numerics::TensorShape(dimensions_));
	}

	~Tensor()
	{
		exatn::destroyTensor(name_);
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

		// FIXME: set data_ as below:

//		for (SizeType x = 0; x < dins; ++x)
//			if (x < douts)
//				data_[x + x*dins] = value;
	}

	void setToRandom()
	{
		//  FIXME: set tensor to random

//		SizeType n = data_.size();
//		ComplexOrRealType sum = 0.0;
//		for (SizeType i = 0; i < n; ++i) {
//			ComplexOrRealType value = rng_();
//			data_[i] = value;
//			sum += value*PsimagLite::conj(value);
//		}

//		RealType tmp = PsimagLite::real(sum);
//		assert(tmp > 1e-6);
//		tmp = 2.0/sqrt(tmp);
//		for (SizeType i = 0; i < n; ++i)
//			data_[i] *= tmp;
	}

	void setToConstant(ComplexOrRealType value)
	{
		// SET TENSOR to a constant value value
		exatn::initTensor(name_, value);
	}

	class SetToMatrix : public exatn::TensorMethod {
	public:

		SetToMatrix(const VectorSizeType& dimensions,
		            SizeType ins,
		            const MatrixType& m) : dimensions_(dimensions), ins_(ins), m_(m)
		{}

		virtual void pack(BytePacket & packet) {}

		virtual void unpack(BytePacket & packet) {}

		virtual const std::string name() const {return "SetToMatrix";}

		virtual const std::string description() const {return "SetToMatrixDescription"; }

	    //Application-defined external tensor method:
	    virtual int apply(talsh::Tensor& local_tensor)
		{
			ComplexOrRealType* ptr = 0;
			bool ret = local_tensor.getDataAccessHost(&ptr);
			checkTalshErrorCode(ret, "getDataAccessHostConst");

			if (dimensions_.size() == ins_) {
				SizeType rows = m_.n_row();
				SizeType cols = m_.n_col();
				for (SizeType i = 0; i < rows; ++i)
					for (SizeType j = 0; j < cols; ++j)
						ptr[i+j*rows] = m_(i,j);
				return 0;
			}

			SizeType dins = 1;
			for (SizeType i = 0; i < ins_; ++i)
				dins *= dimensions_[i];

			SizeType douts = 1;
			for (SizeType i = ins_; i < dimensions_.size(); ++i)
				douts *= dimensions_[i];
			for (SizeType x = 0; x < dins; ++x)
				for (SizeType y = 0; y < douts; ++y)
					ptr[x + y*dins] = m_(x,y);
			// ptr[i0 + i1*iptr[0] + i2*iptr[0]*iptr[1] ] = ...m_(i, j);
			return 0; // error code
		}

	private:

		const VectorSizeType& dimensions_;
		const SizeType ins_;
		const MatrixType& m_;
	};

	void setToMatrix(const MatrixType& m)
	{
		if (ins_ == 0) return;

		if (dimensions_.size() < ins_) return;

		std::shared_ptr<SetToMatrix> setToMatrix1 = std::make_shared<SetToMatrix>(dimensions_,
		                                                                          ins_,
		                                                                          m);

		exatn::transformTensor(name_, setToMatrix1);
	}

	void setSizes(const VectorSizeType& dimensions)
	{
		if (ins_ > dimensions.size())
			throw PsimagLite::RuntimeError("Tensor::setSizes(...): dimensions < ins\n");

		dimensions_ = dimensions;

		SizeType v = volume();
		if (v == 0)
			throw PsimagLite::RuntimeError("Tensor::setSizes(...): dimensions == 0\n");

		exatn::destroyTensor(name_);
		exatn::createTensor(name_, exatn::TensorElementType::REAL64, exatn::numerics::TensorShape(dimensions_));
	}

	SizeType volume() const
	{
		if (dimensions_.size() == 0) return 0;
		SizeType prod = dimensions_[0];
		for (SizeType i = 1; i < dimensions_.size(); ++i)
			prod *= dimensions_[i];
		return prod;
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
		// Return tensor data_ at args, const version

		SizeType index = pack(args);
		TensorBlobType tbt = this->data();
		return tbt[index];
	}

	void setValue(const VectorSizeType& args, const ComplexOrRealType& val)
	{
		SizeType index = pack(args);

		std::shared_ptr<talsh::Tensor> localTensor = exatn::getLocalTensor(name_);
		ComplexOrRealType* ptr = 0;
		bool ret = localTensor->getDataAccessHost(&ptr);
		checkTalshErrorCode(ret, "getDataAccessHostConst");

		ptr[index] = val;
	}

	TensorBlobType data() const
	{
		std::shared_ptr<talsh::Tensor> ptr = exatn::getLocalTensor(name_);
		const ComplexOrRealType** ptr2 = 0;
		bool ret = ptr->getDataAccessHostConst(ptr2);
		checkTalshErrorCode(ret, "getLocalTensor");
		return TensorBlobType(ptr->getVolume(), *ptr2);
	}

	void setData(const TensorBlobType& data)
	{
		std::shared_ptr<talsh::Tensor> ptr = exatn::getLocalTensor(name_);
		ComplexOrRealType** ptr2 = 0;
		bool ret = ptr->getDataAccessHost(ptr2);
		checkTalshErrorCode(ret, "getLocalTensor");

		const SizeType n = ptr->getVolume();
		for (SizeType i = 0; i < n; ++i)
			(*ptr2)[i] = data[i];

		// needs to be done as an add
		// data_ = data;
	}

	PsimagLite::String name() const { return name_; }

private:

	SizeType pack(const VectorSizeType& args) const
	{
		assert(args.size() > 0);
		assert(args.size() == dimensions_.size());
		SizeType index = 0;
		SizeType prod = 1;

		for (SizeType i = 0; i < args.size(); ++i) {
			if (dimensions_[i] == 0) continue;
			assert(args[i] < dimensions_[i]);
			index += args[i]*prod;
			prod *= dimensions_[i];
		}

		return index;
	}

	static void checkTalshErrorCode(bool code, PsimagLite::String what)
	{
		if (code) return;
		throw PsimagLite::RuntimeError("MERA++: TALSH returned false from " + what + "\n");
	}

	static PsimagLite::RandomForTests<ComplexOrRealType> rng_;
	PsimagLite::String name_;
	VectorSizeType dimensions_;
	SizeType ins_;
};

template<typename ComplexOrRealType>
PsimagLite::RandomForTests<ComplexOrRealType> Tensor<ComplexOrRealType>::rng_(1234);

}
#endif // TENSOR_EXATN_H
