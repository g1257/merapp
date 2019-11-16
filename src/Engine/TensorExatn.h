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

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorComplexOrRealType;
	typedef std::pair<PsimagLite::String, SizeType> PairStringSizeType;
	typedef ComplexOrRealType* TensorBlobType;

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
		// FIXME: SET TENSOR to a constant value value

//		// use std:: function here instead of loop, FIXME
//		SizeType n = data_.size();
//		for (SizeType i = 0; i < n; ++i)
//			data_[i] = value;
		exatn::initTensor(name_, value);
	}

	class SetToMatrix : public exatn::TensorMethod {
	public:

		SetToMatrix(const MatrixType& m) : m_(m)
		{}

		virtual void pack(BytePacket & packet) {}
		virtual void unpack(BytePacket & packet) {}
		virtual const std::string name() const {return "SetToMatrix";}
		virtual const std::string description() const {return "SetToMatrixDescription"; }

	    //Application-defined external tensor method:
	    virtual int apply(talsh::Tensor& local_tensor)
		{
			const ComplexOrRealType* ptr = 0;
			bool ret = local_tensor.getDataAccessHostConst(&ptr);
			//if (!ret)
				// check error
			unsigned int numdims = 0;
			const int* iptr = local_tensor.getDimExtents(numdims);

			// ptr[i0 + i1*iptr[0] + i2*iptr[0]*iptr[1] ] = ...m_(i, j);
			return 0;
		}

	private:

		const MatrixType& m_;
	};

	void setToMatrix(const MatrixType& m)
	{
		if (ins_ == 0) return;
		if (dimensions_.size() < ins_) return;

		std::shared_ptr<SetToMatrix> setToMatrix1 = std::make_shared<SetToMatrix>(m);

		exatn::transformTensor(name_, setToMatrix1);
		// FIXME: Set tensor to matrix m as below:

//		if (dimensions_.size() == ins_) {
//			SizeType rows = m.n_row();
//			SizeType cols = m.n_col();
//			for (SizeType i = 0; i < rows; ++i)
//				for (SizeType j = 0; j < cols; ++j)
//					data_[i+j*rows] = m(i,j);
//			return;
//		}

//		SizeType dins = 1;
//		for (SizeType i = 0; i < ins_; ++i)
//			dins *= dimensions_[i];		return data_[0];


//		SizeType douts = 1;
//		for (SizeType i = ins_; i < dimensions_.size(); ++i)
//			douts *= dimensions_[i];
//		for (SizeType x = 0; x < dins; ++x)
//			for (SizeType y = 0; y < douts; ++y)
//				data_[x + y*dins] = m(x,y);
	}

//	void setSizes(const VectorSizeType& dimensions)
//	{
//		if (ins_ > dimensions.size())
//			throw PsimagLite::RuntimeError("Tensor::setSizes(...): dimensions < ins\n");

//		dimensions_ = dimensions;

//		SizeType v = volume();
//		if (v == 0)
//			throw PsimagLite::RuntimeError("Tensor::setSizes(...): dimensions == 0\n");

//		//data_ = std::make_shared<exatn::numerics::Tensor>(name_, exatn::numerics::TensorShape(dimensions_));
//	}

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

//	SizeType index(const VectorSizeType& args) const
//	{
//		return pack(args);
//	}

	const ComplexOrRealType& operator()(const VectorSizeType& args) const
	{
		// FIXME: Return tensor data_ at args, const version

		/* SizeType index = pack(args);
		assert(index < data_.size());
		return data_[index];*/


	}

	// FIXME: GIVES AWAY INTERNALS!!
//	ComplexOrRealType& operator()(const VectorSizeType& args)
//	{
//		// FIXME: Return tensor data_ at args, non const version
//		/*
//		SizeType index = pack(args);
//		assert(index < data_.size());
//		return data_[index];*/


//	}

//	SizeType pack(const VectorSizeType& args) const
//	{
//		assert(args.size() > 0);
//		assert(args.size() == dimensions_.size());
//		SizeType index = 0;
//		SizeType prod = 1;

//		for (SizeType i = 0; i < args.size(); ++i) {
//			if (dimensions_[i] == 0) continue;
//			assert(args[i] < dimensions_[i]);
//			index += args[i]*prod;
//			prod *= dimensions_[i];
//		}

//		return index;
//	}

	const ComplexOrRealType* data() const
	{

		std::shared_ptr<talsh::Tensor> ptr = exatn::getLocalTensor(name_);
		const ComplexOrRealType** ptr2;
		bool ret = ptr->getDataAccessHostConst(ptr2);
		checkTalshErrorCode(ret, "getLocalTensor");
		return *ptr2;
	}

	void setData(const ComplexOrRealType* data)
	{
		std::shared_ptr<talsh::Tensor> ptr = exatn::getLocalTensor(name_);
		ComplexOrRealType** ptr2;
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

	Tensor(const Tensor&) = delete;

	Tensor& operator=(const Tensor&) = delete;

	void checkTalshErrorCode(bool code, PsimagLite::String what) const
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
