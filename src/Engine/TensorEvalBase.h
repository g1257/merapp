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
#ifndef TENSOREVALBASE_H
#define TENSOREVALBASE_H
#include "Tensor.h"
#include "SrepEquation.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorEvalBase {

	class TensorEvalHandle {

	public:

		enum Status {STATUS_IDLE, STATUS_IN_PROGRESS, STATUS_DONE};

		TensorEvalHandle(Status status = STATUS_IDLE)
		    : status_(status)
		{}

		bool done() const
		{
			return (status_ == STATUS_DONE);
		}

	private:

		Status status_;
	};

public:

	typedef TensorEvalHandle HandleType;
	typedef Tensor<ComplexOrRealType> TensorType;
	typedef typename PsimagLite::Vector<TensorType*>::Type VectorTensorType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef SrepEquation<ComplexOrRealType> SrepEquationType;
	typedef typename SrepEquationType::PairStringSizeType PairStringSizeType;
	typedef std::map<PairStringSizeType,SizeType> MapPairStringSizeType;

	virtual HandleType operator()(bool cached) = 0;

};
} // namespace Mera
#endif // TENSOREVALBASE_H
