/*
Copyright (c) 2016-2017, UT-Battelle, LLC

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
#ifndef MERA_MODEL_BASE
#define MERA_MODEL_BASE
#include "Matrix.h"
#include "Vector.h"

namespace  Mera {

template<typename ComplexOrRealType_>
class ModelBase {

public:

	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<MatrixType*>::Type VectorMatrixType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	virtual ~ModelBase() {}

	virtual const RealType& energyShift() const = 0;

	virtual const MatrixType& twoSiteHam(SizeType id) const = 0;

	virtual const VectorSizeType& qOne() const = 0;
};
}
#endif // HEISENBERG_H
