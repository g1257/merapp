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
#ifndef MERA_HUBBARD_H
#define MERA_HUBBARD_H
#include "Matrix.h"
#include "Vector.h"

namespace  Mera {

template<typename ModelBaseType>
class Hubbard : public ModelBaseType {

public:

	typedef typename ModelBaseType::ComplexOrRealType ComplexOrRealType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::RealType RealType;
	typedef typename ModelBaseType::VectorRealType VectorRealType;
	typedef typename ModelBaseType::MatrixType MatrixType;
	typedef typename ModelBaseType::VectorMatrixType VectorMatrixType;
	typedef typename ModelBaseType::VectorSizeType VectorSizeType;

	Hubbard(const VectorType& v)
	    : twoSiteHam_(v.size(),0), shift_(0.0), qOne_(4, 0)
	{
		for (SizeType i = 0; i < qOne_.size(); ++i)
			qOne_[i] = i;
		SizeType h = 4;
		SizeType h2 = h*h;
		SizeType n = twoSiteHam_.size();

		VectorMatrixType cm(2, static_cast<MatrixType*>(0));
		buildCreationMatrices(cm);
		for (SizeType i = 0; i < n; ++i) {
			twoSiteHam_[i] = new MatrixType(h2,h2);
			shift_ += setTwoSiteHam(*(twoSiteHam_[i]),i,v[i],cm);
		}

		std::cout<<"Shift="<<shift_<<"\n";
	}

	~Hubbard()
	{
		for (SizeType i = 0; i < twoSiteHam_.size(); ++i) {
			delete twoSiteHam_[i];
			twoSiteHam_[i] = 0;
		}
	}

	const RealType& energyShift() const
	{
		return shift_;
	}

	const MatrixType& twoSiteHam(SizeType id) const
	{
		assert(id < twoSiteHam_.size());
		assert(twoSiteHam_[id]);
		return *(twoSiteHam_[id]);
	}

	const VectorSizeType& qOne() const { return qOne_; }

private:

	RealType setTwoSiteHam(MatrixType& m,
	                       SizeType,
	                       const ComplexOrRealType& hopping,
	                       const VectorMatrixType& cm)
	{
		SizeType h2 = m.n_row();
		assert(m.n_row() == m.n_col());
		MatrixType cmt;
		MatrixType tmp(h2, h2);
		MatrixType tmpt;
		for (SizeType spin = 0; spin < 2; ++spin) {
			transposeConjugate(cmt, *(cm[spin]));
			cmt *= hopping;
			externalProductDense(tmp, cmt, *(cm[spin]), spin);
			m += tmp;
			transposeConjugate(tmpt, tmp);
			m += tmpt;
		}

		return normalizeHam(m);
	}

	RealType normalizeHam(MatrixType& m2) const
	{
		MatrixType m = m2;
		SizeType n = m.n_row();
		VectorRealType eigs(n,0.0);
		diag(m,eigs,'N');
		assert(n - 1 < eigs.size());
		RealType diagCorrection = eigs[n-1];
		std::cerr<<"MeraSolver: DiagonalCorrection= "<<diagCorrection<<"\n";
		for (SizeType i = 0; i < n; ++i)
			m2(i,i) -= diagCorrection;
		return diagCorrection;
	}

	void buildCreationMatrices(VectorMatrixType& cm) const
	{
		SizeType h = qOne_.size();
		assert(cm.size() == 2); // UP and DOWN

		SizeType i = 0;
		cm[i] = new MatrixType(h, h);
		cm[i]->operator()(1,0) = 1;
		cm[i]->operator()(3,2) = 1;

		i = 1;
		cm[i] = new MatrixType(h, h);
		cm[i]->operator()(2,0) = 1;
		cm[i]->operator()(3,1) = -1;
	}

	void externalProductDense(MatrixType& dest,
	                          const MatrixType& a,
	                          const MatrixType& b,
	                          SizeType spin) const
	{
		SizeType h2 = dest.n_row();
		assert(h2 == dest.n_col());
		SizeType h = a.n_row();
		assert(a.n_col() == h);
		assert(h == b.n_row());
		assert(h == b.n_col());
		assert(h2 == h*h);
		for (SizeType i = 0; i < h2; ++i) {
			div_t qi = div(i, h);
			RealType sign = (i & 2) ? -1.0 : 1.0;
			for (SizeType j = 0; j < h2; ++j) {
				if (spin == 1) sign = (j & 1) ? -1.0 : 1.0;
				div_t qj = div(j, h);
				dest(i,j) = a(qi.rem, qj.rem)*b(qi.quot, qj.quot)*sign;
			}
		}
	}

	VectorMatrixType twoSiteHam_;
	RealType shift_;
	VectorSizeType qOne_;
};
}
#endif // MERA_HUBBARD_H
