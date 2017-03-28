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
#ifndef HEISENBERG_H
#define HEISENBERG_H
#include "Matrix.h"
#include "Vector.h"

namespace  Mera {

template<typename ModelBaseType>
class Heisenberg : public ModelBaseType {

public:

	typedef typename ModelBaseType::ComplexOrRealType ComplexOrRealType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::RealType RealType;
	typedef typename ModelBaseType::VectorRealType VectorRealType;
	typedef typename ModelBaseType::MatrixType MatrixType;
	typedef typename ModelBaseType::VectorMatrixType VectorMatrixType;

	Heisenberg(const VectorType& v)
	    : twoSiteHam_(v.size(),0), shift_(0.0)
	{
		SizeType h = 2; // model dependency here
		SizeType h2 = h*h;
		SizeType n = twoSiteHam_.size();

		for (SizeType i = 0; i < n; ++i) {
			twoSiteHam_[i] = new MatrixType(h2,h2);
			shift_ += setTwoSiteHam(*(twoSiteHam_[i]),i,v);
		}

		std::cout<<"Shift="<<shift_<<"\n";
	}

	~Heisenberg()
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

private:

	RealType setTwoSiteHam(MatrixType& m, SizeType site, const VectorType& v)
	{
		SizeType n = m.n_row();
		assert(n == m.n_col());

		assert(site < v.size());
		ComplexOrRealType scale = v[site];

		// Sz Sz
		for (SizeType i = 0; i < n; ++i)
			m(i,i) = (i == 0 || i ==3) ? 0.25*scale : -0.25*scale;

		// S+S- S-S+
		m(1,2) = m(2,1) = 0.5*scale;

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
		std::cout<<"MeraSolver: DiagonalCorrection= "<<diagCorrection<<"\n";
		for (SizeType i = 0; i < n; ++i)
			m2(i,i) -= diagCorrection;
		return diagCorrection;
	}

	VectorMatrixType twoSiteHam_;
	RealType shift_;
};
}
#endif // HEISENBERG_H
