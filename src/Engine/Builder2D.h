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
#ifndef MERA_BUILDER2D_H
#define MERA_BUILDER2D_H
#include "BuilderBase.h"
#include "TensorSrep.h"

namespace Mera {

class Builder2D : public BuilderBase {

	static const SizeType DIMENSION = 2;

public:

	typedef TensorSrep::VectorSizeType VectorSizeType;
	typedef TensorSrep::VectorPairSizeType VectorPairSizeType;
	typedef TensorSrep::PairSizeType PairSizeType;

	Builder2D(SizeType sites, SizeType arity, bool isPeriodic)
	    : sites_(sites), srep_("")
	{
		if (arity != 4)
			throw PsimagLite::RuntimeError("MeraBuilder2D: arity must be 4 for now\n");

		SizeType x = sqrt(sites);
		if (x*x != sites) sitesNotSupported();
		if (x & 1) sitesNotSupported();
		x /= 2;
		if (x & 1) sitesNotSupported();

		SizeType tensors = sites/4;
		SizeType counter = 0;
		SizeType summed = 0;
		SizeType savedSummedForU = 0;
		SizeType idsU = 0;
		SizeType idsW = 0;
		while (tensors > 1) {
			SizeType savedSummedForW = summed;
			createUlayer(summed,
			             idsU,
			             savedSummedForU,
			             tensors,
			             counter);
			savedSummedForU = summed;
			createWlayer(summed,
			             idsW,
			             savedSummedForW,
			             tensors);
			tensors /= 4;
			counter++;
		}

		assert(summed > 4);
		summed -= 4;
		PsimagLite::String rArgs("");
		for (SizeType i = 0; i < 4; ++i) {
			rArgs += "s" + ttos(summed++);
			if (i < 3) rArgs += ",";
		}

		srep_ += "r0(" + rArgs + ")";
	}

	const PsimagLite::String& srep() const { return srep_; }


	TensorSrep* buildEnergyTerm(SizeType c,
	                            const TensorSrep& tensorSrep) const
	{
		return BuilderBase::energyTerm(c, tensorSrep, DIMENSION, sites_);
	}

private:

	void sitesNotSupported() const
	{
		PsimagLite::String str("Builder2D: Requires sites == (4*x)^2");
		str += " for some integer x\n";
		throw PsimagLite::RuntimeError(str);
	}

	void createUlayer(SizeType& summed,
	                  SizeType& idsU,
	                  SizeType savedSummed,
	                  SizeType n,
	                  SizeType counter)
	{
		VectorSizeType assignment(4);
		SizeType sqrtN = sqrt(n);
		assert(sqrtN*sqrtN == n);
		PsimagLite::String fOrS = (counter == 0) ? "f" : "s";
		for (SizeType i = 0; i < n; ++i) {
			getUAssignment(assignment, i, sqrtN);
			assert(assignment.size() == 4);

			PsimagLite::String inArgs("");
			for (SizeType i = 0; i < 4; ++i) {
				SizeType tmp = assignment[i];
				if (counter > 0) tmp += savedSummed;
				inArgs += fOrS + ttos(tmp);
				if (i < 3) inArgs += ",";
			}

			srep_ += "u" + ttos(idsU++) + "(" + inArgs + "|";

			PsimagLite::String outArgs("");
			for (SizeType i = 0; i < 4; ++i) {
				SizeType tmp = assignment[i];
				if (counter > 0) tmp += summed;
				outArgs += "s" + ttos(tmp);
				if (i < 3) outArgs += ",";
			}

			srep_ += outArgs + ")";
		}

		summed += 4*n;
	}

	void getUAssignment(VectorSizeType& assignment,
	                    SizeType ind,
	                    SizeType l) const
	{
		div_t q = div(ind, l);
		SizeType x = 2*q.rem;
		SizeType y = 2*q.quot;
		SizeType counter = 0;
		SizeType n = 2*l;
		assignment[counter++] = x + y*n;
		assignment[counter++] = x + 1 + y*n;
		assignment[counter++] = x + (y+1)*n;
		assignment[counter++] = x + 1 + (y+1)*n;
	}

	void createWlayer(SizeType& summed,
	                  SizeType& idsW,
	                  SizeType savedSummed,
	                  SizeType n)
	{
		SizeType sqrtN = sqrt(n);
		assert(sqrtN*sqrtN == n);
		VectorSizeType assignment(4, 0);
		for (SizeType i = 0; i < n; ++i) {
			getWAssignment(assignment, i, sqrtN);
			assert(assignment.size() == 4);

			PsimagLite::String inArgs("");
			for (SizeType i = 0; i < 4; ++i) {
				SizeType tmp = assignment[i] + savedSummed;
				inArgs += "s" + ttos(tmp);
				if (i < 3) inArgs += ",";
			}

			srep_ += "w" + ttos(idsW++) +"(" + inArgs + "|";

			srep_ += "s" + ttos(summed++) + ")";
		}
	}

	void getWAssignment(VectorSizeType& assignment,
	                    SizeType ind,
	                    SizeType l) const
	{
		div_t q = div(ind, l);
		int x = 2*q.rem;
		int y = 2*q.quot;
		SizeType counter = 0;
		int n = 2*l;
		SizeType xm1 = this->snapBack(x - 1, n);
		SizeType ym1 = this->snapBack(y - 1, n);

		assignment[counter++] = xm1 + ym1*n;
		assignment[counter++] = x + ym1*n;
		assignment[counter++] = xm1 + y*n;
		assignment[counter++] = x + y*n;
	}

	SizeType sites_;
	PsimagLite::String srep_;
};
}
#endif // MERA_BUILDER2D_H
