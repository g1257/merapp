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
#ifndef MERA_BUILDER1D_H
#define MERA_BUILDER1D_H
#include "BuilderBase.h"
#include "TensorSrep.h"

namespace Mera {

class Builder1D : public BuilderBase {

	static const SizeType DIMENSION = 1;

public:

	typedef TensorSrep::VectorPairSizeType VectorPairSizeType;
	typedef TensorSrep::PairSizeType PairSizeType;

	Builder1D(SizeType sites, SizeType arity, bool isPeriodic)
	    : sites_(sites), srep_("")
	{
		if (arity != 2)
			throw PsimagLite::RuntimeError("MeraBuilder1D: arity must be 2 for now\n");

		SizeType ln = 0;
		if ((ln = ProgramGlobals::logBase2Strict(sites)) == 0)
			throw PsimagLite::RuntimeError("MeraBuilder: sites must be a power of 2\n");

		SizeType tensors = sites/2;
		SizeType counter = 0;
		SizeType summed = 0;
		SizeType savedSummedForU = 0;
		SizeType idsU = 0;
		SizeType idsW = 0;
		while (tensors > 1) {
			SizeType savedSummedForW = summed;
			SizeType pLastIndex = createUlayer(summed,
			                                   idsU,
			                                   savedSummedForU,
			                                   tensors,
			                                   isPeriodic,
			                                   counter);
			savedSummedForU = summed;
			createWlayer(summed,
			             idsW,
			             savedSummedForW,
			             tensors,
			             isPeriodic,
			             pLastIndex,
			             counter);
			tensors /= 2;
			counter++;
		}

		assert(summed > 1);
		srep_ += "r0(s" + ttos(summed-2) + ",s" + ttos(summed-1) + ")";
	}

	const PsimagLite::String& srep() const { return srep_; }


	TensorSrep* buildEnergyTerm(SizeType c,
	                            const TensorSrep& tensorSrep) const
	{
		return BuilderBase::energyTerm(c, tensorSrep, DIMENSION, sites_);
	}

private:

	SizeType createUlayer(SizeType& summed,
	                      SizeType& idsU,
	                      SizeType savedSummed,
	                      SizeType n,
	                      bool isPeriodic,
	                      SizeType counter)
	{
		PsimagLite::String I0("");
		PsimagLite::String I1("");
		PsimagLite::String O0("");
		PsimagLite::String O1("");
		SizeType periodicLastIndex = 0;
		for (SizeType i = 0; i < n; ++i) {
			bool hasTwoOutputs = true;
			const bool oddCounter = (counter & 1);

			if (counter == 0) {
				I0 = "f" + ttos(2*i);
				I1 = "f" + ttos(2*i+1);
			} else {
				I0 = "s" + ttos(savedSummed++);
				I1 = "s" + ttos(savedSummed++);
			}

			if (oddCounter) {
				if (i + 1 == n) hasTwoOutputs = false;
			} else {
				if (i == 0) hasTwoOutputs = false;
			}

			srep_ += "u" + ttos(idsU++) + "("+I0+"," + I1 + "|";

			if (!oddCounter && !hasTwoOutputs && isPeriodic) {
				periodicLastIndex = summed++;
				O1 = "s" + ttos(periodicLastIndex);
				srep_ +=  O1 + ",";
			}

			O0 = "s" + ttos(summed++);
			srep_ += O0;

			if (hasTwoOutputs) {
				O1 = "s" + ttos(summed++);
				srep_ += "," + O1;
			} else if (isPeriodic && oddCounter) {
				// border here
				periodicLastIndex = summed++;
				O1 = "s" + ttos(periodicLastIndex);
				srep_ += "," + O1;
			}

			srep_ += ")";
		}

		return periodicLastIndex;
	}

	void createWlayer(SizeType& summed,
	                  SizeType& idsW,
	                  SizeType savedSummed,
	                  SizeType n,
	                  bool isPeriodic,
	                  SizeType periodicLastIndex,
	                  SizeType counter)
	{
		PsimagLite::String I0("");
		PsimagLite::String I1("");
		PsimagLite::String O0("");
		const bool oddCounter = (counter & 1);
		if (isPeriodic && !oddCounter) ++savedSummed;

		for (SizeType i = 0; i < n; ++i) {
			bool hasTwoInputs = true;

			if (oddCounter) {
				if (i == 0) hasTwoInputs = false;
			} else {
				if (i + 1 == n) hasTwoInputs = false;
			}

			srep_ += "w" + ttos(idsW++) +"(";

			if (oddCounter && !hasTwoInputs && isPeriodic) {
				I1 = "s" + ttos(periodicLastIndex);
				srep_ +=  I1 + ",";
			}

			I0 = "s" + ttos(savedSummed++);
			srep_ += I0;

			if (hasTwoInputs) {
				I1 = "s" + ttos(savedSummed++);
				srep_ += "," + I1;
			} else if (isPeriodic && !oddCounter) {
				I1 = "s" + ttos(periodicLastIndex);
				srep_ += "," + I1;
			}

			O0 = "s" + ttos(summed++);
			srep_ += "|" + O0 + ")";
		}
	}

	SizeType sites_;
	PsimagLite::String srep_;
};
}
#endif // MERA_BUILDER1D_H
