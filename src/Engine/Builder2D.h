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

public:

	typedef TensorSrep::VectorPairSizeType VectorPairSizeType;
	typedef TensorSrep::PairSizeType PairSizeType;

	Builder2D(SizeType sites, SizeType arity, bool isPeriodic)
	    : srep_("")
	{
		if (arity != 4)
			throw PsimagLite::RuntimeError("MeraBuilder2D: arity must be 4 for now\n");

		throw PsimagLite::RuntimeError("2D Mera is work in progress\n");

		SizeType x = sqrt(sites);
		if (x*x != sites) sitesNotSupported();
		if (x & 1) sitesNotSupported();
		x /= 2;
		if (x & 1) sitesNotSupported();

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


	TensorSrep* buildEnergyTerm(SizeType site,
	                            SizeType sites,
	                            const TensorSrep& tensorSrep) const
	{
		TensorSrep tensorSrep2(tensorSrep);
		tensorSrep2.conjugate();
		PsimagLite::String str3("h");
		str3 += ttos(site) + "(f";
		str3 += ttos(site+2) + ",f";
		str3 += ttos(site+3) + "|f";
		str3 += ttos(site) + ",f";
		SizeType x = (site + 1 == sites) ? 0 : site + 1;
		str3 += ttos(x) + ")\n";
		TensorSrep tensorSrep3(str3);
		TensorSrep::VectorSizeType indicesToContract(2,site);
		indicesToContract[1] = x;
		TensorSrep* tensorSrep4 = new TensorSrep(tensorSrep);
		tensorSrep4->contract(tensorSrep3,indicesToContract);
		if (!tensorSrep4->isValid(true))
			throw PsimagLite::RuntimeError("Invalid tensor\n");
		correctFreeIndicesBeforeContraction(*tensorSrep4, site, sites);

		std::cerr<<"LOWER"<<site<<"="<<tensorSrep2.sRep()<<"\n";
		std::cerr<<"UPPER"<<site<<"="<<tensorSrep4->sRep()<<"\n";
		tensorSrep4->contract(tensorSrep2);
		//std::cout<<"e"<<site<<"()="<<tensorSrep4->sRep()<<"\n";
		if (!tensorSrep4->isValid(true))
			throw PsimagLite::RuntimeError("Invalid tensor\n");
		return tensorSrep4;
	}

private:

	void sitesNotSupported() const
	{
		throw PsimagLite::RuntimeError("Builder2D: Requires sites == (4*x)^2 for x integer\n");
	}

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

	void correctFreeIndicesBeforeContraction(TensorSrep& t,
	                                         SizeType site,
	                                         SizeType sites) const
	{
		if (site < 1) return;

		SizeType x = (site + 1 == sites) ? 0 : site + 1;
		t.swapFree(1,x);
		t.swapFree(0,site);
		t.refresh();

		SizeType max = site/2;
		if (site & 1) ++max;
		for (SizeType i = 0; i < t.size(); ++i) {
			TensorStanza& stanza = t(i);
			if (stanza.name() != "u") continue;
			SizeType id = stanza.id();
			if (id > max) continue;
			SizeType ins = stanza.ins();
			SizeType k = 0;
			VectorPairSizeType replacements;
			for (SizeType j = 0; j < ins; ++j) {
				if (stanza.legType(j) != TensorStanza::INDEX_TYPE_FREE)
					continue;
				SizeType legTag = stanza.legTag(j);
				SizeType shouldBe = 2*id + k;
				if (shouldBe > site) break;
				if (legTag != shouldBe)
					replacements.push_back(PairSizeType(legTag, shouldBe));
				++k;
			}

			stanza.replaceSummedOrFrees(replacements, 'f');
			stanza.refresh();
		}

		t.refresh();

		if (site + 1 == sites) {
			t.swapFree(0,site);
			t.swapFree(1,site);
			t.refresh();
		}

		t.isValid(true);
	}

	PsimagLite::String srep_;
};
}
#endif // MERA_BUILDER2D_H
