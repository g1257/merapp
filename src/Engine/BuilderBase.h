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
#ifndef MERA_BUILDERBASE_H
#define MERA_BUILDERBASE_H
#include "AllocatorCpu.h"
#include "TensorSrep.h"

namespace Mera {
class BuilderBase {

	typedef TensorSrep::VectorPairSizeType VectorPairSizeType;
	typedef TensorSrep::PairSizeType PairSizeType;

public:

	virtual ~BuilderBase() { }

	virtual const PsimagLite::String& srep() const = 0;

	virtual TensorSrep* buildEnergyTerm(SizeType site,
	                                    const TensorSrep& tensorSrep) const = 0;

protected:

	TensorSrep* energyTerm(SizeType c,
	                       const TensorSrep& tensorSrep,
	                       SizeType dimension,
	                       SizeType sites) const
	{
		SizeType connections = dimension*sites;
		div_t q = div(c, dimension);
		SizeType direction = q.rem;
		SizeType site = q.quot;
		SizeType sitep = findSitePrime(site, direction, dimension, sites);

		TensorSrep tensorSrep2(tensorSrep);
		tensorSrep2.conjugate();
		PsimagLite::String str3("h");
		str3 += ttos(c) + "(";
		str3 += "f" + ttos(connections) + ",";
		str3 += "f" + ttos(connections + 1) + "|";
		str3 += "f" + ttos(site) + ",";
		str3 += "f" + ttos(sitep) + ")";
		TensorSrep tensorSrep3(str3);
		TensorSrep::VectorSizeType indicesToContract(2,site);
		indicesToContract[1] = sitep;
		TensorSrep* tensorSrep4 = new TensorSrep(tensorSrep);
		bool relabel = false;
		tensorSrep4->contract(tensorSrep3,indicesToContract,relabel);
		bool verbose = false;
		if (!tensorSrep4->isValid(verbose))
			throw PsimagLite::RuntimeError("Invalid tensor\n");
		VectorPairSizeType replacements;
		replacements.push_back(PairSizeType(connections, site));
		replacements.push_back(PairSizeType(connections + 1, sitep));
		correctFreeIndicesBeforeContraction(*tensorSrep4, replacements);

		std::cerr<<"UPPER"<<site<<"_"<<sitep<<"="<<tensorSrep4->sRep()<<"\n";
		std::cerr<<"LOWER"<<site<<"_"<<sitep<<"="<<tensorSrep2.sRep()<<"\n";

		relabel = true;
		tensorSrep4->contract(tensorSrep2, relabel);
		std::cerr<<"e"<<site<<"_"<<sitep<<"="<<tensorSrep4->sRep()<<"\n";
		if (!tensorSrep4->isValid(verbose))
			throw PsimagLite::RuntimeError("Invalid tensor\n");
		return tensorSrep4;
	}

	SizeType snapBack(int x, int l) const
	{
		if (x >= 0 && x < l) return x;
		while (x < 0) x += l;
		while (x >= l) x -= l;
		return x;
	}

private:

	SizeType findSitePrime(SizeType site,
	                       SizeType direction,
	                       SizeType dimension,
	                       SizeType sites) const
	{
		if (dimension == 1) return findSitePrime1d(site, sites);
		if (dimension == 2) return findSitePrime2d(site, direction, sites);
		throw PsimagLite::RuntimeError("findSitePrime: dimension must be 1 or 2\n");
	}

	SizeType findSitePrime1d(SizeType site,
	                         SizeType sites) const
	{
		return (site + 1 == sites) ? 0 : site + 1;
	}

	SizeType findSitePrime2d(SizeType site,
	                       SizeType direction,
	                       SizeType sites) const
	{
		SizeType l = sqrt(sites);
		assert(l*l == sites);
		div_t q = div(site, l);
		int x = q.quot;
		int y = q.rem;
		if (direction == 0)
			return snapBack(x + 1, l)*l + y;
		return x*l + snapBack(y + 1, l);
	}


	void correctFreeIndicesBeforeContraction(TensorSrep& t,
	                                         VectorPairSizeType& replacements) const
	{
		for (SizeType i = 0; i < t.size(); ++i) {
			TensorStanza& stanza = t(i);
			stanza.replaceSummedOrFrees(replacements, 'f');
			stanza.refresh();
		}

		t.refresh();
		t.isValid(true);
	}
};
}

#endif // MERA_BUILDERBASE_H
