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
#ifndef MERABUILDER_H
#define MERABUILDER_H
#include "ProgramGlobals.h"
#include "TensorSrep.h"
#include "Vector.h"
#include "Builder1D.h"
#include "Builder2D.h"

namespace Mera {

template<typename ComplexOrRealType>
class MeraBuilder {

	typedef PsimagLite::Vector<TensorSrep*>::Type VectorTensorSrepType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

public:

	MeraBuilder(SizeType sites,
	            SizeType arity,
	            SizeType dimension,
	            bool isPeriodic,
	            const VectorType& hamTerms)
	    : sites_(sites),
	      arity_(arity),
	      dimension_(dimension),
	      isPeriodic_(isPeriodic),
	      srep_(""),
	      energy_(sites*dimension,0)
	{
		BuilderBase* builder = 0;

		if (dimension == 1) {
			builder = new Builder1D(sites, arity, isPeriodic);
		} else if (dimension == 2) {
			builder = new Builder2D(sites, arity, isPeriodic);
		} else {
			PsimagLite::String msg("MeraBuilder: dimension ");
			throw PsimagLite::RuntimeError(msg + ttos(dimension) + " is not supported\n");
		}

		srep_ = builder->srep();

		buildEnergies(*builder, hamTerms);

		delete builder;
	}

	~MeraBuilder()
	{
		for (SizeType site = 0; site < energy_.size(); ++site) {
			if (energy_[site] == 0) continue;
			delete energy_[site];
			energy_[site] = 0;
		}
	}

	const PsimagLite::String& operator()() const
	{
		return srep_;
	}

	const TensorSrep& energy(SizeType size) const
	{
		assert(size < energy_.size());
		assert(energy_[size]);
		return *(energy_[size]);
	}

	const SizeType& arity() const { return arity_; }

	const SizeType& sites() const { return sites_; }

	bool isPeriodic() const { return isPeriodic_; }

private:

	void buildEnergies(const BuilderBase& builder,
	                   const VectorType& hamTerm)
	{
		TensorSrep tensorSrep(srep_);
		SizeType connections = hamTerm.size();
		assert(connections == sites_*dimension_);
		assert(connections == energy_.size());
		for (SizeType c = 0; c < connections; ++c) {
			assert(hamTerm.size() > c);
			if (hamTerm[c] == 0.0) continue;
			assert(c < energy_.size());
			energy_[c] = builder.buildEnergyTerm(c, tensorSrep);
		}
	}

	SizeType sites_;
	SizeType arity_;
	SizeType dimension_;
	bool isPeriodic_;
	PsimagLite::String srep_;
	VectorTensorSrepType energy_;
}; // class MeraBuilder
} // namespace Mera
#endif // MERABUILDER_H
