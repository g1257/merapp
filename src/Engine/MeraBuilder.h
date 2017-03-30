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
	    : sites_(sites), arity_(arity), srep_(""), energy_(sites,0)
	{
		if (arity != 2)
			throw PsimagLite::RuntimeError("MeraBuilder: arity must be 2 for now\n");

		BuilderBase* builder = 0;

		if (dimension == 1) {
			builder = new Builder1D(sites, isPeriodic);
		} else if (dimension == 2) {
			builder = new Builder2D(sites, isPeriodic);
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

private:

	void buildEnergies(const BuilderBase& builder,
	                   const VectorType& hamTerm)
	{
		TensorSrep tensorSrep(srep_);
		SizeType sites = hamTerm.size();
		for (SizeType site = 0; site < sites; ++site) {
			if (hamTerm[site] == 0.0) continue;
			energy_[site] = builder.buildEnergyTerm(site, sites, tensorSrep);
		}
	}

	SizeType sites_;
	SizeType arity_;
	PsimagLite::String srep_;
	VectorTensorSrepType energy_;
}; // class MeraBuilder
} // namespace Mera
#endif // MERABUILDER_H
