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
#ifndef SrepStatement_H
#define SrepStatement_H
#include "Vector.h"
#include "TensorSrep.h"
#include <map>

namespace  Mera {

template<typename ComplexOrRealType>
class SrepStatement {

public:

	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef PsimagLite::Vector<PairSizeType>::Type VectorPairSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef std::pair<PsimagLite::String,SizeType> PairStringSizeType;
	typedef typename PsimagLite::Vector<PairStringSizeType>::Type VectorPairStringSizeType;
	typedef std::map<PairStringSizeType,SizeType> MapPairStringSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef TensorSrep TensorSrepType;

	SrepStatement(PsimagLite::String str)
	    : lhs_(0),rhs_(0)
	{
		VectorStringType vstr;
		PsimagLite::split(vstr, str, "=");
		if (vstr.size() != 2)
			throw PsimagLite::RuntimeError("SrepStatement:: syntax error " + str + "\n");
		lhs_ = new TensorStanza(vstr[0]);
		rhs_ = new TensorSrepType(vstr[1]);

		nameIdOfOutput_ = PairStringSizeType(lhs_->name(), lhs_->id());
	}

	SrepStatement(const SrepStatement& other)
	    : lhs_(0),rhs_(0)
	{
		lhs_ = new TensorStanza(*(other.lhs_));
		rhs_ = new TensorSrepType(*(other.rhs_));
		nameIdOfOutput_ = PairStringSizeType(lhs_->name(), lhs_->id());
	}

	~SrepStatement()
	{
		delete lhs_;
		delete rhs_;
	}

	void canonicalize() const
	{
		VectorPairSizeType frees;
		computeFrees(frees);

		if (frees.size() == 0) return;
		assert(frees[0].second == 0);
		rhs_->simplify(frees);
		lhs_->replaceSummedOrFrees(frees,'f');
		lhs_->refresh();
	}

	PsimagLite::String sRep() const
	{
		return lhs_->sRep() + "=" + rhs_->sRep();
	}

	const TensorStanza& lhs() const
	{
		assert(lhs_);
		return *lhs_;
	}

	const TensorSrepType& rhs() const
	{
		assert(rhs_);
		return *rhs_;
	}

	const PairStringSizeType& nameIdOfOutput() const
	{
		return nameIdOfOutput_;
	}

	// FIXME: Gives away internals!
	TensorSrepType& rhs()
	{
		assert(rhs_);
		return *rhs_;
	}

private:

	void computeFrees(VectorPairSizeType& replacements) const
	{
		SizeType legs = lhs_->legs();
		SizeType counter = 0;
		for (SizeType j = 0; j < legs; ++j) {
			if (lhs_->legType(j) != TensorStanza::INDEX_TYPE_FREE)
				continue;
			replacements.push_back(PairSizeType(lhs_->legTag(j),counter++));
		}
	}

	SrepStatement& operator=(const SrepStatement&);

	TensorStanza* lhs_;
	TensorSrepType* rhs_;
	PairStringSizeType nameIdOfOutput_;
}; // class SrepStatement

} // namespace Mera

#endif // SrepStatement_H
