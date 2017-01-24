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
#ifndef TENSORLEG_H
#define TENSORLEG_H
#include "Vector.h"

namespace Mera {

class TensorLeg {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;
	typedef std::pair<char, SizeType> PairCharSizeType;

	enum IndexDirectionEnum {INDEX_DIR_IN, INDEX_DIR_OUT};

	TensorLeg(PsimagLite::String tag, IndexDirectionEnum inOrOut)
	    : inOrOut_(inOrOut)
	{
		tag_ = getPairCharInt(tag);
	}

	const char name() const { return tag_.first; }

	char& name() { return tag_.first; }

	const SizeType& numericTag() const { return tag_.second; }

	SizeType& numericTag() { return tag_.second; }

	IndexDirectionEnum dir() const { return inOrOut_; }

private:

	PairCharSizeType getPairCharInt(PsimagLite::String token) const
	{
		SizeType l = token.length();

		if (l == 0) {
			PsimagLite::String str("TensorStanza: malformed stanza, ");
			throw PsimagLite::RuntimeError(str + token + "\n");
		}

		std::size_t index = token.find_first_of("0123456789");

		if (index == 0) {
			// all numeric
			return PairCharSizeType('d',atoi(token.c_str()));
		}

		if (l == 1 && token != "d") {
			PsimagLite::String str("TensorStanza: malformed stanza, expecting a dummy, got ");
			throw PsimagLite::RuntimeError(str + token + "\n");
		}

		PsimagLite::String tmp = (l == 1) ? "0" : token.substr(1,l-1);
		return PairCharSizeType(token[0],atoi(tmp.c_str()));
	}

	PairCharSizeType tag_;
	IndexDirectionEnum inOrOut_;
}; // class TensorLeg
} // namespace Mera
#endif // TENSORLEG_H
