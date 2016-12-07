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
#ifndef PROGRAMGLOBALS_H
#define PROGRAMGLOBALS_H
#include "Vector.h"
#include "RandomForTests.h"

namespace Mera {

struct ProgramGlobals {

	typedef PsimagLite::RandomForTests<double> RngType;

	static SizeType logBase2Strict(SizeType x)
	{
		if (x == 0) return 0;
		SizeType counter = 0;
		SizeType ones = 0;
		while (x > 0) {
			if (x & 1) ones++;
			if (ones > 1) return 0;
			x >>= 1;
			counter++;
		}

		assert(counter > 0);
		return counter - 1;
	}


	static PsimagLite::String addLf(PsimagLite::String str, SizeType each)
	{
		SizeType l = str.length();
		PsimagLite::String str2("");
		for (SizeType i = 0; i < l; ++i) {
			if (i > 0 && i % each == 0) str2 += "\n";
			str2 += str[i];
		}

		return str2;
	}

	static RngType rng;
};

ProgramGlobals::RngType ProgramGlobals::rng(1234);
} // namespace Mera
#endif // PROGRAMGLOBALS_H
