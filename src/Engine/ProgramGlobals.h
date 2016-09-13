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

namespace Mera {

struct ProgramGlobals {

	static SizeType packTimeSpace(SizeType time, SizeType space, SizeType tauMax)
	{
		return time + space*tauMax;
	}

	static void unpackTimeAndSpace(SizeType& time,
	                               SizeType& space,
	                               SizeType id,
	                               SizeType tauMax)
	{
		time = id % tauMax;
		space = id/tauMax;
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
};
} // namespace Mera
#endif // PROGRAMGLOBALS_H
