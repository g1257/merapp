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
#include "TensorBreakup.h"
#include "TensorEval.h"

int main(int argc, char**argv)
{
	PsimagLite::String str = "u0(s0,s1|f0)h0(s0,s1|s2,s3)u0*(s2,s3|f1)";
	if (argc > 1) str = argv[1];

	Mera::TensorSrep srep(str);
	Mera::TensorBreakup tensorBreakup(srep);
	Mera::TensorBreakup::VectorStringType vstr;
	tensorBreakup(vstr);
}
