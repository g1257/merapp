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
#include "MeraToTikz.h"
#include <iostream>
#include "Vector.h"
#include <cstdlib>
#include "Version.h"
#include "IoSimple.h"

int main(int argc, char** argv)
{
	if (argc < 2) {
		std::cerr<<"USAGE: "<<argv[0]<<" filename\n";
		return 1;
	}

	PsimagLite::String srep("");
	PsimagLite::IoSimple::In io(argv[1]);
	SizeType sites = 0;
	io.readline(sites,"Sites=");
	io.rewind();
	io.readline(srep,"Srep=");

	Mera::MeraToTikz<double> obj(srep,sites);
	std::cout<<"%Created by "<<argv[0]<<" version "<<MERA_VERSION<<"\n";
	std::cout<<obj;
}
