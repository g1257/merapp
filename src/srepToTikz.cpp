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

int main(int argc, char** argv)
{
	if (argc == 1) {
		std::cerr<<"USAGE: "<<argv[0]<<" tauMax\n";
		return 1;
	}

	char c = '\0';
	PsimagLite::String srep;
	while (std::cin>>c) {
		srep += c;
	}

	Mera::MeraToTikz<double> obj(srep,atoi(argv[1]));
	std::cout<<"%Created by "<<argv[0]<<" version "<<MERA_VERSION<<"\n";
	std::cout<<obj;
}
