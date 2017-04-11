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
#define USE_PTHREADS_OR_NOT_NG
#include "Vector.h"
#include "MeraSolver.h"
#include "Version.h"
int main(int argc, char** argv)
{
	PsimagLite::String file = "";
	int opt = 0;
	int precision = 6;
	SizeType threads = 1;
	while ((opt = getopt(argc, argv,"f:p:t:")) != -1) {
		switch (opt) {
		case 'f':
			file = optarg;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 't':
			threads = atoi(optarg);
			break;
		}
	}

	if (file == "") {
		std::cerr<<"USAGE: "<<argv[0]<<" -f filename\n";
		return 1;
	}

	PsimagLite::Concurrency c(&argc, &argv, threads);
	std::cout<<"#MERA_VERSION="<<MERA_VERSION<<"\n";
	Mera::MeraSolver<double> meraSolver(file);

	meraSolver.optimize();
}
