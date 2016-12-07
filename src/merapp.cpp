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
/** \ingroup MERA */
/*@{*/

/*! \file main.cpp
 *
 *  The MERA main driver
 *
 */

#include <unistd.h>
#include "MeraEnviron.h"
#include <fstream>
#include "MeraToTikz.h"
#include "Version.h"
#include "DimensionSrep.h"
#include "MeraBuilder.h"

void usageMain(const PsimagLite::String& str)
{
	throw PsimagLite::RuntimeError(str);
}

void main1(PsimagLite::String srep,
           SizeType h,
           SizeType m)
{
	Mera::DimensionSrep dimSrep(srep,h,m);
	PsimagLite::String dsrep = dimSrep();
	PsimagLite::String hString = ttos(h);
	dsrep += "h0(" + hString + "," + hString + "|" + hString + "," + hString + ")";
	dsrep += "i0(" + hString + "|" + hString + ")";
	std::cout<<"DimensionSrep="<<dsrep<<"\n";

	Mera::ParametersForSolver params;
	Mera::MeraEnviron environ(srep,params);

	environ.computeEnvirons();
	std::cerr<<environ;
}

int main(int argc, char **argv)
{
	int opt = 0;
	bool versionOnly = false;
	bool buildOnly = false;
	SizeType sites = 0;
	SizeType arity = 2;
	SizeType dimension = 1;
	SizeType h = 0;
	SizeType m = 0;
	PsimagLite::String srep("");
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -n sites -a arity -d dimension -h hilbertSize [-m m] ";
	strUsage += "| -S srep | -V\n";
	strUsage += "-h hilbertSize is always mandatory\n";

	while ((opt = getopt(argc, argv,"n:a:d:h:m:S:bV")) != -1) {
		switch (opt) {
		case 'n':
			sites = atoi(optarg);
			break;
		case 'a':
			arity = atoi(optarg);
			break;
		case 'd':
			dimension = atoi(optarg);
			break;
		case 'h':
			h = atoi(optarg);
			break;
		case 'm':
			m = atoi(optarg);
			break;
		case 'S':
			srep = optarg;
			break;
		case 'b':
			buildOnly = true;
			break;
		case 'V':
			versionOnly = true;
			break;
		default:
			usageMain(strUsage);
			return 1;
		}
	}

	std::cerr<<"#"<<argv[0]<<" version "<<MERA_VERSION<<"\n";

	if (versionOnly)
		return 0;

	// sanity checks here
	if (h == 0)
		usageMain(strUsage);

	bool buildMode = (sites*arity*dimension > 0);
	if (buildMode && srep != "") {
		strUsage += "Either n*a*d > 0 or srep != "" but not both\n";
		usageMain(strUsage);
	}

	if (!buildMode && srep == "") {
		srep = "u0(f0,f1|s0)u1(f2,f3|s1,s2)w0(s0,s1|s3)w1(s2|s4)r(s3,s4)\n";
		sites = 4;
	}

	if (buildMode) {
		// here build srep
		Mera::MeraBuilder meraBuilder(sites,arity,dimension);
		srep = meraBuilder();
	}

	if (buildOnly) {
		std::cout<<"Sites="<<sites<<"\n";
		std::cout<<"Srep="<<srep<<"\n";
		std::cerr<<argv[0]<<": Stoping here because of -b (build only). ";
		std::cerr<<"Not computing environments\n";
		return 1;
	}

	std::cout<<"#"<<argv[0]<<" version "<<MERA_VERSION<<"\n";
	main1(srep,h,m);
}
