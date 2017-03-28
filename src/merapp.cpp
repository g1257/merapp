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
#include "SymmetryLocal.h"

void usageMain(const PsimagLite::String& str)
{
	throw PsimagLite::RuntimeError(str);
}

template<typename ComplexOrRealType>
void main1(const Mera::MeraBuilder<ComplexOrRealType>& builder,
           const Mera::ParametersForMera<ComplexOrRealType>& params)
{
	PsimagLite::String srep = builder();
	PsimagLite::String meraString = srep;
	PsimagLite::String hString = "D" + ttos(params.qOne.size());
	PsimagLite::String args = "(" + hString + "," + hString + "|" + hString + "," + hString + ")";
	for (SizeType i = 0; i < params.hamiltonianConnection.size(); ++i) {
		if (params.hamiltonianConnection[i] == 0.0) continue;
		srep += "h" + ttos(i) + args;
	}

	Mera::TensorSrep tsrep(srep);
	SizeType maxLegs = 2.0*params.hamiltonianConnection.size();
	Mera::SymmetryLocal symmLocal(tsrep.size(), params.qOne, maxLegs);
	Mera::DimensionSrep<Mera::SymmetryLocal> dimSrep(srep, symmLocal, params.m);
	PsimagLite::String dsrep = dimSrep();

	Mera::MeraEnviron<ComplexOrRealType, Mera::SymmetryLocal> environ(builder,
	                                                                  params,
	                                                                  dsrep,
	                                                                  symmLocal);
	std::cout<<params;
	symmLocal.save(std::cout);
	std::cout<<"DimensionSrep="<<dsrep<<environ.dimensionSrep();
	// add output u1000 to be used by unitary condition checking
	std::cout<<"u1000(1,1)\n";
	SizeType iterMera = 5;
	std::cout<<"NoSymmetryLocal=1\n";
	std::cout<<"IterMera="<<iterMera<<"\n";
	std::cout<<"MERA="<<meraString<<"\n";

	std::cout<<environ.environs();
}

template<typename VectorType>
void fillHamTerms(VectorType& v,
                  PsimagLite::String terms)
{
	PsimagLite::Vector<PsimagLite::String>::Type tokens;
	PsimagLite::tokenizer(terms,tokens,",");
	if (tokens.size() & 1) {
		std::cerr<<"-s site0,value0,site1,value1,... expected\n";
		throw PsimagLite::RuntimeError("fillHamTerms\n");
	}

	for (SizeType i = 0; i < tokens.size(); i += 2) {
		SizeType ind = atoi(tokens[i].c_str());
		assert(ind < v.size());
		assert(i + 1 < tokens.size());
		v[ind] = atof(tokens[i+1].c_str());;
	}
}

int main(int argc, char **argv)
{
	// check for complex or real  here FIXME
	typedef Mera::ParametersForMera<double> MeraParametersType;
	typedef Mera::MeraBuilder<double> MeraBuilderType;

	int opt = 0;
	bool versionOnly = false;
	bool buildOnly = false;
	SizeType sites = 0;
	SizeType arity = 2;
	SizeType dimension = 1;
	double tolerance = 1e-4;
	bool periodic = false;
	MeraParametersType::VectorType hamTerms;
	SizeType m = 0;
	PsimagLite::String evaluator("slow");
	PsimagLite::String strUsage(argv[0]);
	PsimagLite::String model("Heisenberg");
	strUsage += " -n sites -a arity -d dimension [-M model] [-m m] ";
	strUsage += "| -S srep | -V\n";

	while ((opt = getopt(argc, argv,"n:a:d:m:M:s:e:t:PbV")) != -1) {
		switch (opt) {
		case 'n':
			sites = atoi(optarg);
			assert(sites > 1);
			hamTerms.resize(sites,1.0);
			hamTerms[sites-1] = 0.0;
			break;
		case 'a':
			arity = atoi(optarg);
			break;
		case 'd':
			dimension = atoi(optarg);
			break;
		case 'm':
			m = atoi(optarg);
			break;
		case 'M':
			model = optarg;
			break;
		case 's':
			if (hamTerms.size() == 0) {
				std::cerr<<argv[0]<<": option -s must be after -n\n";
				return 1;
			}

			std::fill(hamTerms.begin(),hamTerms.end(),0.0);
			fillHamTerms(hamTerms,optarg);
			break;
		case 'e':
			evaluator = optarg;
			break;
		case 't':
			tolerance = atof(optarg);
			break;
		case 'P':
			periodic = true;
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

	MeraParametersType::VectorSizeType qOne(2,0);
	qOne[1] = 1;

	if (model != "Heisenberg") qOne.clear();

	// sanity checks here
	if (qOne.size() == 0 || sites*arity*dimension == 0)
		usageMain(strUsage);

	if (periodic && sites - 1 < hamTerms.size())
		hamTerms[sites - 1] = 1;

	// here build srep
	MeraBuilderType meraBuilder(sites,arity,dimension,periodic,hamTerms);

	std::cout<<"Sites="<<sites<<"\n";

	if (buildOnly) {
		std::cout<<"Srep="<<meraBuilder()<<"\n";
		std::cerr<<argv[0]<<": Stoping here because of -b (build only). ";
		std::cerr<<"Not computing environments\n";
		return 1;
	}

	std::cout<<"#"<<argv[0]<<" version "<<MERA_VERSION<<"\n";

	MeraParametersType params(hamTerms,qOne,m,model,evaluator,tolerance);
	main1(meraBuilder,params);
}
