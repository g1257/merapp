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
#include "ModelSelector.h"
#include "ModelBase.h"

void usageMain(const PsimagLite::String& str)
{
	throw PsimagLite::RuntimeError(str);
}

template<typename ComplexOrRealType>
void main1(const Mera::MeraBuilder<ComplexOrRealType>& builder,
           const Mera::ParametersForMera<ComplexOrRealType>& params)
{
	typedef Mera::ModelBase<ComplexOrRealType> ModelBaseType;
	Mera::ModelSelector<ModelBaseType> model(params.model, params.hamiltonianConnection);
	PsimagLite::String srep = builder();
	PsimagLite::String meraString = srep;
	PsimagLite::String hString = "D" + ttos(model().qOne().size());
	PsimagLite::String args = "(" + hString + "," + hString + "|" + hString + "," + hString + ")";
	for (SizeType i = 0; i < params.hamiltonianConnection.size(); ++i) {
		if (params.hamiltonianConnection[i] == 0.0) continue;
		srep += "h" + ttos(i) + args;
	}

	Mera::TensorSrep tsrep(srep);
	SizeType maxLegs = 2.0*params.hamiltonianConnection.size();
	Mera::SymmetryLocal symmLocal(tsrep.size(), model().qOne(), maxLegs);
	Mera::DimensionSrep<Mera::SymmetryLocal> dimSrep(srep, symmLocal, params.m);
	PsimagLite::String dsrep = dimSrep();

	Mera::MeraEnviron<ComplexOrRealType, Mera::SymmetryLocal> environ(builder,
	                                                                  params,
	                                                                  dsrep,
	                                                                  symmLocal);
	std::cout<<params;
	std::cout<<"IsMeraPeriodic="<<builder.isPeriodic()<<"\n";
	std::cout<<"NoSymmetryLocal=1\n";
	std::cout<<"IterMera=10\n";
	std::cout<<"IterTensor=100\n";
	std::cout<<"MERA="<<meraString<<"\n";

	Mera::DimensionSrep<Mera::SymmetryLocal>::printOnePerLine(std::cout,
	                                                          dsrep,
	                                                          "CreateTensor=");
	Mera::DimensionSrep<Mera::SymmetryLocal>::printOnePerLine(std::cout,
	                                                          environ.dimensionSrep(),
	                                                          "CreateTensor=");

	// add output u1000 to be used by unitary condition checking
	std::cout<<"DsrepEnvirons=u1000(D1,D1)"<<environ.dimensionSrep()<<"\n";
	std::cout<<environ.environs();
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

	while ((opt = getopt(argc, argv,"n:a:d:m:M:e:t:PbV")) != -1) {
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
		case 'm':
			m = atoi(optarg);
			break;
		case 'M':
			model = optarg;
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

	// sanity checks here
	if (sites*arity*dimension == 0 || sites == 1)
		usageMain(strUsage);

	assert(sites*dimension > 0);
	hamTerms.resize(sites*dimension,1.0);
	if (dimension == 1) {
		assert(sites > 0);
		assert(hamTerms.size() > sites - 1);
		hamTerms[sites-1] = (periodic) ? 1.0 : 0.0;
	}

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

	MeraParametersType params(hamTerms,m,evaluator,model,tolerance);
	main1(meraBuilder,params);
}
