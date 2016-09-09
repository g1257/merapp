
/** \ingroup MERA */
/*@{*/

/*! \file main.cpp
 *
 *  The MERA main driver
 *
 */

#include <unistd.h>
#include "MeraSolver.h"
#include <fstream>
#include "MeraToTikz.h"
#include "Version.h"

void testCase(SizeType num,
              Mera::ParametersForSolver::MeraArchitectureEnum arch,
              SizeType sites,
              SizeType tauMax)
{
	std::cout<<"TEST NUMBER "<<num<<" ";
	std::cout<<Mera::ParametersForSolver::archToString(arch);
	std::cout<<" sites="<<sites<<" tauMax="<<tauMax<<"\n";
	Mera::ParametersForSolver params(arch,sites,tauMax);
	Mera::MeraSolver solver(params);
	solver.computeGroundState();
	std::cout<<solver;
	std::cout<<"-----------------------------------\n\n\n";

	PsimagLite::String file("meraTikzTestNumber");
	file += ttos(num);
	file += ".tex";
	std::ofstream fout(file.c_str());
	Mera::MeraToTikz<double> obj(solver.sRep(),params.tauMax);
	fout<<obj;
	fout.close();
}

int main(int ,char **argv)
{
	SizeType counter = 0;
	std::cout<<argv[0]<<" version "<<MERA_VERSION<<"\n";
	testCase(counter++,Mera::ParametersForSolver::MERA_ARCH_1D_BINARY,16,3);
	//testCase(counter++,Mera::ParametersForSolver::MERA_ARCH_1D_TERNARY,35,3);
}
