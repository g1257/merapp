
/** \ingroup MERA */
/*@{*/

/*! \file main.cpp
 *
 *  The MERA main driver
 *
 */

#include <unistd.h>
#include "MeraSolver.h"

int main(int argc,char *argv[])
{
	Mera::ParametersForSolver params;
	Mera::MeraSolver solver(params);
	solver.computeGroundState();
	std::cout<<solver;
}
