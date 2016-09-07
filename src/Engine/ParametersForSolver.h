#ifndef PARAMETERSFORSOLVER_H
#define PARAMETERSFORSOLVER_H

namespace Mera {

struct ParametersForSolver {

	enum MeraArchitectureEnum {MERA_ARCH_1D_BINARY,MERA_ARCH_1D_TERNARY,
		                       MERA_ARCH_2D_2X2,MERA_ARCH_2D_3X3};

	ParametersForSolver()
	    : meraArch(MERA_ARCH_1D_BINARY),tauMax(3),numberOfSites(16)
	{}

	MeraArchitectureEnum meraArch;
	SizeType tauMax;
	SizeType numberOfSites;
}; // struct ParametersForSolver
} // namespace Mera
#endif // PARAMETERSFORSOLVER_H
