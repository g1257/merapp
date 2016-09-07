#ifndef PARAMETERSFORSOLVER_H
#define PARAMETERSFORSOLVER_H

namespace Mera {

struct ParametersForSolver {

	enum MeraArchitectureEnum {MERA_ARCH_1D_BINARY,MERA_ARCH_1D_TERNARY,
		                       MERA_ARCH_2D_2X2,MERA_ARCH_2D_3X3};

	ParametersForSolver(MeraArchitectureEnum arch,
	                    SizeType numberOfSites_,
	                    SizeType tauMax_)
	    : meraArch(arch),numberOfSites(numberOfSites_),tauMax(tauMax_)
	{}

	static PsimagLite::String archToString(MeraArchitectureEnum arch)
	{
		switch (arch) {
		case MERA_ARCH_1D_BINARY:
			return "MERA_ARCH_1D_BINARY";
		case  MERA_ARCH_1D_TERNARY:
			return "MERA_ARCH_1D_TERNARY";
		case MERA_ARCH_2D_2X2:
			return "MERA_ARCH_2D_2X2";
		case MERA_ARCH_2D_3X3:
			return "MERA_ARCH_2D_3X3";
		default:
			return "MERA_ARCH_UNKNOWN";
		}
	}

	MeraArchitectureEnum meraArch;
	SizeType numberOfSites;
	SizeType tauMax;
}; // struct ParametersForSolver
} // namespace Mera
#endif // PARAMETERSFORSOLVER_H
