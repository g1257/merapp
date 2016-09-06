#ifndef PARAMETERSFORSOLVER_H
#define PARAMETERSFORSOLVER_H

namespace Mera {

struct ParametersForSolver {

	ParametersForSolver()
	    : tauMax(3)
	{}

	SizeType tauMax;
}; // struct ParametersForSolver
} // namespace Mera
#endif // PARAMETERSFORSOLVER_H
