#ifndef PARAMETERSFORSOLVER_H
#define PARAMETERSFORSOLVER_H

namespace Mera {

struct ParametersForSolver {

	ParametersForSolver(SizeType tauMax_)
	    : tauMax(tauMax_)
	{}

	SizeType tauMax;
}; // struct ParametersForSolver
} // namespace Mera
#endif // PARAMETERSFORSOLVER_H
