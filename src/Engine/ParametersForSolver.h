#ifndef PARAMETERSFORSOLVER_H
#define PARAMETERSFORSOLVER_H
#include "Vector.h"

namespace Mera {

struct ParametersForSolver {

	ParametersForSolver(SizeType tauMax_)
	    : tauMax(tauMax_)
	{}

	SizeType tauMax;
}; // struct ParametersForSolver
} // namespace Mera
#endif // PARAMETERSFORSOLVER_H
