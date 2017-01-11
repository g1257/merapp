#ifndef PARAMETERSFORSOLVER_H
#define PARAMETERSFORSOLVER_H
#include "Vector.h"
#include "IoSimple.h"

namespace Mera {

template<typename ComplexOrRealType_>
struct ParametersForSolver {

	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	ParametersForSolver(VectorType& hTerms,
	                    SizeType h1,
	                    SizeType m1)
	    : hamiltonianConnection(hTerms),h(h1),m(m1),verbose(false)
	{}

	ParametersForSolver(PsimagLite::String filename)
	{
		PsimagLite::IoSimple::In io(filename);
		io.read(hamiltonianConnection, "hamiltonianConnection");
		io.readline(h, "h=");
		io.readline(m, "m=");
		int x = 0;
		io.readline(x, "verbose=");
		verbose = (x > 0);
	}

	VectorType hamiltonianConnection;
	SizeType h;
	SizeType m;
	bool verbose;
}; // struct ParametersForSolver

template<typename T>
std::ostream& operator<<(std::ostream& os, const ParametersForSolver<T>& p)
{
	os<<"hamiltonianConnection ";
	os<<p.hamiltonianConnection;
	os<<"h="<<p.h<<"\n";
	os<<"m="<<p.m<<"\n";
	os<<"verbose="<<((p.verbose) ? 1 : 0)<<"\n";

	return os;
}

} // namespace Mera
#endif // PARAMETERSFORSOLVER_H
