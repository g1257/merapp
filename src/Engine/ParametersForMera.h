#ifndef ParametersForMera_H
#define ParametersForMera_H
#include "Vector.h"
#include "IoSimple.h"

namespace Mera {

template<typename ComplexOrRealType_>
struct ParametersForMera {

	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	ParametersForMera(const VectorType& hTerms,
	                  const VectorSizeType& qOne1,
	                  SizeType m1,
	                  PsimagLite::String eval,
	                  RealType tol)
	    : hamiltonianConnection(hTerms),
	      qOne(qOne1),
	      m(m1),
	      verbose(false),
	      evaluator(eval),
	      tolerance(tol)
	{}

	ParametersForMera(PsimagLite::String filename)
	{
		PsimagLite::IoSimple::In io(filename);
		io.readline(options,"MeraOptions=");
		io.read(hamiltonianConnection, "hamiltonianConnection");
		io.read(qOne, "qOne");
		io.readline(m, "m=");
		int x = 0;
		io.readline(x, "verbose=");
		verbose = (x > 0);
		io.readline(evaluator, "evaluator=");
		io.readline(tolerance, "Tolerance=");
	}

	PsimagLite::String options;
	VectorType hamiltonianConnection;
	VectorSizeType qOne;
	SizeType m;
	bool verbose;
	PsimagLite::String evaluator;
	RealType tolerance;
}; // struct ParametersForMera

template<typename T>
std::ostream& operator<<(std::ostream& os, const ParametersForMera<T>& p)
{
	os<<"MeraOptions=none\n";
	os<<"hamiltonianConnection ";
	os<<p.hamiltonianConnection;
	os<<"qOne\n";
	os<<p.qOne<<"\n";
	os<<"m="<<p.m<<"\n";
	os<<"verbose="<<((p.verbose) ? 1 : 0)<<"\n";
	os<<"evaluator="<<p.evaluator<<"\n";
	os<<"Tolerance="<<p.tolerance<<"\n";
	return os;
}

} // namespace Mera
#endif // ParametersForMera_H
