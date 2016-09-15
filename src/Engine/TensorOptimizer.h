#ifndef TENSOROPTIMIZER_H
#define TENSOROPTIMIZER_H
#include "Vector.h"
#include "IoSimple.h"
#include "TensorSrep.h"

namespace Mera {

template<typename ComplexOrRealType>
class TensorOptimizer {

	//typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::Vector<TensorSrep*>::Type VectorTensorSrepType;

public:

	typedef PsimagLite::IoSimple::In IoInType;

	TensorOptimizer(IoInType& io, SizeType number)
	{
		PsimagLite::String str("Y_FOR_TENSOR" + ttos(number));
		io.advance(str);
		SizeType terms = 0;
		io.readline(terms,"TERMS=");
		std::cerr<<"Read "<<terms<<" for tensor number "<<number<<"\n";
		tensorSrep_.resize(terms,0);
		for (SizeType i = 0; i < terms; ++i) {
			PsimagLite::String srep;
			io.readline(srep,"STRING=");
			std::cerr<<"Read string "<<srep<<"\n";
			tensorSrep_[i] = new TensorSrep(srep);
		}

		std::cerr<<"TensorOptimizer::ctor() done\n";
	}

	~TensorOptimizer()
	{
		SizeType terms = tensorSrep_.size();
		for (SizeType i = 0; i < terms; ++i) {
			delete tensorSrep_[i];
			tensorSrep_[i] = 0;
		}
	}

	void optimize()
	{}

private:

	TensorOptimizer(const TensorOptimizer&);

	TensorOptimizer& operator=(const TensorOptimizer&);

	VectorTensorSrepType tensorSrep_;

}; // class TensorOptimizer
} // namespace Mera
#endif // TENSOROPTIMIZER_H
