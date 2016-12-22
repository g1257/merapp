#ifndef SREPEQUATION_H
#define SREPEQUATION_H
#include "Tensor.h"
#include "Vector.h"
#include "TensorSrep.h"
#include <map>

namespace  Mera {

template<typename ComplexOrRealType>
class SrepEquation {

public:

	typedef Tensor<ComplexOrRealType> TensorType;
	typedef typename TensorType::VectorSizeType VectorSizeType;
	typedef typename PsimagLite::Vector<TensorType*>::Type VectorTensorType;
	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef PsimagLite::Vector<PairSizeType>::Type VectorPairSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef std::pair<PsimagLite::String,SizeType> PairStringSizeType;
	typedef PsimagLite::Vector<PairStringSizeType>::Type VectorPairStringSizeType;
	typedef std::map<PairStringSizeType,SizeType> MapPairStringSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef TensorSrep TensorSrepType;

	SrepEquation(PsimagLite::String str,
	             const VectorTensorType& vt,
	             const VectorPairStringSizeType& tensorNameIds,
	             MapPairStringSizeType& nameIdsTensor)
	    : lhs_(0),rhs_(0),outputTensor_(0)
	{
		VectorStringType vstr;
		PsimagLite::tokenizer(str,vstr,"=");
		if (vstr.size() != 2)
			throw PsimagLite::RuntimeError("SrepEquation:: syntax error " + str + "\n");
		lhs_ = new TensorSrepType(vstr[0]);
		rhs_ = new TensorSrepType(vstr[1]);

		if (lhs_->size() != 1)
			PsimagLite::RuntimeError("SrepEquation:: LHS should have exactly 1 tensor\n");

		PairStringSizeType nameIdOfOutput(lhs_->operator ()(0).name(),
		                                  lhs_->operator ()(0).id());
		SizeType outputTensorIndex = nameIdsTensor[nameIdOfOutput];
		if (tensorNameIds[outputTensorIndex] != nameIdOfOutput) {
			PsimagLite::String msg("SrepEquation: Could not find ");
			msg += "output tensor " + nameIdOfOutput.first + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

		assert(outputTensorIndex < vt.size());
		outputTensor_ = vt[outputTensorIndex];
	}

	~SrepEquation()
	{
		delete lhs_;
		delete rhs_;
	}

	void fillOutput(const VectorSizeType& free,
	                ComplexOrRealType value)
	{
		outputTensor_->operator()(free) = value;
	}

	const TensorSrepType& lhs() const
	{
		assert(lhs_);
		return *lhs_;
	}

	const TensorSrepType& rhs() const
	{
		assert(rhs_);
		return *rhs_;
	}

	// FIXME: Gives away internals!
	TensorSrepType& rhs()
	{
		assert(rhs_);
		return *rhs_;
	}

	const TensorType& outputTensor() const
	{
		assert(outputTensor_);
		return *outputTensor_;
	}

	// FIXME: Gives away internals!
	TensorType& outputTensor()
	{
		assert(outputTensor_);
		return *outputTensor_;
	}

private:

	SrepEquation(const SrepEquation&);

	SrepEquation& operator=(const SrepEquation&);

	TensorSrepType* lhs_;
	TensorSrepType* rhs_;
	TensorType* outputTensor_;
}; // class SrepEquation

} // namespace Mera

#endif // SREPEQUATION_H
