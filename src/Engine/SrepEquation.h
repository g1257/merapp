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
		lhs_ = new TensorStanza(vstr[0]);
		rhs_ = new TensorSrepType(vstr[1]);


		PairStringSizeType nameIdOfOutput(lhs_->name(), lhs_->id());
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

	void canonicalize() const
	{
		VectorPairSizeType frees;
		computeFrees(frees);
		rhs_->simplifyFrees(frees);
	}

	void fillOutput(const VectorSizeType& free,
	                ComplexOrRealType value)
	{
		outputTensor_->operator()(free) = value;
	}

	const TensorStanza& lhs() const
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

	void computeFrees(VectorPairSizeType& replacements) const
	{
		TensorStanza::IndexDirectionEnum in = TensorStanza::INDEX_DIR_IN;
		TensorStanza::IndexDirectionEnum out = TensorStanza::INDEX_DIR_OUT;

		SizeType legs = lhs_->ins();
		SizeType counter = 0;
		for (SizeType j = 0; j < legs; ++j) {
			if (lhs_->legType(j,in) != TensorStanza::INDEX_TYPE_FREE)
				continue;
			replacements.push_back(PairSizeType(lhs_->legTag(j,in),counter++));
		}

		legs = lhs_->outs();
		for (SizeType j = 0; j < legs; ++j) {
			if (lhs_->legType(j,out) != TensorStanza::INDEX_TYPE_FREE)
				continue;
			replacements.push_back(PairSizeType(lhs_->legTag(j,out),counter++));
		}
	}

	SrepEquation(const SrepEquation&);

	SrepEquation& operator=(const SrepEquation&);

	TensorStanza* lhs_;
	TensorSrepType* rhs_;
	TensorType* outputTensor_;
}; // class SrepEquation

} // namespace Mera

#endif // SREPEQUATION_H
