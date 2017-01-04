#ifndef MERA_SREP_DIMENSION_H
#define MERA_SREP_DIMENSION_H
#include "TensorSrep.h"
#include "Vector.h"

namespace  Mera {

class DimensionSrep {

	typedef TensorSrep TensorSrepType;
	typedef TensorSrepType::TensorStanzaType TensorStanzaType;
	typedef TensorSrepType::VectorSizeType VectorSizeType;

public:

	DimensionSrep(PsimagLite::String srep, SizeType d, SizeType m)
	    :srep_(srep),m_(m),dsrep_(srep)
	{
		alterFrees(d);
		alterSummed();
	}

	const PsimagLite::String& operator()() const
	{
		return dsrep_.sRep();
	}

private:

	void alterFrees(SizeType d)
	{
		for (SizeType i = 0; i < dsrep_.size(); ++i) {
			TensorStanzaType ts = dsrep_(i);
			if (ts.type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;

			SizeType legs = ts.legs();
			for (SizeType j = 0; j < legs; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j);
				if (t == TensorStanzaType::INDEX_TYPE_FREE) {
					dsrep_.legTypeChar(i,j) = 'D';
					dsrep_.legTag(i,j) = d;
				}
			}
		}

		dsrep_.refresh();
	}

	void alterSummed()
	{
		for (SizeType i = 0; i < dsrep_.size(); ++i) {
			TensorStanzaType ts = dsrep_(i);
			if (ts.type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;

			SizeType ins = ts.ins();
			VectorSizeType dim(ins,0);
			SizeType counter = 0;
			for (SizeType j = 0; j < ins; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j);
				if (t != TensorStanzaType::INDEX_TYPE_DIM) continue;
				dim[j] = ts.legTag(j);
				counter++;
			}

			if (counter != ins) continue;

			SizeType outs = ts.outs();
			for (SizeType j = 0; j < outs; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j + ins);
				if (t != TensorStanzaType::INDEX_TYPE_SUMMED) continue;
				dsrep_.legTypeChar(i,j + ins) = 'D';
				SizeType s = dsrep_.legTag(i,j + ins);
				if (outs == 1) {
					dsrep_.legTag(i,j + ins) = truncateDimension(productOf(dim));
				} else if (outs == dim.size()) {
					dsrep_.legTag(i,j + ins) = dim[j];
				} else {
					throw PsimagLite::RuntimeError("DimensionSrep: outs > ins not supported\n");
				}

				replaceSummed(s,dsrep_.legTag(i,j + ins));
			}
		}

		dsrep_.refresh();
	}

	SizeType productOf(const VectorSizeType& dim) const
	{
		SizeType n = dim.size();
		if (n == 0) return 0;
		SizeType ret = dim[0];
		for (SizeType i = 1; i < n; ++i)
			ret *= dim[i];

		return ret;
	}

	void replaceSummed(SizeType s, SizeType val)
	{
		for (SizeType i = 0; i < dsrep_.size(); ++i) {
			TensorStanzaType ts = dsrep_(i);
			if (ts.type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;

			SizeType legs = ts.legs();
			for (SizeType j = 0; j < legs; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j);
				if (t != TensorStanzaType::INDEX_TYPE_SUMMED) continue;
				if (ts.legTag(j) != s) continue;
				dsrep_.legTypeChar(i,j) = 'D';
				dsrep_.legTag(i,j) = val;
			}
		}
	}

	SizeType truncateDimension(SizeType x) const
	{
		if (m_ == 0) return x;
		return std::min(m_,x);
	}

	TensorSrepType srep_;
	SizeType m_;
	TensorSrepType dsrep_;
}; //

} // namespace  Mera
#endif // MERA_SREP_DIMENSION_H
