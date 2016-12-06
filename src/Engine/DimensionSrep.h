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
		TensorStanzaType::IndexDirectionEnum in = TensorStanzaType::INDEX_DIR_IN;
		TensorStanzaType::IndexDirectionEnum out = TensorStanzaType::INDEX_DIR_OUT;

		for (SizeType i = 0; i < dsrep_.size(); ++i) {
			TensorStanzaType ts = dsrep_(i);
			if (ts.type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;

			SizeType ins = ts.ins();
			for (SizeType j = 0; j < ins; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j,in);
				if (t == TensorStanzaType::INDEX_TYPE_FREE) {
					dsrep_.legTypeChar(i,j,in) = 'D';
					dsrep_.legTag(i,j,in) = d;
				}
			}

			SizeType outs = ts.outs();
			for (SizeType j = 0; j < outs; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j,out);
				if (t == TensorStanzaType::INDEX_TYPE_FREE) {
					dsrep_.legTypeChar(i,j,out) = 'D';
					dsrep_.legTag(i,j,out) = d;
				}
			}
		}

		dsrep_.refresh();
	}

	void alterSummed()
	{
		TensorStanzaType::IndexDirectionEnum in = TensorStanzaType::INDEX_DIR_IN;
		TensorStanzaType::IndexDirectionEnum out = TensorStanzaType::INDEX_DIR_OUT;

		for (SizeType i = 0; i < dsrep_.size(); ++i) {
			TensorStanzaType ts = dsrep_(i);
			if (ts.type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;

			SizeType ins = ts.ins();
			VectorSizeType dim(ins,0);
			SizeType counter = 0;
			for (SizeType j = 0; j < ins; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j,in);
				if (t != TensorStanzaType::INDEX_TYPE_DIM) continue;
				dim[j] = ts.legTag(j,in);
				counter++;
			}

			if (counter != ins) continue;

			SizeType outs = ts.outs();
			for (SizeType j = 0; j < outs; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j,out);
				if (t != TensorStanzaType::INDEX_TYPE_SUMMED) continue;
				dsrep_.legTypeChar(i,j,out) = 'D';
				SizeType s = dsrep_.legTag(i,j,out);
				if (outs == 1) {
					dsrep_.legTag(i,j,out) = truncateDimension(productOf(dim));
				} else if (outs == dim.size()) {
					dsrep_.legTag(i,j,out) = dim[j];
				} else {
					throw PsimagLite::RuntimeError("DimensionSrep: outs > ins not supported\n");
				}

				replaceSummed(s,dsrep_.legTag(i,j,out));
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
		TensorStanzaType::IndexDirectionEnum in = TensorStanzaType::INDEX_DIR_IN;
		TensorStanzaType::IndexDirectionEnum out = TensorStanzaType::INDEX_DIR_OUT;

		for (SizeType i = 0; i < dsrep_.size(); ++i) {
			TensorStanzaType ts = dsrep_(i);
			if (ts.type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;

			SizeType ins = ts.ins();
			for (SizeType j = 0; j < ins; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j,in);
				if (t != TensorStanzaType::INDEX_TYPE_SUMMED) continue;
				if (ts.legTag(j,in) != s) continue;
				dsrep_.legTypeChar(i,j,in) = 'D';
				dsrep_.legTag(i,j,in) = val;
			}

			SizeType outs = ts.outs();
			for (SizeType j = 0; j < outs; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j,out);
				if (t != TensorStanzaType::INDEX_TYPE_SUMMED) continue;
				if (ts.legTag(j,out) != s) continue;
				dsrep_.legTypeChar(i,j,out) = 'D';
				dsrep_.legTag(i,j,out) = val;
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
