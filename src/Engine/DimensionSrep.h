#ifndef MERA_SREP_DIMENSION_H
#define MERA_SREP_DIMENSION_H
#include "TensorSrep.h"
#include "Vector.h"

namespace  Mera {

template<typename SymmetryLocalType>
class DimensionSrep {

	typedef TensorSrep TensorSrepType;
	typedef TensorSrepType::TensorStanzaType TensorStanzaType;
	typedef TensorSrepType::VectorSizeType VectorSizeType;

public:

	typedef typename SymmetryLocalType::VectorVectorSizeType VectorVectorSizeType;

	DimensionSrep(PsimagLite::String srep,
	              SymmetryLocalType& symmLocal,
	              SizeType m)
	    : m_(m),dsrep_(srep),symmLocal_(symmLocal)
	{
		dStoSymm(symmLocal_.qOne().size());
		alterFrees(symmLocal_.qOne().size());
		alterSummed();
	}

	const PsimagLite::String& operator()() const
	{
		return dsrep_.sRep();
	}

	const SymmetryLocalType& symmLocal() const { return symmLocal_; }

private:

	void alterFrees(SizeType d)
	{
		for (SizeType i = 0; i < dsrep_.size(); ++i) {

			TensorStanzaType ts = dsrep_(i);
			if (ts.type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;

			symmLocal_.setNameId(i, ts.name() + ttos(ts.id()));

			SizeType legs = ts.legs();
			for (SizeType j = 0; j < legs; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j);
				if (t == TensorStanzaType::INDEX_TYPE_FREE) {
					dsrep_.legTypeChar(i,j) = 'D';
					dsrep_.legTag(i,j) = d;
					symmLocal_.setQ(i,j);
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
			VectorVectorSizeType q(ins,0);
			SizeType counter = 0;

			for (SizeType j = 0; j < ins; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j);
				if (t != TensorStanzaType::INDEX_TYPE_DIM) continue;
				dim[j] = ts.legTag(j);
				q[j] = symmLocal_.q(i,j);

				counter++;
			}

			if (counter != ins) continue;

			SizeType outs = ts.outs();
			SizeType prodDim = SymmetryLocalType::truncateDimension(dim, m_);
			for (SizeType j = 0; j < outs; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j + ins);
				if (t != TensorStanzaType::INDEX_TYPE_SUMMED) continue;
				dsrep_.legTypeChar(i,j + ins) = 'D';
				SizeType s = dsrep_.legTag(i,j + ins);
				if (outs == 1) {
					dsrep_.legTag(i,j + ins) = prodDim;
					symmLocal_.setQ(i,j + ins, q, dim, m_);
				} else if (outs == dim.size()) {
					dsrep_.legTag(i,j + ins) = dim[j];
					symmLocal_.setQ(i,j + ins, q[j]);
				} else {
					throw PsimagLite::RuntimeError("DimensionSrep: outs > ins not supported\n");
				}

				replaceSummed(s,dsrep_.legTag(i,j + ins), symmLocal_.q(i,j + ins));
			}
		}

		dsrep_.refresh();
	}

	void replaceSummed(SizeType s, SizeType val, VectorSizeType* q)
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
				symmLocal_.setQ(i,j, q);
			}
		}
	}

	void dStoSymm(SizeType d)
	{
		for (SizeType i = 0; i < dsrep_.size(); ++i) {
			TensorStanzaType ts = dsrep_(i);
			if (ts.type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;

			SizeType legs = ts.legs();
			for (SizeType j = 0; j < legs; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j);
				if (t != TensorStanzaType::INDEX_TYPE_DIM) continue;
				if (ts.legTag(j) != d) continue;
				symmLocal_.setQ(i,j);
			}
		}
	}

	SizeType m_;
	TensorSrepType dsrep_;
	SymmetryLocalType& symmLocal_;
}; //

} // namespace  Mera
#endif // MERA_SREP_DIMENSION_H
