#ifndef MERA_SREP_DIMENSION_H
#define MERA_SREP_DIMENSION_H
#include "TensorSrep.h"
#include "Vector.h"

namespace  Mera {

class DimensionSrep {

	typedef TensorSrep TensorSrepType;
	typedef TensorSrepType::TensorStanzaType TensorStanzaType;

public:

	DimensionSrep(PsimagLite::String srep, SizeType d)
	    :srep_(srep),dsrep_(srep)
	{
		alterFrees();
	}

	const PsimagLite::String& operator()() const
	{
		return dsrep_.sRep();
	}

private:

	void alterFrees()
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
				}
			}

			SizeType outs = ts.outs();
			for (SizeType j = 0; j < outs; ++j) {
				TensorStanzaType::IndexTypeEnum t = ts.legType(j,out);
				if (t == TensorStanzaType::INDEX_TYPE_FREE) {
					dsrep_.legTypeChar(i,j,out) = 'D';
				}
			}
		}

		dsrep_.refresh();
	}

	TensorSrepType srep_;
	TensorSrepType dsrep_;
}; //

} // namespace  Mera
#endif // MERA_SREP_DIMENSION_H
