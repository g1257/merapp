#ifndef TENSORSREP_H
#define TENSORSREP_H
#include "Vector.h"
#include "TensorStanza.h"

namespace Mera {

class TensorSrep {

	typedef PsimagLite::Vector<TensorStanza*>::Type VectorTensorStanza;

public:

	TensorSrep(PsimagLite::String srep)
	    : srep_(srep)
	{
		cleanWhiteSpace(srep_);
		parseIt();
	}

	~TensorSrep()
	{
		for (SizeType i = 0; i < data_.size(); ++i) {
			delete data_[i];
			data_[i] = 0;
		}
	}

	SizeType size() const { return data_.size(); }

	const TensorStanza& operator()(SizeType ind) const
	{
		assert(ind < data_.size());
		return *data_[ind];
	}

private:

	void cleanWhiteSpace(PsimagLite::String& srep) const
	{
		SizeType l = srep.length();
		PsimagLite::String tmp("");

		for (SizeType i = 0; i < l; ++i) {
			char c = srep[i];
			if (isWhiteSpace(c)) continue;
			tmp += c;
		}

		srep = tmp;
	}

	bool isWhiteSpace(unsigned char c) const
	{
		return (c == ' ' || c == '\t' || c=='\n');
	}

	void parseIt()
	{
		SizeType l = srep_.length();
		SizeType counter = 0;
		SizeType loc = 0;
		while (loc < l) {
			std::size_t index = srep_.find(")",loc,l);
			if (index == PsimagLite::String::npos) {
				PsimagLite::String str("TensorSrep: stanza without closing brace?!\n");
				throw PsimagLite::RuntimeError(str + " at offset " + ttos(loc) + "\n");
			}

			PsimagLite::String stanza = srep_.substr(loc,index+1);
			data_[counter++] = new TensorStanza(stanza);
			loc = index + 1;
		}
	}

	PsimagLite::String srep_;
	VectorTensorStanza data_;
}; // class TensorSrep

} // namespace Mera
#endif // TENSORSREP_H
