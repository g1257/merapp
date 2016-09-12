#ifndef TENSORSREP_H
#define TENSORSREP_H
#include "Vector.h"
#include "TensorStanza.h"
#include <algorithm>

namespace Mera {

class TensorSrep {

	typedef PsimagLite::Vector<TensorStanza*>::Type VectorTensorStanza;

public:

	typedef TensorStanza::VectorSizeType VectorSizeType;
	typedef TensorStanza TensorStanzaType;

	TensorSrep(PsimagLite::String srep)
	    : srep_(srep)
	{
		cleanWhiteSpace(srep_);
		parseIt();
	}

	TensorSrep(const TensorSrep& other)
	    : srep_(other.srep_),data_(other.data_.size())
	{
		for (SizeType i = 0; i < data_.size(); ++i)
			data_[i] = new TensorStanzaType(*(other.data_[i]));
	}

	~TensorSrep()
	{
		for (SizeType i = 0; i < data_.size(); ++i) {
			delete data_[i];
			data_[i] = 0;
		}
	}

	void shiftSummedBy(SizeType ms)
	{
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			data_[i]->shiftSummedBy(ms);
			srep_ += data_[i]->sRep();
		}
	}

	// FIXME: should it be a member?
	void contract(const TensorSrep& other,
	              const VectorSizeType& indicesToContract)
	{
		TensorSrep copy(other);
		SizeType ms = maxSummed();
		copy.shiftSummedBy(ms + 1);
		append(copy);
		contract(indicesToContract);
	}

	const PsimagLite::String& sRep() const { return srep_; }

	void conjugate()
	{
		srep_ = "";
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			data_[i]->conjugate();
			srep_ += data_[i]->sRep();
		}
	}

	SizeType size() const { return data_.size(); }

	const TensorStanza& operator()(SizeType ind) const
	{
		assert(ind < data_.size());
		return *data_[ind];
	}

	friend std::ostream& operator<<(std::ostream& os, const TensorSrep& ts)
	{
		os<<"tensorSrep.size="<<ts.size()<<"\n";
		for (SizeType i = 0; i < ts.size(); ++i) {
			os<<*(ts.data_[i])<<"\n";
		}

		return os;
	}

private:

	TensorSrep& operator=(const TensorSrep&);

	void parseIt()
	{
		SizeType l = srep_.length();
		SizeType counter = 0;
		SizeType loc = 0;
		SizeType parensOpen = std::count(srep_.begin(),srep_.end(),'(');

		SizeType parensClosed = std::count(srep_.begin(),srep_.end(),')');

		if (parensOpen != parensClosed) {
			PsimagLite::String str("TensorSrep: unbalanced parens\n");
			throw PsimagLite::RuntimeError(str + " at offset " + ttos(loc) + "\n");
		}

		data_.resize(parensOpen,0);
		while (loc < l) {
			std::size_t index = srep_.find(")",loc);
			if (index == PsimagLite::String::npos) {
				PsimagLite::String str("TensorSrep: stanza without closing brace?!\n");
				throw PsimagLite::RuntimeError(str);
			}

			PsimagLite::String stanza = srep_.substr(loc,index + 1 - loc);
			assert(counter < data_.size());
			data_[counter++] = new TensorStanza(stanza);
			loc = index + 1;
		}
	}

	void contract(const VectorSizeType& indicesToContract)
	{
		if (indicesToContract.size() == 0) return;
		SizeType ms = 1 + maxSummed();
		srep_ = "";
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			data_[i]->contract(indicesToContract,ms);
			srep_ += data_[i]->sRep();
		}
	}

	void append(const TensorSrep& other)
	{
		SizeType ntensors = data_.size();
		SizeType add = other.size();
		data_.resize(ntensors + add);
		for (SizeType i = 0; i < add; ++i) {
			data_[ntensors + i] = new TensorStanza(*other.data_[i]);
			srep_ += data_[ntensors + i]->sRep();
		}
	}

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

	SizeType maxSummed() const
	{
		SizeType ntensors = data_.size();
		SizeType max = 0;
		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType tmp = data_[i]->maxSummed();
			if (max < tmp) max = tmp;
		}

		return max;
	}

	PsimagLite::String srep_;
	VectorTensorStanza data_;
}; // class TensorSrep

} // namespace Mera
#endif // TENSORSREP_H
