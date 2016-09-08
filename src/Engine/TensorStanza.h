#ifndef TENSORSTANZA_H
#define TENSORSTANZA_H
#include "Vector.h"
#include "Tokenizer.h"
#include "TypeToString.h"
#include <cassert>
#include <iostream>

namespace Mera {

class TensorStanza {

	typedef std::pair<char,SizeType> PairCharSizeType;
	typedef PsimagLite::Vector<PairCharSizeType>::Type VectorPairCharSizeType;

public:

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	enum TensorTypeEnum {TENSOR_TYPE_UNKNOWN,TENSOR_TYPE_W,TENSOR_TYPE_U};

	enum IndexDirectionEnum {INDEX_DIR_IN, INDEX_DIR_OUT};

	enum IndexTypeEnum {INDEX_TYPE_SUMMED, INDEX_TYPE_FREE};

	TensorStanza(PsimagLite::String srep)
	    :  conjugate_(false),
	      srep_(srep),name_(""),
	      type_(TENSOR_TYPE_UNKNOWN)
	{
		VectorStringType tokens;
		PsimagLite::tokenizer(srep_,tokens,":");
		SizeType ts = tokens.size();
		if (ts < 6) {
			PsimagLite::String str("TensorStanza: malformed partial srep ");
			throw PsimagLite::RuntimeError(str + srep_ + "\n");
		}

		name_ = tokens[0];
		SizeType l = name_.length();
		if (l == 0) {
			PsimagLite::String str("TensorStanza: malformed partial srep, empty name ");
			throw PsimagLite::RuntimeError(str + srep_ + "\n");
		}

		if (name_.substr(l-1,l) == "*") conjugate_ = true;
		SizeType l2 = (conjugate_) ? l - 1 : l;
		if (name_.substr(0,l2) == "u") type_ = TENSOR_TYPE_U;
		if (name_.substr(0,l2) == "w") type_ = TENSOR_TYPE_W;

		if (type_ == TENSOR_TYPE_UNKNOWN) {
			PsimagLite::String str("TensorStanza: partial srep, tensor type");
			throw PsimagLite::RuntimeError(str + srep_ + "\n");
		}

		timeIndex_ = atoi(tokens[1].c_str());
		spaceIndex_ = atoi(tokens[2].c_str());
		SizeType ins = atoi(tokens[3].c_str());
		SizeType outs = atoi(tokens[4].c_str());

		if (ins + outs + 5 != ts) {
			PsimagLite::String str("TensorStanza: malformed partial srep ");
			throw PsimagLite::RuntimeError(str + srep_ + "\n");
		}

		for (SizeType i = 0; i < ins; ++i) {
			PsimagLite::String item = tokens[5+i];
			l = item.length();
			if (l == 0) {
				PsimagLite::String str("TensorStanza: malformed srep, in s-index ");
				throw PsimagLite::RuntimeError(str + ttos(i) + ", " + srep_ + "\n");
			}

			if (l == 1 && item != "d") {
				PsimagLite::String str("TensorStanza: malformed srep, in s-index ");
				throw PsimagLite::RuntimeError(str + ttos(i) + ", " + srep_ + "\n");
			}

			PsimagLite::String tmp = item.substr(1,l);
			insSi_.push_back(PairCharSizeType(item[0],atoi(tmp.c_str())));
		}

		for (SizeType i = 0; i < outs; ++i) {
			PsimagLite::String item = tokens[5+ins+i];
			l = item.length();

			if (l == 0) {
				PsimagLite::String str("TensorStanza: malformed srep, out s-index ");
				throw PsimagLite::RuntimeError(str + ttos(i) + ", " + srep_ + "\n");
			}

			if (l == 1 && item != "d") {
				PsimagLite::String str("TensorStanza: malformed srep, out s-index ");
				throw PsimagLite::RuntimeError(str + ttos(i) + ", " + srep_ + "\n");
			}

			PsimagLite::String tmp = (l == 1) ? "0" : item.substr(1,l);
			outsSi_.push_back(PairCharSizeType(item[0],atoi(tmp.c_str())));
		}
	}

	SizeType x() const { return spaceIndex_; }

	SizeType y() const { return timeIndex_; }

	SizeType ins() const { return insSi_.size(); }

	SizeType outs() const { return outsSi_.size(); }

	IndexTypeEnum indexType(SizeType ind, IndexDirectionEnum dir) const
	{
		char c = (dir == INDEX_DIR_IN) ? insSi_[ind].first : outsSi_[ind].first;
		return (c == 's') ? INDEX_TYPE_SUMMED : INDEX_TYPE_FREE;
	}

	TensorTypeEnum type() const { return type_; }

private:

	SizeType timeIndex_;
	SizeType spaceIndex_;
	bool conjugate_;
	PsimagLite::String srep_;
	PsimagLite::String name_;
	TensorTypeEnum type_;
	VectorPairCharSizeType insSi_;
	VectorPairCharSizeType outsSi_;
};

} // namespace Mera
#endif // TENSORSTANZA_H
