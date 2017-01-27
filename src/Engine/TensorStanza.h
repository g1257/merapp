/*
Copyright (c) 2016, UT-Battelle, LLC

MERA++, Version 0.

This file is part of MERA++.
MERA++ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
MERA++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with MERA++. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef TENSORSTANZA_H
#define TENSORSTANZA_H
#include "Vector.h"
#include "Tokenizer.h"
#include "TypeToString.h"
#include "TensorLeg.h"
#include <cassert>
#include <iostream>
#include <algorithm>

namespace Mera {

class TensorStanza {

	typedef std::pair<char,SizeType> PairCharSizeType;
	typedef PsimagLite::Vector<PairCharSizeType>::Type VectorPairCharSizeType;
	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef PsimagLite::Vector<PairSizeType>::Type VectorPairSizeType;

public:

	typedef TensorLeg TensorLegType;
	typedef PsimagLite::Vector<TensorLegType>::Type VectorTensorLegType;
	typedef TensorLegType::IndexDirectionEnum IndexDirectionEnum;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	static const IndexDirectionEnum INDEX_DIR_IN = TensorLegType::INDEX_DIR_IN;
	static const IndexDirectionEnum INDEX_DIR_OUT = TensorLegType::INDEX_DIR_OUT;

	enum TensorTypeEnum {TENSOR_TYPE_UNKNOWN,
		                 TENSOR_TYPE_ERASED,
		                 TENSOR_TYPE_GP};

	enum IndexTypeEnum {INDEX_TYPE_SUMMED, INDEX_TYPE_FREE, INDEX_TYPE_DUMMY, INDEX_TYPE_DIM};

	struct NonPointerOpaque {
		NonPointerOpaque(PsimagLite::String srep)
		    : id_(0),
		      conjugate_(false),
		      srep_(srep),
		      name_(""),
		      type_(TENSOR_TYPE_GP),
		      maxSummed_(0),
		      maxFree_(0)
		{}

		SizeType id_;
		bool conjugate_;
		PsimagLite::String srep_;
		PsimagLite::String name_;
		TensorTypeEnum type_;
		SizeType maxSummed_;
		SizeType maxFree_;
	};

	explicit TensorStanza(PsimagLite::String srep)
	    : opaque_(srep)
	{
		VectorStringType tokens;
		PsimagLite::tokenizer(opaque_.srep_,tokens,"|");
		SizeType ts = tokens.size();
		if (ts == 0 || ts > 2) {
			PsimagLite::String str("TensorStanza: malformed stanza ");
			throw PsimagLite::RuntimeError(str + opaque_.srep_ + "\n");
		}

		PsimagLite::String nameAndId = getNameFromToken(tokens[0]);
		cleanToken(tokens[0]);
		setArgVector(legs_, tokens[0], INDEX_DIR_IN);

		if (ts == 2) {
			cleanToken(tokens[1]);
			setArgVector(legs_,tokens[1], INDEX_DIR_OUT);
		}

		SizeType l = nameAndId.length();
		if (l == 0) {
			PsimagLite::String str("TensorStanza: malformed partial srep, empty name/id ");
			throw PsimagLite::RuntimeError(str + opaque_.srep_ + "\n");
		}

		SizeType tmp = std::count(nameAndId.begin(),nameAndId.end(),'*');
		if (tmp > 1) {
			PsimagLite::String str("TensorStanza: too many * in stanza ");
			throw PsimagLite::RuntimeError(str + opaque_.srep_ + "\n");
		}

		opaque_.conjugate_ = (tmp == 1);

		std::size_t index = nameAndId.find("*");
		if (index != PsimagLite::String::npos) nameAndId.erase(index,1);

		l = nameAndId.length();
		assert(l > 0);
		index = nameAndId.find_first_of("0123456789");

		opaque_.name_ = nameAndId.substr(0,index);

		if (index == PsimagLite::String::npos)
			opaque_.id_ = 0;
		else
			opaque_.id_ = atoi(nameAndId.substr(index,l-index).c_str());

		if (index == PsimagLite::String::npos) {
			PsimagLite::String str("TensorStanza: no digit for token ");
			throw PsimagLite::RuntimeError(str + nameAndId + "\n");
		}

		if (opaque_.type_ == TENSOR_TYPE_UNKNOWN) {
			PsimagLite::String str("TensorStanza: unknown tensor type ");
			throw PsimagLite::RuntimeError(str + opaque_.srep_ + "\n");
		}

		opaque_.maxSummed_ = maxIndex('s');
		opaque_.maxFree_ = maxIndex('f');
	}

	void conjugate()
	{
		opaque_.conjugate_ = (!opaque_.conjugate_);
		opaque_.srep_ = srepFromObject();
	}

	void shiftSummedBy(SizeType ms)
	{
		if (ms == 0) return;
		SizeType legs = legs_.size();
		for (SizeType i = 0; i < legs; ++i) {
			if (legs_[i].name() != 's') continue;
			legs_[i].numericTag() += ms;
		}

		opaque_.maxSummed_ += ms;
		opaque_.srep_ = srepFromObject();
	}

	bool replaceSummedOrFrees(const VectorPairSizeType& replacements,
	                          char type)
	{
		if (opaque_.type_ == TENSOR_TYPE_ERASED) return false;

		bool simplificationHappended = false;
		SizeType r = replacements.size();
		SizeType legs = legs_.size();
		for (SizeType i = 0; i < legs; ++i) {
			if (legs_[i].name() != type)
				continue;

			SizeType s1 = legs_[i].numericTag();
			bool replace = false;
			for (SizeType j = 0; j < r; ++j) {
				if (s1 != replacements[j].first) continue;
				s1 = replacements[j].second;
				replace = true;
				break;
			}

			if (!replace) continue;

			simplificationHappended = true;
			legs_[i].numericTag() = s1;
		}

		refresh();
		return simplificationHappended;
	}

	void contract(const VectorSizeType* indicesToContract)
	{
		if (indicesToContract && indicesToContract->size() == 0)
			return;

		SizeType legs = legs_.size();

		for (SizeType i = 0; i < legs; ++i) {
			if (legs_[i].name() != 'f') continue;
			if (indicesToContract && std::find(indicesToContract->begin(),
			                                   indicesToContract->end(),
			                                   legs_[i].numericTag()) == indicesToContract->end())
				continue;

			legs_[i].name() = 's';
		}

		opaque_.maxSummed_ = maxIndex('s');
		opaque_.maxFree_ = maxIndex('f');
		opaque_.srep_ = srepFromObject();
	}

	void eraseTensor(VectorSizeType& s)
	{
		SizeType total = legs_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (legs_[i].name() != 's') continue;
			s.push_back(legs_[i].numericTag());
		}

		setAsErased();
	}

	void setAsErased()
	{
		opaque_.type_ = TENSOR_TYPE_ERASED;
		opaque_.maxSummed_ = 0;
		opaque_.maxFree_ = 0;
		legs_.clear();
		opaque_.srep_ = "";
	}

	SizeType uncontract(const VectorSizeType& erased,
	                    SizeType count,
	                    VectorSizeType* mapping = 0)
	{
		if (opaque_.type_ == TENSOR_TYPE_ERASED) return count;

		SizeType total = legs_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (legs_[i].name() != 's') continue;
			if (std::find(erased.begin(),erased.end(),legs_[i].numericTag()) ==
			        erased.end()) continue;
			legs_[i].name() = 'f';
			legs_[i].numericTag() = count++;
			if (mapping)
				mapping->operator[](legs_[i].numericTag()) = legs_[i].numericTag();
		}

		opaque_.maxSummed_ = maxIndex('s');
		opaque_.maxFree_ = maxIndex('f');
		opaque_.srep_ = srepFromObject();
		return count;
	}

	SizeType relabelFrees(SizeType count)
	{
		SizeType total = legs_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (legs_[i].name() != 'f') continue;
			legs_[i].numericTag() = count++;
		}

		opaque_.maxFree_ = maxIndex('f');
		opaque_.srep_ = srepFromObject();
		return count;
	}

	void setIndices(const VectorSizeType& summed, char c)
	{
		if (opaque_.type_ == TENSOR_TYPE_ERASED) return;

		SizeType total = legs_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (legs_[i].name() != c) continue;
			SizeType ind = legs_[i].numericTag();
			assert(ind < summed.size());
			legs_[i].numericTag() = summed[ind];
		}

		opaque_.maxFree_ = maxIndex('f');
		opaque_.maxSummed_ = maxIndex('s');
		opaque_.srep_ = srepFromObject();
	}

	void loadSummedOrFree(VectorSizeType& summed,
	                      char c) const
	{
		if (opaque_.type_ == TENSOR_TYPE_ERASED) return;

		SizeType total = legs_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (legs_[i].name() != c) continue;
			summed.push_back(legs_[i].numericTag());
		}
	}

	bool isConjugate() const { return opaque_.conjugate_; }

	const PsimagLite::String& sRep() const { return opaque_.srep_; }

	const PsimagLite::String& name() const { return opaque_.name_; }

	SizeType id() const { return opaque_.id_; }

	SizeType legs() const { return legs_.size(); }

	SizeType ins() const { return countLegsWithDir(INDEX_DIR_IN); }

	SizeType outs() const { return countLegsWithDir(INDEX_DIR_OUT); }

	IndexTypeEnum legType(SizeType ind) const
	{
		assert(ind < legs_.size());
		char c = legs_[ind].name();

		if (c == 's') return INDEX_TYPE_SUMMED;
		if (c == 'f') return INDEX_TYPE_FREE;
		if (c == 'D') return INDEX_TYPE_DIM;
		return INDEX_TYPE_DUMMY;
	}

	char& legTypeChar(SizeType ind)
	{
		assert(ind < legs_.size());
		return legs_[ind].name();
	}

	const SizeType& legTag(SizeType ind) const
	{
		assert(ind < legs_.size());
		return legs_[ind].numericTag();
	}

	SizeType& legTag(SizeType ind)
	{
		assert(ind < legs_.size());
		return legs_[ind].numericTag();
	}

	TensorTypeEnum type() const { return opaque_.type_; }

	const SizeType& maxTag(char c) const
	{
		return (c == 's') ? opaque_.maxSummed_ : opaque_.maxFree_;
	}

	static PsimagLite::String indexTypeToString(IndexTypeEnum t)
	{
		if (t == INDEX_TYPE_SUMMED) return "s";
		if (t == INDEX_TYPE_FREE) return "f";
		if (t == INDEX_TYPE_DIM) return "D";
		return "d";
	}

	void refresh()
	{
		opaque_.srep_ = srepFromObject();
		opaque_.maxFree_ = maxIndex('f');
		opaque_.maxSummed_ = maxIndex('s');
	}

	bool hasLegType(char c) const
	{
		SizeType total = legs_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (legs_[i].name() != c) continue;
			return true;
		}

		return false;
	}

	static TensorStanza* newStanza(const TensorStanza& other)
	{
		TensorStanza* intercept = new TensorStanza(other);
		return intercept;
	}

private:

	SizeType countLegsWithDir(IndexDirectionEnum dir) const
	{
		SizeType count = 0;
		for (SizeType i = 0; i < legs_.size(); ++i) {
			if (legs_[i].dir() == dir) ++count;
		}

		return count;
	}

	PsimagLite::String srepFromObject() const
	{
		PsimagLite::String srep = opaque_.name_;
		srep += ttos(opaque_.id_);
		if (opaque_.conjugate_) srep += "*";
		srep += "(";
		SizeType nins = countLegsWithDir(INDEX_DIR_IN);
		SizeType nouts = countLegsWithDir(INDEX_DIR_OUT);

		for (SizeType i = 0; i < nins; ++i) {
			assert(legs_[i].dir() == INDEX_DIR_IN);
			srep += legs_[i].name() + ttos(legs_[i].numericTag());
			if (i == nins - 1) continue;
			srep += ",";
		}

		if (nouts == 0) {
			srep += ")";
			return srep;
		}

		srep += "|";
		for (SizeType i = 0; i < nouts; ++i) {
			assert(i + nins < legs_.size());
			if (legs_[i+nins].dir() != INDEX_DIR_OUT) continue;
			srep += legs_[i+nins].name() + ttos(legs_[i+nins].numericTag());
			if (i == nouts - 1) continue;
			srep += ",";
		}

		srep += ")";
		return srep;
	}

	void setArgVector(VectorTensorLegType& si,
	                  PsimagLite::String part,
	                  IndexDirectionEnum inOrOut) const
	{
		if (part.length() == 0 || part == ")") return;
		VectorStringType tokens;
		PsimagLite::tokenizer(part,tokens,",");
		SizeType total = tokens.size();
		for (SizeType i = 0; i < total; ++i)
			si.push_back(TensorLegType(tokens[i],inOrOut));
	}

	PsimagLite::String getNameFromToken(PsimagLite::String token) const
	{
		std::size_t index = token.find("(");
		return token.substr(0,index);
	}

	void cleanToken(PsimagLite::String& token) const
	{
		SizeType l = token.length();
		assert(l > 0);
		std::size_t index = token.find("(");
		PsimagLite::String tmp = token;
		if (index != PsimagLite::String::npos)
			tmp = token.substr(index+1,token.length());

		if (token[l-1] == ')') tmp.substr(0,l-1);

		token = tmp;
	}

	SizeType maxIndex(char c) const
	{
		SizeType max = 0;
		SizeType legs = legs_.size();
		for (SizeType i = 0; i < legs; ++i) {
			if (legs_[i].name() != c) continue;
			if (max < legs_[i].numericTag()) max = legs_[i].numericTag();
		}

		return max;
	}

	NonPointerOpaque opaque_;
	VectorTensorLegType legs_;
};

} // namespace Mera
#endif // TENSORSTANZA_H
