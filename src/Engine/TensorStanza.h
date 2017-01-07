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

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	enum TensorTypeEnum {TENSOR_TYPE_UNKNOWN,
		                 TENSOR_TYPE_ERASED,
		                 TENSOR_TYPE_GP};

	enum IndexDirectionEnum {INDEX_DIR_IN, INDEX_DIR_OUT};

	enum IndexTypeEnum {INDEX_TYPE_SUMMED, INDEX_TYPE_FREE, INDEX_TYPE_DUMMY, INDEX_TYPE_DIM};

	explicit TensorStanza(PsimagLite::String srep)
	    : conjugate_(false),
	      srep_(srep),
	      name_(""),
	      type_(TENSOR_TYPE_GP),
	      maxSummed_(0),
	      maxFree_(0)
	{
		VectorStringType tokens;
		PsimagLite::tokenizer(srep_,tokens,"|");
		SizeType ts = tokens.size();
		if (ts == 0 || ts > 2) {
			PsimagLite::String str("TensorStanza: malformed stanza ");
			throw PsimagLite::RuntimeError(str + srep_ + "\n");
		}

		PsimagLite::String nameAndId = getNameFromToken(tokens[0]);
		cleanToken(tokens[0]);
		setArgVector(insSi_,tokens[0]);

		if (ts == 2) {
			cleanToken(tokens[1]);
			setArgVector(outsSi_,tokens[1]);
		}

		SizeType l = nameAndId.length();
		if (l == 0) {
			PsimagLite::String str("TensorStanza: malformed partial srep, empty name/id ");
			throw PsimagLite::RuntimeError(str + srep_ + "\n");
		}

		SizeType tmp = std::count(nameAndId.begin(),nameAndId.end(),'*');
		if (tmp > 1) {
			PsimagLite::String str("TensorStanza: too many * in stanza ");
			throw PsimagLite::RuntimeError(str + srep_ + "\n");
		}

		conjugate_ = (tmp == 1);

		std::size_t index = nameAndId.find("*");
		if (index != PsimagLite::String::npos) nameAndId.erase(index,1);

		l = nameAndId.length();
		assert(l > 0);
		index = nameAndId.find_first_of("0123456789");

		name_ = nameAndId.substr(0,index);

		if (index == PsimagLite::String::npos)
			id_ = 0;
		else
			id_ = atoi(nameAndId.substr(index,l-index).c_str());

		if (index == PsimagLite::String::npos) {
			PsimagLite::String str("TensorStanza: no digit for token ");
			throw PsimagLite::RuntimeError(str + nameAndId + "\n");
		}

		if (type_ == TENSOR_TYPE_UNKNOWN) {
			PsimagLite::String str("TensorStanza: unknown tensor type ");
			throw PsimagLite::RuntimeError(str + srep_ + "\n");
		}

		maxSummed_ = maxIndex('s');
		maxFree_ = maxIndex('f');
	}

	void conjugate()
	{
		conjugate_ = (!conjugate_);
		srep_ = srepFromObject();
	}

	void shiftSummedBy(SizeType ms)
	{
		if (ms == 0) return;
		SizeType ins = insSi_.size();
		for (SizeType i = 0; i < ins; ++i) {
			if (insSi_[i].first != 's') continue;
			insSi_[i].second += ms;
		}

		SizeType outs = outsSi_.size();
		for (SizeType i = 0; i < outs; ++i) {
			if (outsSi_[i].first != 's') continue;
			outsSi_[i].second += ms;
		}

		maxSummed_ += ms;
		srep_ = srepFromObject();
	}

	bool replaceSummedOrFrees(const VectorPairSizeType& replacements,
	                          char type)
	{
		if (type_ == TENSOR_TYPE_ERASED) return false;

		bool simplificationHappended = false;
		SizeType r = replacements.size();
		SizeType ins = insSi_.size();
		for (SizeType i = 0; i < ins; ++i) {
			if (insSi_[i].first != type)
				continue;

			SizeType s1 = insSi_[i].second;
			bool replace = false;
			for (SizeType j = 0; j < r; ++j) {
				if (s1 != replacements[j].first) continue;
				s1 = replacements[j].second;
				replace = true;
				break;
			}

			if (!replace) continue;

			simplificationHappended = true;
			insSi_[i].second = s1;
		}

		SizeType outs = outsSi_.size();
		for (SizeType i = 0; i < outs; ++i) {
			if (outsSi_[i].first != type)
				continue;

			SizeType s1 = outsSi_[i].second;
			bool replace = false;
			for (SizeType j = 0; j < r; ++j) {
				if (s1 != replacements[j].first) continue;
				s1 = replacements[j].second;
				replace = true;
				break;
			}

			if (!replace) continue;

			simplificationHappended = true;
			outsSi_[i].second = s1;
		}

		refresh();
		return simplificationHappended;
	}

	void contract(const VectorSizeType* indicesToContract)
	{
		if (indicesToContract && indicesToContract->size() == 0)
			return;

		SizeType ins = insSi_.size();

		for (SizeType i = 0; i < ins; ++i) {
			if (insSi_[i].first != 'f') continue;
			if (indicesToContract && std::find(indicesToContract->begin(),
			                                   indicesToContract->end(),
			                                   insSi_[i].second) == indicesToContract->end())
				continue;

			insSi_[i].first = 's';
		}

		SizeType outs = outsSi_.size();
		for (SizeType i = 0; i < outs; ++i) {
			if (outsSi_[i].first != 'f') continue;
			if (indicesToContract && std::find(indicesToContract->begin(),
			                                   indicesToContract->end(),
			                                   outsSi_[i].second) == indicesToContract->end())
				continue;

			outsSi_[i].first = 's';
		}

		maxSummed_ = maxIndex('s');
		maxFree_ = maxIndex('f');
		srep_ = srepFromObject();
	}

	void eraseTensor(VectorSizeType& s)
	{
		SizeType total = insSi_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (insSi_[i].first != 's') continue;
			s.push_back(insSi_[i].second);
		}

		total = outsSi_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (outsSi_[i].first != 's') continue;
			s.push_back(outsSi_[i].second);
		}

		setAsErased();
	}

	void setAsErased()
	{
		type_ = TENSOR_TYPE_ERASED;
		maxSummed_ = 0;
		maxFree_ = 0;
		insSi_.clear();
		outsSi_.clear();
		srep_ = "";
	}

	SizeType uncontract(const VectorSizeType& erased,
	                    SizeType count,
	                    VectorSizeType* mapping = 0)
	{
		if (type_ == TENSOR_TYPE_ERASED) return count;

		SizeType total = insSi_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (insSi_[i].first != 's') continue;
			if (std::find(erased.begin(),erased.end(),insSi_[i].second) ==
			        erased.end()) continue;
			insSi_[i].first = 'f';
			insSi_[i].second = count++;
			if (mapping)
				mapping->operator[](insSi_[i].second) = insSi_[i].second;
		}

		total = outsSi_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (outsSi_[i].first != 's') continue;
			if (std::find(erased.begin(),erased.end(),outsSi_[i].second) ==
			        erased.end()) continue;
			outsSi_[i].first = 'f';
			outsSi_[i].second = count++;
			if (mapping)
				mapping->operator[](outsSi_[i].second) = outsSi_[i].second;
		}

		maxSummed_ = maxIndex('s');
		maxFree_ = maxIndex('f');
		srep_ = srepFromObject();
		return count;
	}

	SizeType relabelFrees(SizeType count)
	{
		SizeType total = insSi_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (insSi_[i].first != 'f') continue;
			insSi_[i].second = count++;
		}

		total = outsSi_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (outsSi_[i].first != 'f') continue;
			outsSi_[i].second = count++;
		}

		maxFree_ = maxIndex('f');
		srep_ = srepFromObject();
		return count;
	}

	void setIndices(const VectorSizeType& summed, char c)
	{
		if (type_ == TENSOR_TYPE_ERASED) return;

		SizeType total = insSi_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (insSi_[i].first != c) continue;
			SizeType ind = insSi_[i].second;
			assert(ind < summed.size());
			insSi_[i].second = summed[ind];
		}

		total = outsSi_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (outsSi_[i].first != c) continue;
			SizeType ind = outsSi_[i].second;
			assert(ind < summed.size());
			outsSi_[i].second = summed[ind];
		}

		maxFree_ = maxIndex('f');
		maxSummed_ = maxIndex('s');
		srep_ = srepFromObject();
	}

	void loadSummedOrFree(VectorSizeType& summed,
	                      char c) const
	{
		if (type_ == TENSOR_TYPE_ERASED) return;

		SizeType total = insSi_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (insSi_[i].first != c) continue;
			summed.push_back(insSi_[i].second);
		}

		total = outsSi_.size();
		for (SizeType i = 0; i < total; ++i) {
			if (outsSi_[i].first != c) continue;
			summed.push_back(outsSi_[i].second);
		}
	}

	bool isConjugate() const { return conjugate_; }

	const PsimagLite::String& sRep() const { return srep_; }

	const PsimagLite::String& name() const { return name_; }

	SizeType id() const { return id_; }

	SizeType legs() const { return insSi_.size() + outsSi_.size(); }

	SizeType ins() const { return insSi_.size(); }

	SizeType outs() const { return outsSi_.size(); }

	IndexTypeEnum legType(SizeType ind) const
	{
		char c = '\0';
		SizeType ins = insSi_.size();
		if (ind < ins) {
			assert(ind < ins);
			c = insSi_[ind].first;
		}  else {
			assert(ind - ins < outsSi_.size());
			c = outsSi_[ind - ins].first;
		}

		if (c == 's') return INDEX_TYPE_SUMMED;
		if (c == 'f') return INDEX_TYPE_FREE;
		if (c == 'D') return INDEX_TYPE_DIM;
		return INDEX_TYPE_DUMMY;
	}

	char& legTypeChar(SizeType ind)
	{
		SizeType ins = insSi_.size();
		if (ind < ins) {
			assert(ind < insSi_.size());
			return insSi_[ind].first;
		}

		assert(ind - ins < outsSi_.size());
		return outsSi_[ind - ins].first;
	}

	const SizeType& legTag(SizeType ind) const
	{
		SizeType ins = insSi_.size();
		if (ind < ins) {
			assert(ind < insSi_.size());
			return insSi_[ind].second;
		}

		assert(ind - ins < outsSi_.size());
		return outsSi_[ind - ins].second;
	}

	SizeType& legTag(SizeType ind)
	{
		SizeType ins = insSi_.size();
		if (ind < ins) {
			assert(ind < insSi_.size());
			return insSi_[ind].second;
		}

		assert(ind - ins < outsSi_.size());
		return outsSi_[ind - ins].second;
	}

	TensorTypeEnum type() const { return type_; }

	const SizeType& maxTag(char c) const
	{
		return (c == 's') ? maxSummed_ : maxFree_;
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
		srep_ = srepFromObject();
		maxFree_ = maxIndex('f');
		maxSummed_ = maxIndex('s');
	}

//	friend std::ostream& operator<<(std::ostream& os, const TensorStanza& ts)
//	{
//		os<<"id= "<<ts.id_<<"\n";
//		os<<"conjugate= "<<ts.conjugate_<<"\n";
//		os<<"stanza= "<<ts.srep_<<"\n";
//		os<<"name= "<<ts.name_<<"\n";
//		os<<"type= "<<ts.type_<<"\n";
//		os<<"needs ins and outs FIXME\n";
//		return os;
//	}

private:

	PsimagLite::String srepFromObject() const
	{
		PsimagLite::String srep = name_;
		srep += ttos(id_);
		if (conjugate_) srep += "*";
		srep += "(";
		for (SizeType i = 0; i < insSi_.size(); ++i) {
			srep += insSi_[i].first + ttos(insSi_[i].second);
			if (i == insSi_.size() - 1) continue;
			srep += ",";
		}

		if (outsSi_.size() == 0) {
			srep += ")";
			return srep;
		}

		srep += "|";
		for (SizeType i = 0; i < outsSi_.size(); ++i) {
			srep += outsSi_[i].first + ttos(outsSi_[i].second);
			if (i == outsSi_.size() - 1) continue;
			srep += ",";
		}

		srep += ")";
		return srep;
	}

	void setArgVector(VectorPairCharSizeType& si, PsimagLite::String part) const
	{
		if (part.length() == 0 || part == ")") return;
		VectorStringType tokens;
		PsimagLite::tokenizer(part,tokens,",");
		SizeType total = tokens.size();
		for (SizeType i = 0; i < total; ++i)
			si.push_back(getPairCharInt(tokens[i]));
	}

	PairCharSizeType getPairCharInt(PsimagLite::String token) const
	{
		SizeType l = token.length();

		if (l == 0) {
			PsimagLite::String str("TensorStanza: malformed stanza, ");
			throw PsimagLite::RuntimeError(str + token + "\n");
		}

		std::size_t index = token.find_first_of("0123456789");

		if (index == 0) {
			// all numeric
			return PairCharSizeType('d',atoi(token.c_str()));
		}

		if (l == 1 && token != "d") {
			PsimagLite::String str("TensorStanza: malformed stanza, expecting a dummy, got ");
			throw PsimagLite::RuntimeError(str + token + "\n");
		}

		PsimagLite::String tmp = (l == 1) ? "0" : token.substr(1,l-1);
		return PairCharSizeType(token[0],atoi(tmp.c_str()));
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
		SizeType ins = insSi_.size();
		for (SizeType i = 0; i < ins; ++i) {
			if (insSi_[i].first != c) continue;
			if (max < insSi_[i].second) max = insSi_[i].second;
		}

		SizeType outs = outsSi_.size();
		for (SizeType i = 0; i < outs; ++i) {
			if (outsSi_[i].first != c) continue;
			if (max < outsSi_[i].second) max = outsSi_[i].second;
		}

		return max;
	}

	SizeType id_;
	bool conjugate_;
	PsimagLite::String srep_;
	PsimagLite::String name_;
	TensorTypeEnum type_;
	SizeType maxSummed_;
	SizeType maxFree_;
	VectorPairCharSizeType insSi_;
	VectorPairCharSizeType outsSi_;
};

} // namespace Mera
#endif // TENSORSTANZA_H
