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

	TensorSrep& operator=(const TensorSrep& other)
	{
		if (this == &other) return *this;

		for (SizeType i = 0; i < data_.size(); ++i)
			delete data_[i];

		srep_ = other.srep_;

		data_.clear();
		data_.resize(other.data_.size());
		for (SizeType i = 0; i < data_.size(); ++i)
			data_[i] = new TensorStanzaType(*(other.data_[i]));

		return *this;
	}

	~TensorSrep()
	{
		for (SizeType i = 0; i < data_.size(); ++i) {
			delete data_[i];
			data_[i] = 0;
		}
	}

	void contract(const TensorSrep& other)
	{
		TensorSrep copy(other);
		SizeType ms = maxTag('s');
		copy.shiftSummedBy(ms + 1);
		append(copy);
		contract(0);
	}

	// FIXME: should it be a member?
	void contract(const TensorSrep& other,
	              const VectorSizeType& indicesToContract)
	{
		TensorSrep copy(other);
		SizeType ms = maxTag('s');
		copy.shiftSummedBy(ms + 1);
		append(copy);
		contract(&indicesToContract);
		relabelFrees(size() - copy.size());
	}

	// FIXME: Erase tensor by name/id not index
	void eraseTensor(SizeType index)
	{
		assert(index < data_.size());
		VectorSizeType sErased;
		data_[index]->eraseTensor(sErased);
		std::cerr<<sErased;
		SizeType ntensors = data_.size();
		SizeType count = maxTag('f') + 1;
		srep_ = "";
		for (SizeType i = 0; i < ntensors; ++i) {
			count = data_[i]->uncontract(sErased,count);
			srep_ += data_[i]->sRep();
		}
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

	bool isValid(bool verbose) const
	{
		TensorStanzaType::IndexDirectionEnum in = TensorStanzaType::INDEX_DIR_IN;
		TensorStanzaType::IndexDirectionEnum out = TensorStanzaType::INDEX_DIR_OUT;

		VectorSizeType summedIndices(maxTag('s')+1,0);
		VectorSizeType freeIndices(maxTag('f')+1,0);

		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType legs = data_[i]->ins();
			for (SizeType j = 0; j < legs; ++j) {
				if (data_[i]->legType(j,in) == TensorStanzaType::INDEX_TYPE_SUMMED)
					summedIndices[data_[i]->legTag(j,in)]++;
				else if (data_[i]->legType(j,in) == TensorStanzaType::INDEX_TYPE_FREE)
					freeIndices[data_[i]->legTag(j,in)]++;
			}

			legs = data_[i]->outs();
			for (SizeType j = 0; j < legs; ++j) {
				if (data_[i]->legType(j,out) == TensorStanzaType::INDEX_TYPE_SUMMED)
					summedIndices[data_[i]->legTag(j,out)]++;
				else if (data_[i]->legType(j,out) == TensorStanzaType::INDEX_TYPE_FREE)
					freeIndices[data_[i]->legTag(j,out)]++;
			}
		}

		bool b = shouldAppear(summedIndices,2,"s",verbose);
		if (!verbose && !b) return false;

		b = shouldAppear(freeIndices,1,"f",verbose);
		return b;
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

	void contract(const VectorSizeType* indicesToContract)
	{
		if (indicesToContract && indicesToContract->size() == 0)
			return;

		SizeType ms = 1 + maxTag('s');
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

	void relabelFrees(SizeType start)
	{
		SizeType ntensors = data_.size();
		SizeType count = 0;
		srep_ = "";
		for (SizeType i = start; i < ntensors; ++i)
			count = data_[i]->relabelFrees(count);

		for (SizeType i = 0; i < start; ++i)
			count = data_[i]->relabelFrees(count);

		for (SizeType i = 0; i < ntensors; ++i)
			srep_ += data_[i]->sRep();
	}

	void shiftSummedBy(SizeType ms)
	{
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			data_[i]->shiftSummedBy(ms);
			srep_ += data_[i]->sRep();
		}
	}

	bool shouldAppear(const VectorSizeType& indices,
	                  SizeType times,
	                  PsimagLite::String c,
	                  bool verbose) const
	{
		SizeType total = indices.size();
		SizeType count = 0;
		bool flag = true;
		for (SizeType i = 0; i < total; ++i) {
			if (indices[i] == 0) continue;
			if (indices[i] == times) {
				count++;
				continue;
			}

			if (verbose) {
				PsimagLite::String str("Index " + c + ttos(i));
				str += " should appear " + ttos(times) + " times, not ";
				str += ttos(indices[i]) + " times\n";
				std::cerr<<str;
				flag = false;
			} else {
				return false;
			}
		}

		if (verbose)
			std::cerr<<"Found " + ttos(count) + " indices of type " + c + "\n";

		return flag;
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

	SizeType maxTag(char c) const
	{
		SizeType ntensors = data_.size();
		SizeType max = 0;
		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType tmp = data_[i]->maxTag(c);
			if (max < tmp) max = tmp;
		}

		return max;
	}

	PsimagLite::String srep_;
	VectorTensorStanza data_;
}; // class TensorSrep

} // namespace Mera
#endif // TENSORSREP_H
