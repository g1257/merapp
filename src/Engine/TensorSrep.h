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
#include "Sort.h"
#include "IrreducibleIdentity.h"

namespace Mera {

class TensorSrep {

	typedef PsimagLite::Vector<TensorStanza*>::Type VectorTensorStanza;

public:

	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef PsimagLite::Vector<PairSizeType>::Type VectorPairSizeType;
	typedef TensorStanza::VectorSizeType VectorSizeType;
	typedef TensorStanza TensorStanzaType;

	explicit TensorSrep(PsimagLite::String srep)
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
	void eraseTensor(IrreducibleIdentity& irrIdentity,
	                 SizeType index,
	                 VectorSizeType* mapping)
	{
		assert(index < data_.size());

		bool identityIdIncreased = false;
		if (data_[index]->name() == "r")
			identityIdIncreased = addIrreducibleIdentity(irrIdentity);
		if (identityIdIncreased) irrIdentity.increase();

		VectorSizeType sErased;
		data_[index]->eraseTensor(sErased);
		SizeType ntensors = data_.size();
		SizeType count = maxTag('f') + 1;
		if (mapping) mapping->resize(maxTag('s') + 1,1000);

		srep_ = "";
		for (SizeType i = 0; i < ntensors; ++i) {
			count = data_[i]->uncontract(sErased,count,mapping);
			srep_ += data_[i]->sRep();
		}

		canonicalize();
	}

	void setAsErased(SizeType index)
	{
		assert(index < data_.size());
		data_[index]->setAsErased();
	}

	void replaceStanza(SizeType index, const TensorStanzaType& stanza)
	{
		assert(index < data_.size());
		*(data_[index]) = stanza;
		refresh();
	}

	const PsimagLite::String& sRep() const { return srep_; }

	char& legTypeChar(SizeType i,
	                  SizeType ind,
	                  TensorStanzaType::IndexDirectionEnum dir)
	{
		assert(i < data_.size());
		return data_[i]->legTypeChar(ind,dir);
	}

	SizeType& legTag(SizeType i,
	                 SizeType ind,
	                 TensorStanzaType::IndexDirectionEnum dir)
	{
		assert(i < data_.size());
		return data_[i]->legTag(ind,dir);
	}

	void refresh()
	{
		srep_ = "";
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			if (data_[i]->type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;
			data_[i]->refresh();
			srep_ += data_[i]->sRep();
		}
	}

	void conjugate()
	{
		srep_ = "";
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			data_[i]->conjugate();
			srep_ += data_[i]->sRep();
		}
	}

	void simplify(VectorPairSizeType& replacements)
	{
		simplifyFrees(replacements);
		simplifySummed();
	}

	void swapFree(SizeType ind, SizeType jnd)
	{
		TensorStanzaType::IndexDirectionEnum in = TensorStanzaType::INDEX_DIR_IN;
		TensorStanzaType::IndexDirectionEnum out = TensorStanzaType::INDEX_DIR_OUT;
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType legs = data_[i]->ins();
			for (SizeType j = 0; j < legs; ++j) {
				if (data_[i]->legType(j,in) != TensorStanzaType::INDEX_TYPE_FREE)
					continue;

				SizeType index = data_[i]->legTag(j,in);
				if (index == ind)
					data_[i]->legTag(j,in) = jnd;
				if (index == jnd)
					data_[i]->legTag(j,in) = ind;
			}

			legs = data_[i]->outs();
			for (SizeType j = 0; j < legs; ++j) {
				if (data_[i]->legType(j,out) != TensorStanzaType::INDEX_TYPE_FREE)
					continue;

				SizeType index = data_[i]->legTag(j,out);
				if (index == ind)
					data_[i]->legTag(j,out) = jnd;
				if (index == jnd)
					data_[i]->legTag(j,out) = ind;
			}
		}

		refresh();
	}

	SizeType size() const { return data_.size(); }

	const TensorStanza& operator()(SizeType ind) const
	{
		assert(ind < data_.size());
		return *data_[ind];
	}

	// FIXME: EXPOSES INTERNALS!!
	TensorStanza& operator()(SizeType ind)
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
			if (data_[i]->type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;
			SizeType ins = data_[i]->ins();
			SizeType outs = data_[i]->outs();
			if (ins + outs == 0) return false;
			for (SizeType j = 0; j < ins; ++j) {
				if (data_[i]->legType(j,in) == TensorStanzaType::INDEX_TYPE_SUMMED)
					summedIndices[data_[i]->legTag(j,in)]++;
				else if (data_[i]->legType(j,in) == TensorStanzaType::INDEX_TYPE_FREE)
					freeIndices[data_[i]->legTag(j,in)]++;
			}

			for (SizeType j = 0; j < outs; ++j) {
				if (data_[i]->legType(j,out) == TensorStanzaType::INDEX_TYPE_SUMMED)
					summedIndices[data_[i]->legTag(j,out)]++;
				else if (data_[i]->legType(j,out) == TensorStanzaType::INDEX_TYPE_FREE)
					freeIndices[data_[i]->legTag(j,out)]++;
			}
		}

		bool b1 = shouldAppear(summedIndices,2,"s",verbose);
		if (!verbose && !b1) return false;

		bool b2 = shouldAppear(freeIndices,1,"f",verbose);
		return (b1 & b2);
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

	SizeType findConjugate(SizeType ind) const
	{
		SizeType ntensors = data_.size();
		SizeType id = data_[ind]->id();
		PsimagLite::String name = data_[ind]->name();
		for (SizeType i = 0; i < ntensors; ++i) {
			if (data_[i]->type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;
			if (data_[i]->isConjugate() && data_[i]->name()==name && data_[i]->id() == id)
				return i;
		}

		return ntensors;
	}

	//	friend std::ostream& operator<<(std::ostream& os, const TensorSrep& ts)
	//	{
	//		os<<"tensorSrep.size="<<ts.size()<<"\n";
	//		for (SizeType i = 0; i < ts.size(); ++i) {
	//			os<<*(ts.data_[i])<<"\n";
	//		}

	//		return os;
	//	}

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

		SizeType mf = maxTag('f') + 1;
		shiftSummedBy(mf);

		srep_ = "";
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			data_[i]->contract(indicesToContract);
			srep_ += data_[i]->sRep();
		}

		verifySummed(0);

		fixDuplicatedFrees();

		canonicalize();
	}

	// we assume here that all stanzas are valid
	// but the srep_ might not be, due to duplicated frees
	void fixDuplicatedFrees()
	{
		SizeType mf = maxTag('f') + 1;
		VectorSizeType counter(mf,0);
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			VectorSizeType newfrees;
			data_[i]->loadSummedOrFree(newfrees,'f');
			if (newfrees.size() == 0) continue;
			// if newfrees in frees calculate replacements and replace
			VectorPairSizeType replacements;
			for (SizeType j = 0; j < newfrees.size(); ++j) {
				if (counter[newfrees[j]] == 1) {
					replacements.push_back(PairSizeType(newfrees[j],mf++));
				} else {
					assert(counter[newfrees[j]] == 0);
					counter[newfrees[j]]++;
				}
			}

			if (replacements.size() > 0) data_[i]->replaceSummedOrFrees(replacements, 'f');
		}

		srep_ = "";
		for (SizeType i = 0; i < ntensors; ++i) {
			if (data_[i]->type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;
			data_[i]->refresh();
			srep_ += data_[i]->sRep();
		}
	}

	bool verifySummed(VectorSizeType* usummed) const
	{
		VectorSizeType summed;
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i)
			data_[i]->loadSummedOrFree(summed,'s');

		SizeType n = summed.size();
		if (n&1)
			return false;

		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType iperm(n);
		sort.sort(summed,iperm);

		if (usummed) usummed->resize(1 + maxTag('s'),1000);
		for (SizeType i = 0; i < n; i += 2) {
			if (i + 1 >= summed.size()) break;
			if (summed[i] != summed[i + 1])
				return false;
			if (i > 0 && summed[i] == summed[i - 1])
				return false;
			if (!usummed) continue;
			assert(summed[i] < usummed->size());
			assert(summed[i + 1] < usummed->size());
			usummed->operator[](summed[i]) = usummed->operator[](summed[i + 1]) = i/2;
		}

		return true;
	}

	void canonicalize()
	{
		simplify();

		TensorStanzaType::IndexDirectionEnum in = TensorStanzaType::INDEX_DIR_IN;
		TensorStanzaType::IndexDirectionEnum out = TensorStanzaType::INDEX_DIR_OUT;

		VectorSizeType frees(maxTag('f') + 1,0);
		SizeType ntensors = data_.size();
		SizeType counter = 0;
		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType legs = data_[i]->ins();
			for (SizeType j = 0; j < legs; ++j) {
				if (data_[i]->legType(j,in) != TensorStanzaType::INDEX_TYPE_FREE)
					continue;
				frees[data_[i]->legTag(j,in)] = counter++;
			}

			legs = data_[i]->outs();
			for (SizeType j = 0; j < legs; ++j) {
				if (data_[i]->legType(j,out) != TensorStanzaType::INDEX_TYPE_FREE)
					continue;
				frees[data_[i]->legTag(j,out)] = counter++;
			}
		}

		srep_ = "";
		for (SizeType i = 0; i < ntensors; ++i) {
			if (data_[i]->type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;
			data_[i]->setIndices(frees,'f');
			data_[i]->refresh();
			srep_ += data_[i]->sRep();
		}

		simplifySummed();
	}

	void simplify()
	{
		while (simplifyOnce());
	}

	bool simplifyOnce()
	{
		bool simplificationHappended = false;
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			if (data_[i]->isConjugate()) continue;
			SizeType j = findConjugate(i);
			if (j >= data_.size()) continue; // no conjugate
			if (!inputsMatch(i,j)) continue;
			if (!simplify(i,j)) continue;
			simplificationHappended = true;
			break; // only one simplification
		}

		return simplificationHappended;
	}

	bool simplify(SizeType ind, SizeType jnd)
	{
		bool simplificationHappended = false;
		SizeType outs = data_[ind]->outs();
		VectorPairSizeType replacements(outs);
		if (!computeReplacements(replacements,ind,jnd))
			return simplificationHappended;

		data_[ind]->setAsErased();
		data_[jnd]->setAsErased();

		SizeType ntensors = data_.size();
		srep_ = "";
		for (SizeType i = 0; i < ntensors; ++i) {
			if (data_[i]->replaceSummedOrFrees(replacements,'s'))
				simplificationHappended = true;
			srep_ += data_[i]->sRep();
		}

		return simplificationHappended;
	}

	void simplifyFrees(VectorPairSizeType& replacements)
	{
		if (replacements.size() == 0) return;
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i)
			data_[i]->replaceSummedOrFrees(replacements,'f');

		refresh();
	}

	void simplifySummed()
	{
		VectorSizeType usummed;
		if (!verifySummed(&usummed))
			throw PsimagLite::RuntimeError("simplifySummed: Invalid Srep\n");

		srep_ = "";
		SizeType ntensors = data_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			data_[i]->setIndices(usummed,'s');
			srep_ += data_[i]->sRep();
		}
	}

	bool computeReplacements(VectorPairSizeType& replacements,
	                         SizeType ind,
	                         SizeType jnd) const
	{
		SizeType outs = data_[ind]->outs();
		if (outs != data_[jnd]->outs()) return false;

		for (SizeType i = 0; i < outs; ++i) {
			if (data_[ind]->legType(i,TensorStanzaType::INDEX_DIR_OUT) !=
			        TensorStanzaType::INDEX_TYPE_SUMMED) return false;
			if (data_[jnd]->legType(i,TensorStanzaType::INDEX_DIR_OUT) !=
			        TensorStanzaType::INDEX_TYPE_SUMMED) return false;
			SizeType s1 = data_[ind]->legTag(i,TensorStanzaType::INDEX_DIR_OUT);
			SizeType s2 = data_[jnd]->legTag(i,TensorStanzaType::INDEX_DIR_OUT);
			replacements[i] = PairSizeType((s1 < s2) ? s2 : s1,(s1 < s2) ? s1 : s2);
		}

		return true;
	}

	bool inputsMatch(SizeType ind, SizeType jnd) const
	{
		SizeType ins = data_[ind]->ins();
		if (ins != data_[jnd]->ins()) return false;
		for (SizeType i = 0; i < ins; ++i) {
			if (data_[ind]->legType(i,TensorStanzaType::INDEX_DIR_IN) !=
			        TensorStanzaType::INDEX_TYPE_SUMMED) return false;
			if (data_[jnd]->legType(i,TensorStanzaType::INDEX_DIR_IN) !=
			        TensorStanzaType::INDEX_TYPE_SUMMED) return false;
			if (data_[ind]->legTag(i,TensorStanzaType::INDEX_DIR_IN) !=
			        data_[jnd]->legTag(i,TensorStanzaType::INDEX_DIR_IN))
				return false;
		}

		return true;
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

	bool addIrreducibleIdentity(IrreducibleIdentity& irrIdentity)
	{
		SizeType loc1 = 0;
		VectorSizeType summed1;
		PsimagLite::String str("");
		if ((str = findRandRifContracted(loc1,summed1,irrIdentity.maxIndex())) == "")
			return false;

		data_[loc1]->setIndices(summed1,'s');
		data_.push_back(new TensorStanzaType(str));

		SizeType ntensors = data_.size();
		srep_ = "";
		for (SizeType i = 0; i < ntensors; ++i)
			srep_ += data_[i]->sRep();
		return true;
	}

	PsimagLite::String findRandRifContracted(SizeType& loc1,
	                                         VectorSizeType& usummed1,
	                                         SizeType identityId) const
	{
		SizeType ntensors = data_.size();
		SizeType flag0 = 0;
		SizeType flag1 = 0;
		SizeType loc0 = 0;
		for (SizeType i = 0; i < ntensors; ++i) {
			if (data_[i]->name() != "r") continue;
			bool b = data_[i]->isConjugate();
			if (b) {
				flag0++;
				loc0 = i;
				continue;
			}

			flag1++;
			loc1 = i;
		}

		assert(flag0 < 2 && flag1 < 2);
		// r and r* both need to be present
		if (flag0 != 1 || flag1 != 1) return "";

		VectorSizeType summed0;
		assert(loc0 < data_.size());
		data_[loc0]->loadSummedOrFree(summed0,'s');

		VectorSizeType summed1;
		assert(loc1 < data_.size());
		data_[loc1]->loadSummedOrFree(summed1,'s');

		// check that r or r* haven't been deleted
		if (summed1.size() == 0 || summed0.size() == 0) return "";

		SizeType max = *std::max_element(summed1.begin(),summed1.end());
		usummed1.resize(max + 1,0);
		for (SizeType i = 0; i < usummed1.size(); ++i)
			usummed1[i] = i;

		// one index summed must be common to r and r*
		flag0 = 0;
		PsimagLite::String str = "i" + ttos(identityId) + "(s";
		for (SizeType i = 0; i < summed0.size(); ++i) {
			for (SizeType j = 0; j < summed1.size(); ++j) {
				if (summed0[i] == summed1[j]) {
					flag0 = 1;
					assert(summed0[i] < usummed1.size());
					usummed1[summed0[i]] = maxTag('s') + 1;
					str += ttos(summed0[i]) + "|s" + ttos(usummed1[summed0[i]]);
					break;
				}
			}
		}

		return (flag0 == 0) ? "" : str + ")";
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

		for (SizeType i = 0; i < ntensors; ++i) {
			if (data_[i]->type() == TensorStanzaType::TENSOR_TYPE_ERASED)
				continue;
			data_[i]->refresh();
			srep_ += data_[i]->sRep();
		}
	}

	void shiftSummedBy(SizeType ms)
	{
		SizeType ntensors = data_.size();
		srep_ = "";
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

	PsimagLite::String srep_;
	VectorTensorStanza data_;
}; // class TensorSrep

} // namespace Mera
#endif // TENSORSREP_H
