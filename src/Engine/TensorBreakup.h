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
#ifndef TENSORBREAKUP_H
#define TENSORBREAKUP_H
#include "TensorSrep.h"
#include "Vector.h"

namespace Mera {

class TensorBreakup {

	typedef TensorStanza::VectorSizeType VectorSizeType;

public:

	TensorBreakup(PsimagLite::String str) : srep_(str),tid_(0)
	{}

	void operator()()
	{
		while (true) {
			VectorSizeType setS;
			TensorSrep::PairSizeType pair;
			bool found = findPairForBreakUp(pair,setS);
			if (!found) break;
			if (!breakUpTensor(pair,setS))
				break;
		}
	}

private:

	void getAllStags(VectorSizeType& v0,
	                 const TensorStanza& stanza0) const
	{
		TensorStanza::IndexDirectionEnum in = TensorStanza::INDEX_DIR_IN;
		TensorStanza::IndexDirectionEnum out = TensorStanza::INDEX_DIR_OUT;
		for (SizeType i = 0; i < stanza0.ins(); ++i) {
			if (stanza0.legType(i,in) != TensorStanza::INDEX_TYPE_SUMMED)
				continue;
			v0.push_back(stanza0.legTag(i,in));
		}

		for (SizeType i = 0; i < stanza0.outs(); ++i) {
			if (stanza0.legType(i,out) != TensorStanza::INDEX_TYPE_SUMMED)
				continue;
			v0.push_back(stanza0.legTag(i,out));
		}
	}

	PsimagLite::String stringT0Part(const VectorSizeType& setS,
	                                const TensorStanza& stanza) const
	{
		TensorStanza::IndexDirectionEnum in = TensorStanza::INDEX_DIR_IN;
		TensorStanza::IndexDirectionEnum out = TensorStanza::INDEX_DIR_OUT;
		PsimagLite::String str = stanza.name() + ttos(stanza.id());
		SizeType counter = 0;
		for (SizeType i = 0; i < stanza.ins(); ++i) {
			TensorStanza::IndexTypeEnum legType = stanza.legType(i,in);
			SizeType legTag = stanza.legTag(i,in);
			bool isSummed = (legType == TensorStanza::INDEX_TYPE_SUMMED);
			bool notInSet = (isSummed && std::find(setS.begin(),setS.end(),legTag) == setS.end());
			if (!isSummed || notInSet) {
				if (counter++ == 0)
					str += "(";
				else
					str += ",";
				str += TensorStanza::indexTypeToString(legType) + ttos(legTag);
			}
		}

		if (counter == 0) str += "(";
		counter = 0;
		for (SizeType i = 0; i < stanza.outs(); ++i) {
			TensorStanza::IndexTypeEnum legType = stanza.legType(i,out);
			SizeType legTag = stanza.legTag(i,out);
			bool isSummed = (legType == TensorStanza::INDEX_TYPE_SUMMED);
			bool notInSet = (isSummed && std::find(setS.begin(),setS.end(),legTag) == setS.end());
			if (!isSummed || notInSet) {
				if (counter++ == 0)
					str += "|";
				else
					str += ",";
				str += TensorStanza::indexTypeToString(legType) + ttos(legTag);
			}
		}

		return str + ")";
	}

	PsimagLite::String buildT0Part(const TensorStanza& stanza,
	                               SizeType legs,
	                               TensorStanza::IndexDirectionEnum inOrOut) const
	{
		PsimagLite::String str("");
		SizeType counter = 0;
		for (SizeType i = 0; i < legs; ++i) {
			TensorStanza::IndexTypeEnum legType = stanza.legType(i,inOrOut);
			SizeType legTag = stanza.legTag(i,inOrOut);
			if (counter++ > 0) str += ",";
			str += TensorStanza::indexTypeToString(legType) + ttos(legTag);
		}

		return str;
	}

	PsimagLite::String buildTid(PsimagLite::String str,
	                           SizeType id) const
	{
		TensorStanza::IndexDirectionEnum in = TensorStanza::INDEX_DIR_IN;
		TensorStanza::IndexDirectionEnum out = TensorStanza::INDEX_DIR_OUT;
		TensorSrep srep(str);
		SizeType ntensors = srep.size();
		PsimagLite::String t0In("t" + ttos(id) + "(");
		SizeType counter = 0;
		for (SizeType i = 0; i < ntensors; ++i) {
			const TensorStanza& stanza = srep(i);
			SizeType ins = stanza.ins();
			SizeType outs = stanza.outs();
			PsimagLite::String tmp = (stanza.isConjugate()) ?
			            buildT0Part(stanza,outs,out) : buildT0Part(stanza,ins,in);
			if (tmp == "") continue;
			if (counter++ > 0) t0In += ",";
			t0In += tmp;
		}

		counter = 0;
		PsimagLite::String t0Out("");
		for (SizeType i = 0; i < ntensors; ++i) {
			const TensorStanza& stanza = srep(i);
			SizeType ins = stanza.ins();
			SizeType outs = stanza.outs();
			PsimagLite::String tmp = (stanza.isConjugate()) ?
			            buildT0Part(stanza,ins,in) : buildT0Part(stanza,outs,out);
			if (tmp == "") continue;
			if (counter++ > 0) t0Out += ",";
			t0Out += tmp;
		}

		if (t0Out != "") t0In += "|";
		return t0In + t0Out + ")";
	}

	bool breakUpTensor(const TensorSrep::PairSizeType& pair,
	                   const VectorSizeType& setS)
	{
		SizeType ind = pair.first;
		SizeType jnd = pair.second;
		if (ind >= jnd) return false;
		if (srep_.size() == 2) return false;

		TensorStanza stanza0 = srep_(ind);
		TensorStanza stanza1 = srep_(jnd);
		std::cerr<<"We're going to break tensor ";
		PsimagLite::String c = (stanza0.isConjugate()) ? "*" : "";
		std::cerr<<stanza0.name()<<stanza0.id()<<c<<" and ";
		c = (stanza1.isConjugate()) ? "*" : "";
		std::cerr<<stanza1.name()<<stanza1.id()<<c;
		std::cerr<<" from "<<srep_.sRep()<<"\n";

		// build t0
		PsimagLite::String str0 = stringT0Part(setS,stanza0);
		PsimagLite::String str1 = stringT0Part(setS,stanza1);
		str0 += str1;
		PsimagLite::String t0 = buildTid(str0,tid_);
		std::cout<<t0<<"="<<str0<<"\n";

		// build t1
		// mark tensor ind as erased
		srep_.setAsErased(ind);
		// tensor jnd becomes t0
		TensorStanza t0stanza(t0);
		srep_.replaceStanza(jnd,t0stanza);
		std::cout<<"t1="<<srep_.sRep()<<"\n";
		tid_++;
		return true;
	}

	void findCommonSetS(VectorSizeType& setS,
	                    const TensorStanza& stanza0,
	                    const TensorStanza& stanza1) const
	{
		// find commont setS
		VectorSizeType v0;
		getAllStags(v0,stanza0);
		VectorSizeType v1;
		getAllStags(v1,stanza1);

		setS.clear();
		for (SizeType i = 0; i < v0.size(); ++i)
			for (SizeType j = 0; j < v1.size(); ++j)
				if (v0[i] == v1[j]) setS.push_back(v0[i]);
	}

	bool findPairForBreakUp(TensorSrep::PairSizeType& pair,
	                        VectorSizeType& setS) const
	{
		if (srep_.size() < 2) return false;
		for (SizeType start = 0; start < srep_.size()-1; ++start) {
			setS.clear();
			findPairForBreakUp(pair,setS,start);
			if (pair.first < pair.second && setS.size() > 0)
				return true;
		}

		return false;
	}

	void findPairForBreakUp(TensorSrep::PairSizeType& pair,
	                        VectorSizeType& setS,
	                        SizeType start) const
	{
		pair.first = pair.second = 0;
		SizeType ntensors = srep_.size();
		SizeType counter = 0;
		for (SizeType i = start; i < ntensors; ++i) {
			const TensorStanza& stanza = srep_(i);
			if (stanza.type() == TensorStanza::TENSOR_TYPE_ERASED)
				continue;
			if (stanza.name() == "t") continue;
			if (counter == 0) {
				pair.first = i;
				counter++;
			} else if (counter == 1) {
					findCommonSetS(setS,srep_(pair.first),stanza);
					if (setS.size() > 0) {
						pair.second = i;
						break;
					}
			}
		}
	}

	TensorSrep srep_;
	SizeType tid_;
};

}
#endif // TENSORBREAKUP_H
