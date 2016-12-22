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

public:

	typedef TensorStanza::VectorSizeType VectorSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	TensorBreakup(const TensorSrep& srep)
	    : srep_(srep), // deep copy of srep
	      tid_(computeInitialTid())
	{}

	void operator()(VectorStringType& vstr)
	{
		PsimagLite::String lhs;
		PsimagLite::String rhs;

		while (true) {
			VectorSizeType setS;
			TensorSrep::PairSizeType pair;
			bool found = findPairForBreakUp(pair,setS);
			if (!found) break;
			if (!breakUpTensor(lhs,rhs,pair,setS))
				break;
			vstr.push_back(lhs);
			vstr.push_back(rhs);
		}

		vstr.push_back("result");
		vstr.push_back(srep_.sRep());
	}

private:

	SizeType computeInitialTid() const
	{
		SizeType maxIdOfTtensors = 0;
		SizeType ntensors = srep_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			if (srep_(i).name() != "t") continue;
			SizeType tmp = srep_(i).id();
			if (tmp > maxIdOfTtensors) maxIdOfTtensors = tmp;
		}

		return maxIdOfTtensors;
	}

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

	void stringT0Part(PsimagLite::String& str,
	                  PsimagLite::String& actualStr,
	                  const VectorSizeType& setS,
	                  const TensorStanza& stanza) const
	{
		TensorStanza::IndexDirectionEnum in = TensorStanza::INDEX_DIR_IN;
		TensorStanza::IndexDirectionEnum out = TensorStanza::INDEX_DIR_OUT;
		str = stanza.name() + ttos(stanza.id());
		if (stanza.isConjugate()) str += "*";
		actualStr = str;
		bool seen = false;
		bool seenActual = false;
		SizeType addedDummyIndex = 0;
		for (SizeType i = 0; i < stanza.ins(); ++i) {
			TensorStanza::IndexTypeEnum legType = stanza.legType(i,in);
			SizeType legTag = stanza.legTag(i,in);
			bool isSummed = (legType == TensorStanza::INDEX_TYPE_SUMMED);
			bool notInSet = (isSummed && std::find(setS.begin(),setS.end(),legTag) == setS.end());
			if (!isSummed || notInSet) {
				str += (!seen) ? "(" : ",";
				actualStr += (!seenActual) ? "(" : ",";

				PsimagLite::String tmp = TensorStanza::indexTypeToString(legType) + ttos(legTag);
				str += tmp;
				actualStr += tmp;

				seen = seenActual = true;
			} else {
				if (!seenActual) actualStr += "(";
				else actualStr += ",";
				actualStr += "d" + ttos(addedDummyIndex++);
				seenActual = true;
			}
		}

		if (!seen) str += "(";
		if (!seen && !seenActual) actualStr += "(";
		seen = seenActual = false;
		for (SizeType i = 0; i < stanza.outs(); ++i) {
			TensorStanza::IndexTypeEnum legType = stanza.legType(i,out);
			SizeType legTag = stanza.legTag(i,out);
			bool isSummed = (legType == TensorStanza::INDEX_TYPE_SUMMED);
			bool notInSet = (isSummed && std::find(setS.begin(),setS.end(),legTag) == setS.end());
			if (!isSummed || notInSet) {
				str += (!seen) ? "|" : ",";
				actualStr += (!seenActual) ? "|" : ",";
				PsimagLite::String tmp = TensorStanza::indexTypeToString(legType) + ttos(legTag);
				str += tmp;
				actualStr += tmp;

				seen = seenActual = true;
			} else {
				if (!seenActual) actualStr += "|";
				else actualStr += ",";
				actualStr += "d" + ttos(addedDummyIndex++);
				seenActual = true;
			}
		}

		str += ")";
		actualStr += ")";
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

	bool breakUpTensor(PsimagLite::String& lhs,
	                   PsimagLite::String& rhs,
	                   const TensorSrep::PairSizeType& pair,
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
		PsimagLite::String str0("");
		PsimagLite::String actualStr0("");
		stringT0Part(str0,actualStr0,setS,stanza0);
		PsimagLite::String str1("");
		PsimagLite::String actualStr1("");
		stringT0Part(str1,actualStr1,setS,stanza1);
		str0 += str1;
		actualStr0 += actualStr1;

		PsimagLite::String t0 = buildTid(str0,tid_);
		TensorStanza t0stanzaActual(t0);
		VectorSizeType mapping;
		freeTheSummed(t0stanzaActual,mapping);
		TensorSrep rightHandSide(actualStr0);
		for (SizeType i =0; i < rightHandSide.size(); ++i)
			freeTheSummed(rightHandSide(i),mapping);

		rightHandSide.refresh();
		lhs = t0stanzaActual.sRep();
		rhs = rightHandSide.sRep();
		// build t1
		// mark tensor ind as erased
		srep_.setAsErased(ind);
		// tensor jnd becomes t0
		TensorStanza t0stanza(t0);
		srep_.replaceStanza(jnd,t0stanza);
		std::cerr<<"Remaining="<<srep_.sRep()<<"\n";
		tid_++;
		return true;
	}

	void freeTheSummed(TensorStanza& stanza,
	                   VectorSizeType& mapping) const
	{
		bool leftHandSide = (mapping.size() == 0);
		TensorStanza::IndexDirectionEnum in = TensorStanza::INDEX_DIR_IN;
		TensorStanza::IndexDirectionEnum out = TensorStanza::INDEX_DIR_OUT;
		SizeType counter = stanza.maxTag('f') + 1;
		if (leftHandSide)
			mapping.resize(stanza.maxTag('s')+1,1000);
		SizeType ins = stanza.ins();
		for (SizeType i = 0; i < ins; ++i) {
			if (!leftHandSide &&
			        stanza.legType(i,in) == TensorStanza::INDEX_TYPE_DUMMY) {
				stanza.legTypeChar(i,in) = 's';
				continue;
			}

			if (stanza.legType(i,in) != TensorStanza::INDEX_TYPE_SUMMED)
				continue;
			SizeType index = stanza.legTag(i,in);
			stanza.legTypeChar(i,in) = 'f';
			if (leftHandSide) mapping[index] = counter++;
			stanza.legTag(i,in) = mapping[index];
		}

		SizeType outs = stanza.outs();
		for (SizeType i = 0; i < outs; ++i) {
			if (!leftHandSide &&
			        stanza.legType(i,out) == TensorStanza::INDEX_TYPE_DUMMY) {
				stanza.legTypeChar(i,out) = 's';
				continue;
			}

			if (stanza.legType(i,out) != TensorStanza::INDEX_TYPE_SUMMED)
				continue;
			SizeType index = stanza.legTag(i,out);
			stanza.legTypeChar(i,out) = 'f';
			if (leftHandSide) mapping[index] = counter++;
			stanza.legTag(i,out) = mapping[index];
		}

		stanza.refresh();
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
