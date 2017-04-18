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
#ifndef MERA_ENVIRON_H
#define MERA_ENVIRON_H
#include <iostream>
#include "ParametersForMera.h"
#include "TensorSrep.h"
#include "MeraBuilder.h"

namespace Mera {

template<typename ComplexOrRealType, typename SymmetryLocalType>
class MeraEnviron {

	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef ParametersForMera<ComplexOrRealType> ParametersForMeraType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::Vector<TensorSrep*>::Type VectorTensorSrepType;
	typedef MeraBuilder<ComplexOrRealType> MeraBuilderType;

public:

	MeraEnviron(const MeraBuilderType& builder,
	            const ParametersForMeraType& params,
	            PsimagLite::String dimensionSrep,
	            SymmetryLocalType& symmLocal)
	    : builder_(builder),
	      params_(params),
	      dimensionSrep_(dimensionSrep),
	      symmLocal_(symmLocal),
	      tensorSrep_(builder()),
	      envs_(""),
	      dsrep_("")
	{
		sizeOfRoot_ = findSizeOfRoot();
		SizeType counterForOutput = 100;
		VectorSizeType limits;
		findLimits(limits);
		for (SizeType i = 0; i < tensorSrep_.size(); ++i) {
			counterForOutput += environForTensor(i,
			                                     counterForOutput,
			                                     limits);
		}

		energies();
	}

	const PsimagLite::String& environs() const
	{
		return envs_;
	}

	const PsimagLite::String& dimensionSrep() const
	{
		return dsrep_;
	}

private:

	// find Y (environment) for this tensor
	SizeType environForTensor(SizeType ind,
	                          SizeType counterForOutput,
	                          const VectorSizeType& limits)
	{
		SizeType id = tensorSrep_(ind).id();
		PsimagLite::String name = tensorSrep_(ind).name();
		SizeType connections = params_.hamiltonianConnection.size();
		VectorStringType vstr(connections,"");
		VectorStringType argForOutput(connections,"");
		VectorStringType vdsrep(connections,"");
		SizeType terms = 0;

		for (SizeType c = 0; c < connections; ++c) {
			if (params_.hamiltonianConnection[c] == 0.0) continue;
			TensorSrep tmp = environForTensorOneSite(ind, c);
			vstr[c] = tmp.sRep();
			argForOutput[c] = calcArgForOutput(vdsrep[c],tmp);
			if (vstr[c] != "") ++terms;
		}

		if (terms == 0) return terms;

		PsimagLite::String thisEnv("TensorId=" + name + "," + ttos(id) + "\n");
		thisEnv += "Terms=" + ttos(terms) + "\n";
		thisEnv += "IgnoreTerm=" + ttos(2*connections+1) + "\n";
		PairSizeType layer = findLayerNumber(name, id, limits);
		thisEnv += "Layer=" + ttos(layer.first) + "\n";
		thisEnv += "FirstOfLayer=" + ttos(layer.second) + "\n";
		bool isRootTensor = (tensorSrep_(ind).name() == "r");
		for (SizeType c = 0; c < connections; ++c) {
			if (vstr[c] == "") continue;
			PsimagLite::String tmp = "u" + ttos(counterForOutput++);
			thisEnv += "Environ=" + tmp + argForOutput[c] + "=" + vstr[c] + "\n";
			dsrep_ += tmp + vdsrep[c];
			if (isRootTensor)
				irreducibleIdentityDsrep(tmp + argForOutput[c], vstr[c]);
		}

		thisEnv += "\n";
		envs_ += thisEnv;

		return terms;
	}

	TensorSrep environForTensorOneSite(SizeType ind,
	                                   SizeType site) const
	{
		const TensorSrep& energySrep = builder_.energy(site);
		TensorSrep tensorSrep4(energySrep);
		tensorSrep4.eraseTensor(irreducibleIdentity_,ind,0);
		bool verbose = false;
		if (!tensorSrep4.isValid(verbose))
			throw PsimagLite::RuntimeError("Invalid tensor\n");
		SizeType jnd = tensorSrep4.findConjugate(ind);
		bool hasConjugate = (jnd < tensorSrep4.size());
		bool isRootTensor = (tensorSrep_(ind).name() == "r");
		if (!hasConjugate) {
			if (isRootTensor) {
				throw PsimagLite::RuntimeError("Environ for root: INTERNAL ERROR\n");
			} else {
				std::cerr<<"EMPTY_ENVIRON="<<tensorSrep4.sRep()<<"\n";
				return TensorSrep("");
			}
		} else if (isRootTensor) {
			tensorSrep4.eraseTensor(irreducibleIdentity_,jnd,0);
			// use energySrep to compute r size
			//
			if (!tensorSrep4.isValid(verbose))
				throw PsimagLite::RuntimeError("Invalid tensor\n");
		}

		return tensorSrep4;
	}

	PsimagLite::String calcArgForOutput(PsimagLite::String& dsrep,
	                                    const TensorSrep& srep) const
	{
		VectorStringType ins;
		VectorStringType outs;
		SizeType ntensors = srep.size();
		for (SizeType i = 0; i < ntensors; ++i)
			calcArgForOutput(ins,outs,srep(i));

		PsimagLite::String ret = "(";
		dsrep = "(";
		for (SizeType i = 0; i < ins.size(); ++i) {
			ret += ins[i];
			dsrep += "D1";
			if (i + 1 >= ins.size()) continue;
			ret += ",";
			dsrep += ",";
		}

		if (outs.size() > 0) {
			ret += "|";
			dsrep += "|";
		}

		for (SizeType i = 0; i < outs.size(); ++i) {
			ret += outs[i];
			dsrep += "D1";
			if (i + 1 >= outs.size()) continue;
			ret += ",";
			dsrep += ",";
		}

		dsrep += ")";
		return ret + ")";
	}

	void calcArgForOutput(VectorStringType& ins,
	                      VectorStringType& outs,
	                      const TensorStanza& stanza) const
	{
		SizeType insNumber = stanza.ins();
		for (SizeType i = 0; i < insNumber; ++i) {
			if (stanza.legType(i) != TensorStanza::INDEX_TYPE_FREE) continue;
			PsimagLite::String tmp = "f" + ttos(stanza.legTag(i));
			if (stanza.isConjugate())
				outs.push_back(tmp);
			else
				ins.push_back(tmp);
		}

		for (SizeType i = 0; i < stanza.outs(); ++i) {
			if (stanza.legType(i + insNumber) != TensorStanza::INDEX_TYPE_FREE) continue;
			PsimagLite::String tmp = "f" + ttos(stanza.legTag(i + insNumber));
			if (stanza.isConjugate())
				ins.push_back(tmp);
			else
				outs.push_back(tmp);
		}
	}

	void irreducibleIdentityDsrep(PsimagLite::String left, PsimagLite::String right)
	{
		TensorSrep rightSrep(right);
		SizeType ntensors = rightSrep.size();
		int indexOfIdentity = -1;
		for (SizeType i = 0; i < ntensors; ++i) {
			if (rightSrep(i).name() == "i") {
				indexOfIdentity = i;
				break;
			}
		}

		if (indexOfIdentity < 0) return;

		TensorStanza leftSrep(left);
		SizeType ins = rightSrep(indexOfIdentity).ins();
		VectorSizeType frees;
		for (SizeType i = 0; i < ins; ++i) {
			if (rightSrep(indexOfIdentity).legType(i) != TensorStanza::INDEX_TYPE_FREE)
				continue;
			frees.push_back(rightSrep(indexOfIdentity).legTag(i));
		}

		VectorSizeType frees2;
		ins = leftSrep.ins();
		for (SizeType i = 0; i < ins; ++i) {
			if (leftSrep.legType(i) != TensorStanza::INDEX_TYPE_FREE)
				continue;
			SizeType tag = leftSrep.legTag(i);
			if (std::find(frees.begin(), frees.end(), tag) != frees.end())
				continue;
			frees2.push_back(tag);
		}

		// find dimensions of frees2
		SizeType sizeWithoutIrrIdentity = 1;
		for (SizeType i = 0; i < frees2.size(); ++i) {
			sizeWithoutIrrIdentity *= findDimension(rightSrep, frees2[i]);
		}

		if (sizeOfRoot_ % sizeWithoutIrrIdentity != 0)
			throw PsimagLite::RuntimeError("irreducibleIdentityDsrep\n");

		SizeType tmp = sizeOfRoot_/sizeWithoutIrrIdentity;
		SizeType id = rightSrep(indexOfIdentity).id();
		dsrep_ += "i" + ttos(id) + "(D" + ttos(tmp) + "|D" + ttos(tmp) + ")";
		symmLocal_.addIdentity(id, tmp);

	}

	SizeType findDimension(const TensorSrep& srep, SizeType indexOfFree) const
	{
		SizeType ntensors = srep.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType legs = srep(i).legs();
			for (SizeType j = 0; j < legs; ++j) {
				if (srep(i).legType(j) != TensorStanza::INDEX_TYPE_FREE)
					continue;
				if (srep(i).legTag(j) != indexOfFree)
					continue;
				return findDimension(srep(i).name(), srep(i).id(), j);
			}
		}

		throw PsimagLite::RuntimeError("findDimension(1)\n");
	}

	SizeType findDimension(PsimagLite::String name,
	                       SizeType id,
	                       SizeType legIndex) const
	{
		SizeType ntensors = dimensionSrep_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			if (dimensionSrep_(i).name() != name) continue;
			if (dimensionSrep_(i).id() != id) continue;

			return dimensionSrep_(i).legTag(legIndex);
		}

		throw PsimagLite::RuntimeError("findDimension(2)\n");
	}

	SizeType findSizeOfRoot() const
	{
		SizeType ntensors = dimensionSrep_.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			if (dimensionSrep_(i).name() != "r") continue;

			SizeType ins = dimensionSrep_(i).ins();
			SizeType ret = 1;
			for (SizeType j = 0; j < ins; ++j)
				ret *= dimensionSrep_(i).legTag(j);
			return ret;
		}

		throw PsimagLite::RuntimeError("findSizeOfRoot\n");
	}

	void energies()
	{
		PsimagLite::String e("TensorId=E,0\n");
		PsimagLite::String d("");
		SizeType terms = params_.hamiltonianConnection.size();
		SizeType effectiveTerms = 0;
		for (SizeType i = 0; i < terms; ++i) {
			if (params_.hamiltonianConnection[i] == 0.0) continue;
			++effectiveTerms;
		}

		e += "Terms=" + ttos(effectiveTerms) + "\n";
		e += "IgnoreTerm=" + ttos(terms+1) + "\n";
		for (SizeType i = 0; i < terms; ++i) {
			if (params_.hamiltonianConnection[i] == 0.0) continue;
			PsimagLite::String rhs = builder_.energy(i).sRep();
			assert(rhs != "");
			e += "Environ=e" + ttos(i) + "()=" + rhs + "\n";
			d += "e" + ttos(i) + "()";
		}

		envs_ += e;
		dsrep_ += d;
	}

	PairSizeType findLayerNumber(PsimagLite::String name,
	                         SizeType id,
	                         const VectorSizeType& limits) const
	{
		if (name != "w" && name != "u")
			return PairSizeType(0,1); // FIXME: think about layer number for other tensors

		SizeType n = limits.size();
		for (SizeType i = 0; i < n; ++i) {
			SizeType firstOfLayer = (i > 0) ? limits[i-1] : 0;
			if (id < limits[i]) return PairSizeType(i,firstOfLayer);
		}

		throw PsimagLite::RuntimeError("findLayerNumber: not found\n");
	}

	void findLimits(VectorSizeType& limits) const
	{
		SizeType arity = builder_.arity();
		SizeType sites = builder_.sites();

		if (arity != 2) {
			PsimagLite::RuntimeError("findLimits: supported only for arity=2\n");
		}

		limits.push_back(sites/arity);
		SizeType tmp = sites/arity;
		SizeType i = 1;

		while (tmp > 1) {
			limits.push_back(limits[i-1] + tmp/arity);
			tmp /= arity;
			++i;
		}
	}

	MeraEnviron(const MeraEnviron&);

	MeraEnviron& operator=(const MeraEnviron&);

	const MeraBuilderType& builder_;
	const ParametersForMeraType& params_;
	TensorSrep dimensionSrep_;
	SymmetryLocalType& symmLocal_;
	SizeType sizeOfRoot_;
	TensorSrep tensorSrep_;
	PsimagLite::String envs_;
	PsimagLite::String dsrep_;
	mutable IrreducibleIdentity irreducibleIdentity_;
}; //class

} //namespace

#endif
