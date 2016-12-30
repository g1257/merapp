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
#include "ParametersForSolver.h"
#include "TensorSrep.h"

namespace Mera {

class MeraEnviron {

	typedef ParametersForSolver ParametersForSolverType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::Vector<TensorSrep*>::Type VectorTensorSrepType;
public:

	MeraEnviron(const TensorSrep& srep, const ParametersForSolver& params)
	    : params_(params), tensorSrep_(srep), envs_(""), dsrep_("")
	{
		SizeType sites = tensorSrep_.maxTag('f');
		assert(params_.hamiltonianTerm.size() == sites);
		VectorTensorSrepType energy(sites, 0);
		for (SizeType site = 0; site < sites; ++site) {
			if (!params_.hamiltonianTerm[site]) continue;
			energy[site] = buildEnergyTerm(site);
		}

		SizeType counterForOutput = 100;
		for (SizeType i = 0; i < tensorSrep_.size(); ++i) {
			counterForOutput += environForTensor(i, counterForOutput, energy);
		}

		for (SizeType site = 0; site < sites; ++site) {
			if (!params_.hamiltonianTerm[site]) continue;
			delete energy[site];
			energy[site] = 0;
		}
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

	MeraEnviron(const MeraEnviron&);

	MeraEnviron& operator=(const MeraEnviron&);

	// find Y (environment) for this tensor
	SizeType environForTensor(SizeType ind,
	                          SizeType counterForOutput,
	                          const VectorTensorSrepType& energy)
	{
		SizeType id = tensorSrep_(ind).id();
		PsimagLite::String name = tensorSrep_(ind).name();
		SizeType sites = tensorSrep_.maxTag('f');
		VectorStringType vstr(sites,"");
		VectorStringType argForOutput(sites,"");
		VectorStringType vdsrep(sites,"");
		SizeType terms = 0;
		assert(params_.hamiltonianTerm.size() == sites);
		for (SizeType site = 0; site < sites; ++site) {
			if (!params_.hamiltonianTerm[site]) continue;
			assert(energy[site]);
			TensorSrep tmp = environForTensorOneSite(ind, site, *(energy[site]));
			vstr[site] = tmp.sRep();
			argForOutput[site] = calcArgForOutput(vdsrep[site],tmp);
			if (vstr[site] != "") ++terms;
		}

		PsimagLite::String thisEnv("TensorId=" + name + "," + ttos(id) + "\n");
		thisEnv += "Terms=" + ttos(terms) + "\n";
		thisEnv += "IgnoreTerm=" + ttos(2*sites+1) + "\n";
		thisEnv += "Layer=0\n"; // FIXME
		for (SizeType site = 0; site < sites; ++site) {
			if (vstr[site] == "") continue;
			PsimagLite::String tmp = "u" + ttos(counterForOutput++);
			thisEnv += "Environ=" + tmp + argForOutput[site] + "=" + vstr[site] + "\n";
			dsrep_ += tmp + vdsrep[site];
		}

		thisEnv += "\n";
		envs_ += thisEnv;
		return terms;
	}

	TensorSrep* buildEnergyTerm(SizeType site) const
	{
		TensorSrep tensorSrep2(tensorSrep_);
		tensorSrep2.conjugate();
		tensorSrep2.swapFree(0,site);
		tensorSrep2.swapFree(1,site+1);
		PsimagLite::String str3("h0(f");
		str3 += ttos(site+2) + ",f";
		str3 += ttos(site+3) + "|f";
		str3 += ttos(site) + ",f";
		str3 += ttos(site+1) + ")\n";
		TensorSrep tensorSrep3(str3);
		TensorSrep::VectorSizeType indicesToContract(2,site);
		indicesToContract[1] = site + 1;
		TensorSrep* tensorSrep4 = new TensorSrep(tensorSrep_);
		tensorSrep4->contract(tensorSrep3,indicesToContract);
		if (!tensorSrep4->isValid(true))
			throw PsimagLite::RuntimeError("Invalid tensor\n");
		std::cerr<<"LOWER="<<tensorSrep2.sRep()<<"\n";
		std::cerr<<"UPPER="<<tensorSrep4->sRep()<<"\n";
		tensorSrep4->contract(tensorSrep2);
		std::cerr<<"ENERGY="<<tensorSrep4->sRep()<<"\n";
		if (!tensorSrep4->isValid(true))
			throw PsimagLite::RuntimeError("Invalid tensor\n");
		return tensorSrep4;
	}

	TensorSrep environForTensorOneSite(SizeType ind,
	                                   SizeType,
	                                   const TensorSrep& energySrep) const
	{
		TensorSrep tensorSrep4(energySrep);
		tensorSrep4.eraseTensor(ind);
		if (!tensorSrep4.isValid(true))
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
			tensorSrep4.eraseTensor(jnd);
			if (!tensorSrep4.isValid(true))
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
			dsrep += "1";
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
			dsrep += "1";
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
		TensorStanza::IndexDirectionEnum in = TensorStanza::INDEX_DIR_IN;
		TensorStanza::IndexDirectionEnum out = TensorStanza::INDEX_DIR_OUT;

		for (SizeType i = 0; i < stanza.ins(); ++i) {
			if (stanza.legType(i,in) != TensorStanza::INDEX_TYPE_FREE) continue;
			PsimagLite::String tmp = "f" + ttos(stanza.legTag(i,in));
			if (stanza.isConjugate())
				outs.push_back(tmp);
			else
				ins.push_back(tmp);
		}

		for (SizeType i = 0; i < stanza.outs(); ++i) {
			if (stanza.legType(i,out) != TensorStanza::INDEX_TYPE_FREE) continue;
			PsimagLite::String tmp = "f" + ttos(stanza.legTag(i,out));
			if (stanza.isConjugate())
				ins.push_back(tmp);
			else
				outs.push_back(tmp);
		}
	}

	const ParametersForSolver& params_;
	const TensorSrep& tensorSrep_;
	PsimagLite::String envs_;
	PsimagLite::String dsrep_;
}; //class

} //namespace

#endif
