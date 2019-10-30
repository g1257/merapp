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
#ifndef MERATOTIKZ_H
#define MERATOTIKZ_H
#include "Vector.h"
#include "TypeToString.h"
#include <cassert>
#include <iostream>
#include "TensorSrep.h"
#include "ProgramGlobals.h"

namespace Mera {

template<typename RealType>
class MeraToTikz {

	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<PairSizeType>::Type VectorPairSizeType;

public:

	typedef std::pair<PsimagLite::String, SizeType> PairStringSizeType;

	MeraToTikz(PsimagLite::String srep, SizeType sites)
	    : srep_(srep),tauMax_(ProgramGlobals::logBase2Strict(sites)-1)
	{
		buildPacking(sites);
		fillBuffer();
	}

	friend std::ostream& operator<<(std::ostream& os, const MeraToTikz& mt)
	{
		printHeader(os);
		os<<buffer_;
		printFooter(os);
		return os;
	}

private:

	void fillBuffer()
	{
		TensorSrep tensorSrep(srep_);

		buffer_ = "%" + ProgramGlobals::addLf(srep_,80);
		buffer_ += "\n";

		RealType dx = 1.;
		RealType dy0 = 1.;

		SizeType ntensors = tensorSrep.size();
		VectorRealType x(ntensors);
		VectorRealType y(ntensors);
		computeCoordinates(x,y,dx,tensorSrep);
		SizeType layer = 0;
		PsimagLite::String lastSeen = "u";
		SizeType fcounter = 0; // counter of free
		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType ins = tensorSrep(i).ins();
			SizeType outs = tensorSrep(i).outs();
			RealType ysign = (tensorSrep(i).isConjugate()) ? -1.0 : 1.0;
			RealType dy = ysign*dy0;
			PairStringSizeType mypair = TensorSrep::splitIntoNameAndId(tensorSrep(i).fullName());
			PsimagLite::String label = mypair.first;

			if (lastSeen != label) ++layer;
			lastSeen = label;

			if (label == "u" || label == "h") {
				RealType factor = (label == "u") ? 1.0 : 0.5;
				PsimagLite::String style = (label == "u") ? "disen" : "ham";
				buffer_ += "\\coordinate (A";
				buffer_ += ttos(i) + ") at (" + ttos(x[i]) + "," + ttos(y[i]) + ");\n";
				buffer_ += "\\coordinate (B";
				buffer_ += ttos(i) + ") at (" + ttos(x[i] + dx);
				buffer_ += "," + ttos(y[i] + factor*dy) + ");\n";
				buffer_ += "\\draw[" + style + "] (A" + ttos(i) + ") rectangle (B";
				buffer_ += ttos(i)+ ");\n";
				// ins for u
				assert(ins > 1);
				RealType a = dx/(ins - 1);
				for (SizeType j = 0; j < ins; ++j) {
					SizeType k = absoluteLegNumber(i,j,ntensors);
					RealType xtmp = a*j + x[i];
					buffer_ += "\\coordinate (I" + label;
					buffer_ += ttos(k) + ") at (" + ttos(xtmp) + "," + ttos(y[i]) + ");\n";
					if (tensorSrep(i).legType(j) == TensorStanza::INDEX_TYPE_FREE) {
						buffer_ += ("\\coordinate (I" + label + "F");
						buffer_ += ttos(k) + ") at (" + ttos(xtmp) + ",";
						buffer_ += ttos(y[i]-0.5*dy) + ");\n";
						buffer_ += "\\draw[myfreelink] (I" + label + ttos(k);
						buffer_ += ") -- (I" + label + "F" + ttos(k) + ");\n";
						buffer_ += "\\node at ($(I" + label + "F" + ttos(k) +
						        ") + (0, -0.5)$) {\\large $\\sigma_{" +
						        ttos(fcounter++) + "}$};";
					}
				}

				// outs for u
				assert(outs > 0);
				RealType divisor = (outs == 1) ? 1 : outs - 1;
				a = dx/divisor;
				RealType doy = (label == "u") ? dy : 0.5*dy;
				for (SizeType j = 0; j < outs; ++j) {
					SizeType k = absoluteLegNumber(i, j, ntensors);
					RealType xtmp = a*j + x[i];
					if (i == 0) xtmp += dx;
					buffer_ += "\\coordinate (O" + label;
					buffer_ += ttos(k) + ") at (" + ttos(xtmp) + ",";
					buffer_ += ttos(y[i] + doy) + ");\n";
					if (tensorSrep(i).legType(j + ins) == TensorStanza::INDEX_TYPE_FREE) {
						buffer_ += ("\\coordinate (O" + label + "F");
						buffer_ += ttos(k) + ") at (" + ttos(xtmp) + ",";
						buffer_ += ttos(y[i]+1.5*dy) + ");\n";
						buffer_ += "\\draw[myfreelink] (O" + label + ttos(k);
						buffer_ += ") -- (O" + label + "F" + ttos(k) + ");\n";
					}
				}
			} else if (label == "w") {
				buffer_ += "\\coordinate (A";
				buffer_ += ttos(i) + ") at (" + ttos(x[i]) + "," + ttos(y[i]) + ");\n";
				buffer_ += "\\coordinate (B";
				buffer_ += ttos(i) + ") at (" + ttos(x[i]+dx);
				buffer_ += "," + ttos(y[i]) + ");\n";
				buffer_ += "\\coordinate (C";
				buffer_ += ttos(i) + ") at (" + ttos(x[i]+0.5*dx) + ",";
				buffer_ += ttos(y[i]+dy) + ");\n";
				buffer_ += "\\draw[isom] (A" + ttos(i) + ") -- (B" + ttos(i) + ") -- ";
				buffer_ += "(C" + ttos(i) + ") -- cycle;\n";

				// ins for w
				assert(ins > 0);
				RealType divisor = (ins == 1) ? 1 : ins - 1;
				RealType a = dx/divisor;
				for (SizeType j = 0; j < ins; ++j) {
					SizeType k = absoluteLegNumber(i,j,ntensors);
					RealType xtmp = a*j + x[i];
					buffer_ += "\\coordinate (I" + label;
					buffer_ += ttos(k) + ") at (" + ttos(xtmp) + "," + ttos(y[i]) + ");\n";
					if (tensorSrep(i).legType(j) == TensorStanza::INDEX_TYPE_FREE) {
						buffer_ += ("\\coordinate (I" + label + "F");
						buffer_ += ttos(k) + ") at (" + ttos(xtmp) + ",";
						buffer_ += ttos(y[i]-0.5*dy-1.5*ysign) + ");\n";
						buffer_ += "\\draw[myfreelink] (I" + label + ttos(k) + ") -- ";
						buffer_ += "(I" + label + "F" + ttos(k) + ");\n";
					}
				}

				// outs for w
				if (outs == 0) continue;
				if (outs > 1) {
					PsimagLite::String str("Dont' know how to draw a w with outs > 1\n");
					throw PsimagLite::RuntimeError(str);
				}

				RealType xtmp = x[i] + 0.5*dx;
				buffer_ += "\\coordinate (O" + label;
				SizeType k = absoluteLegNumber(i,0,ntensors);
				buffer_ +=  ttos(k) + ") at (" + ttos(xtmp) + "," + ttos(y[i]+dy) + ");\n";
			} else if (label == "r") {
				buffer_ += "\\coordinate (A";
				buffer_ += ttos(i) + ") at (" + ttos(x[i]) + "," + ttos(y[i]) + ");\n";
				buffer_ += "\\draw[fill=myblue] (A" + ttos(i) + ") circle (0.5);\n";
				for (SizeType j = 0; j < ins; ++j) {
					SizeType k = absoluteLegNumber(i,j,ntensors);
					buffer_ += "\\coordinate (I" + label;
					buffer_ += ttos(k) + ") at (" + ttos(x[i]) + "," + ttos(y[i]) + ");\n";
				}
			}

			//buffer_ += "\\node at (A" + ttos(i) + ") {" + ttos(i) + ",";
			//buffer_ += ttos(tensorSrep(i).id()) + "};\n";
		}

		drawConnections(tensorSrep);
	}

	void computeCoordinates(VectorRealType& x,
	                        VectorRealType& y,
	                        RealType dx,
	                        const TensorSrep& tensorSrep) const
	{
		SizeType ntensors = tensorSrep.size();

		SizeType yoffset0 = dx;
		VectorRealType savedXForR(2,0);
		SizeType modeForSavingForR = 0;
		for (SizeType i = 0; i < ntensors; ++i) {
			if (tensorSrep(i).type() == TensorStanza::TENSOR_TYPE_ERASED)
				continue;
			PairStringSizeType mypair = TensorSrep::splitIntoNameAndId(tensorSrep(i).fullName());
			PsimagLite::String name = mypair.first;
			if (mypair.second >= unpackTimeAndSpace_.size())
				continue;
			PairSizeType tensorXY = unpackTimeAndSpace_[mypair.second];
			SizeType tensorX = tensorXY.first;
			SizeType tensorY = tensorXY.second;

			RealType xwoffset = 0.5*dx - tensorY*dx*1.2;
			if (tensorY >= 2) xwoffset -= 1.1*tensorX*dx + 3*dx;
			RealType xuoffset = (tensorY >= 2) ? -6*dx + 3*dx*tensorX : 0.0;
			RealType ysign = (tensorSrep(i).isConjugate()) ? -1.0 : 1.0;
			RealType xsep = 3.0*dx*(1+tensorY);
			RealType xwsign = (tensorY & 1) ? -1 : 1;
			RealType xoffset = 3.0*pow(2,tensorY);
			if (tensorX == 0 && tensorY > 0 && name == "u") {
				SizeType j = findTensor(tensorSrep, tensorSrep(i).fullName());
				if (tensorSrep(j).type() == TensorStanza::TENSOR_TYPE_ERASED)
					continue;
			}

			if (name == "u") {
				x[i] = xsep*dx*tensorX + xoffset + xuoffset;
				y[i] = 3.5*tensorY*ysign + yoffset0*ysign;
			} else if (name == "w") {
				x[i] = xsep*dx*tensorX  + xoffset + pow(2,tensorY)*xwsign + xwoffset;
				y[i] = ysign*(3.5*tensorY + 1.5) + yoffset0*ysign;
				savedXForR[modeForSavingForR] = x[i];
				modeForSavingForR = (modeForSavingForR == 0) ? 1 : 0;
			} else if (name == "r") {
				// x[i] = 1.5*tauMax_*xsep;
				x[i] = (savedXForR[0] + savedXForR[1])*0.5;
				y[i] = ysign*3.5*tauMax_ + yoffset0*ysign;
			} else if (name == "h") {
				RealType tmp = findXForH(i,tensorSrep);
				x[i] = xsep*dx*tmp + xoffset;
				y[i] = 0;
			}
		}
	}

	SizeType tensorsAtLayer(SizeType layer,
	                        SizeType type,
	                        const TensorSrep& tensorSrep) const
	{
		SizeType counter = 0;
		SizeType ntensors = tensorSrep.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			TensorStanza::TensorTypeEnum t = tensorSrep(i).type();
			if (t != type) continue;
			TensorSrep::PairStringSizeType mypair = TensorSrep::
			        splitIntoNameAndId(tensorSrep(i).fullName());

			assert(mypair.second < unpackTimeAndSpace_.size());
			PairSizeType layerXY = unpackTimeAndSpace_[mypair.second];
			if (layerXY.second != layer) continue;
			counter++;
		}

		return counter;
	}

	SizeType findTensor(const TensorSrep& tensorSrep,
	                    PsimagLite::String fullName) const
	{
		SizeType ntensors = tensorSrep.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			if (tensorSrep(i).fullName() != fullName) continue;
			return i;
		}

		return ntensors;
	}

	SizeType absoluteLegNumber(SizeType ind, SizeType jnd, SizeType ntensors) const
	{
		return ind + jnd*ntensors;
	}

	void drawConnections(const TensorSrep& tensorSrep)
	{
		SizeType ntensors = tensorSrep.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType ins = tensorSrep(i).ins();
			for (SizeType j = 0; j < ins; ++j) {
				if (tensorSrep(i).legType(j) != TensorStanza::INDEX_TYPE_SUMMED) continue;
				TensorSrep::PairStringSizeType mypair = TensorSrep::
				        splitIntoNameAndId(tensorSrep(i).fullName());
				SizeType k1 = absoluteLegNumber(i,j,ntensors);
				PsimagLite::String label1 = "(I";
				label1 += mypair.first;
				label1 += ttos(k1) + ")";

				SizeType what = tensorSrep(i).legTag(j);
				PairSizeType k = findTarget(tensorSrep,i,what,TensorStanza::INDEX_DIR_OUT);
				if (k.first < ntensors) {
					PsimagLite::String label2 = "(O";
					TensorSrep::PairStringSizeType mypair2 = TensorSrep::
					        splitIntoNameAndId(tensorSrep(k.first).fullName());
					label2 += mypair2.first;
					SizeType k2 = absoluteLegNumber(k.first,k.second,ntensors);
					label2 += ttos(k2) + ")";
					buffer_ += "\\draw[mycon] " + label1 + " -- " + label2 + ";\n";
				}

				k = findTarget(tensorSrep,i,what,TensorStanza::INDEX_DIR_IN);
				if (k.first >= ntensors) continue;
				TensorSrep::PairStringSizeType mypair3 = TensorSrep::
				        splitIntoNameAndId(tensorSrep(k.first).fullName());
				PsimagLite::String label2 = "(I";
				label2 +=  mypair3.first;
				SizeType k2 = absoluteLegNumber(k.first,k.second,ntensors);
				label2 += ttos(k2) + ")";
				buffer_ += "\\draw[mycon] " + label1 + " -- " + label2 + ";\n";
			}
		}
	}

	PairSizeType findTarget(const TensorSrep& tensorSrep,
	                        SizeType selfInd,
	                        SizeType what,
	                        TensorStanza::IndexDirectionEnum dir) const
	{
		SizeType ntensors = tensorSrep.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			if (i == selfInd) continue;
			SizeType ins = tensorSrep(i).ins();
			SizeType s = (dir == TensorStanza::INDEX_DIR_IN) ? ins : tensorSrep(i).outs();
			SizeType offset = (dir == TensorStanza::INDEX_DIR_IN) ? 0 : ins;
			for (SizeType j = 0; j < s; ++j) {
				if (tensorSrep(i).legType(j + offset) != TensorStanza::INDEX_TYPE_SUMMED)
					continue;
				if (tensorSrep(i).legTag(j + offset) == what)
					return PairSizeType(i,j);
			}
		}

		return PairSizeType(ntensors,0);
	}

	void buildPacking(SizeType sites)
	{
		SizeType y = 0;
		SizeType stage = sites/2;
		while (stage > 1) {
			for (SizeType i = 0; i < stage; ++i)
				unpackTimeAndSpace_.push_back(PairSizeType(i,y));
			stage /= 2;
			y++;
		}
	}

	RealType findXForH(SizeType indexOfH,const TensorSrep& tSrep) const
	{
		VectorSizeType legs;
		findAllSummedLegTagOf(legs,tSrep(indexOfH));
		if (legs.size() == 0) {
			std::cerr<<"WARNING: findXForH failed\n";
			return 0;
		}

		RealType summedIds = 0.0;
		for (SizeType i = 0; i < legs.size(); ++i) {
			SizeType index = findTensorWithSummedLeg(legs[i],indexOfH,tSrep);
			TensorSrep::PairStringSizeType mypair = TensorSrep::
			        splitIntoNameAndId(tSrep(index).fullName());
			SizeType tmp =  mypair.second;
			summedIds += unpackTimeAndSpace_[tmp].first;
		}

		return summedIds/legs.size();
	}

	void findAllSummedLegTagOf(VectorSizeType& legs, const TensorStanza& stanza) const
	{
		SizeType legsNumber = stanza.legs();
		for (SizeType i = 0; i < legsNumber; ++i) {
			if (stanza.legType(i) != TensorStanza::INDEX_TYPE_SUMMED)
				continue;
			legs.push_back(stanza.legTag(i));
		}
	}

	SizeType findTensorWithSummedLeg(SizeType legTag,
	                                 SizeType indexOfH,
	                                 const TensorSrep& tSrep) const
	{
		SizeType ntensors = tSrep.size();
		for (SizeType j = 0; j < ntensors; ++j) {
			if (j == indexOfH) continue;
			const TensorStanza& stanza = tSrep(j);
			SizeType legs = stanza.legs();
			for (SizeType i = 0; i < legs; ++i) {
				if (stanza.legType(i) != TensorStanza::INDEX_TYPE_SUMMED)
					continue;
				if (legTag != stanza.legTag(i)) continue;
				return j;
			}
		}

		std::cerr<<"WARNING: findTensorWithSummedLeg failed\n";
		return 0;
	}

	static void printHeader(std::ostream& os)
	{
		PsimagLite::String str = "\\documentclass{standalone}\n";
		str += "\\usepackage{tikz}\n";
		str += "\\usepackage{pgfplots}\n";
		str += "\\pgfplotsset{width=7cm,compat=1.8}\n";
		str += "\\usetikzlibrary{calc}\n";
		str += "\\usetikzlibrary{positioning}\n";
		str += "\\usetikzlibrary{shadows}\n";
		str += "\\usetikzlibrary{shapes}\n";
		str += "\\usetikzlibrary{arrows}\n";
		str += "\\usepackage{xcolor}\n";
		str += "\\definecolor{myfuchsia}{HTML}{FF12BE}\n";
		str += "\\definecolor{myblue}{HTML}{0074D9}\n";
		str += "\\definecolor{mygreen}{HTML}{2ECC40}\n";
		str += "\\definecolor{myyellow}{HTML}{FFDC00}\n";
		str += "\\definecolor{mygray}{HTML}{AAAAAA}\n";
		str += "\\definecolor{notbisque}{HTML}{007733}\n";
		str += "\\definecolor{notdarksalmon}{HTML}{B0E97A}\n";
		str += "\\definecolor{notcrimson}{HTML}{FF4136}\n";
		str += "\\definecolor{notfirebrick}{HTML}{48B221}\n";
		str += "\\definecolor{notdarkred}{HTML}{004000}\n";
		str += "\\definecolor{notCoral}{HTML}{00C0C0}\n";

		str += "\\begin{document}";
		str += "\\begin{tikzpicture}[\n";
		str += "mycon/.style={draw=mygray,thick},\n";
		str += "disen/.style={fill=notfirebrick,draw=notfirebrick},\n";
		str += "isom/.style={fill=notcrimson,draw=notcrimson},\n";
		str += "myfreelink/.style={very thick,notCoral},\n";
		str += "ham/.style={fill=myyellow}\n";
		str += "]";
		os<<str<<"\n";
	}

	static void printFooter(std::ostream& os)
	{
		PsimagLite::String str = "\\end{tikzpicture}\n";
		str += "\\end{document}\n";
		os<<str;
	}

	static PsimagLite::String buffer_;
	PsimagLite::String srep_;
	SizeType tauMax_;
	VectorPairSizeType unpackTimeAndSpace_;
}; // class MeraToTikz

template<typename T>
PsimagLite::String MeraToTikz<T>::buffer_;
} // namespace Mera
#endif // MERATOTIKZ_H
