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
#include "Tokenizer.h"
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

public:

	MeraToTikz(PsimagLite::String srep, SizeType tauMax)
	    : srep_(srep),tauMax_(tauMax)
	{
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

		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType ins = tensorSrep(i).ins();
			SizeType outs = tensorSrep(i).outs();
			RealType ysign = (tensorSrep(i).isConjugate()) ? -1.0 : 1.0;
			RealType dy = ysign*dy0;
			PsimagLite::String label = tensorSrep(i).label();

			if (tensorSrep(i).type() == TensorStanza::TENSOR_TYPE_U ||
			       tensorSrep(i).type() == TensorStanza::TENSOR_TYPE_H) {
				RealType factor = (tensorSrep(i).type() == TensorStanza::TENSOR_TYPE_U) ?
				            1.0 : 0.5;
				PsimagLite::String style = (tensorSrep(i).type() == TensorStanza::TENSOR_TYPE_U) ?
				            "disen" : "ham";
				buffer_ += "\\coordinate (A";
				buffer_ += ttos(i) + ") at (" + ttos(x[i]) + "," + ttos(y[i]) + ");\n";
				buffer_ += "\\coordinate (B";
				buffer_ += ttos(i) + ") at (" + ttos(x[i] + dx);
				buffer_ += "," + ttos(y[i] + factor*dy) + ");\n";
				buffer_ += "\\draw[" + style + "] (A" + ttos(i) + ") rectangle (B";
				buffer_ += ttos(i)+ ");\n";
				// ins for u
				assert(ins > 0);
				RealType a = dx/(ins - 1);
				for (SizeType j = 0; j < ins; ++j) {
					SizeType k = absoluteLegNumber(i,j,ntensors);
					RealType xtmp = a*j + x[i];
					buffer_ += "\\coordinate (I" + label;
					buffer_ += ttos(k) + ") at (" + ttos(xtmp) + "," + ttos(y[i]) + ");\n";
					if (tensorSrep(i).legType(j,TensorStanza::INDEX_DIR_IN) ==
					        TensorStanza::INDEX_TYPE_FREE) {
						buffer_ += ("\\coordinate (I" + label + "F");
						buffer_ += ttos(k) + ") at (" + ttos(xtmp) + ",";
						buffer_ += ttos(y[i]-0.5*dy) + ");\n";
						buffer_ += "\\draw[myfreelink] (I" + label + ttos(k);
						buffer_ += ") -- (I" + label + "F" + ttos(k) + ");\n";
					}
				}

				// outs for u
				a = dx/(outs - 1);
				for (SizeType j = 0; j < outs; ++j) {
					SizeType k = absoluteLegNumber(i,j,ntensors);
					RealType xtmp = a*j + x[i];
					buffer_ += "\\coordinate (O" + label;
					buffer_ += ttos(k) + ") at (" + ttos(xtmp) + "," + ttos(y[i] + dy) + ");\n";
					if (tensorSrep(i).legType(j,TensorStanza::INDEX_DIR_OUT) ==
					        TensorStanza::INDEX_TYPE_FREE) {
						buffer_ += ("\\coordinate (O" + label + "F");
						buffer_ += ttos(k) + ") at (" + ttos(xtmp) + ",";
						buffer_ += ttos(y[i]+1.5*dy) + ");\n";
						buffer_ += "\\draw[myfreelink] (O" + label + ttos(k);
						buffer_ += ") -- (O" + label + "F" + ttos(k) + ");\n";
					}
				}
			} else if (tensorSrep(i).type() == TensorStanza::TENSOR_TYPE_W) {
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
				RealType a = dx/(ins - 1);
				for (SizeType j = 0; j < ins; ++j) {
					SizeType k = absoluteLegNumber(i,j,ntensors);
					RealType xtmp = a*j + x[i];
					buffer_ += "\\coordinate (I" + label;
					buffer_ += ttos(k) + ") at (" + ttos(xtmp) + "," + ttos(y[i]) + ");\n";
					if (tensorSrep(i).legType(j,TensorStanza::INDEX_DIR_IN) ==
					        TensorStanza::INDEX_TYPE_FREE) {
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
			} else if (tensorSrep(i).type() == TensorStanza::TENSOR_TYPE_ROOT) {
				buffer_ += "\\coordinate (A";
				buffer_ += ttos(i) + ") at (" + ttos(x[i]) + "," + ttos(y[i]) + ");\n";
				buffer_ += "\\draw[fill=myblue] (A" + ttos(i) + ") circle (0.5);\n";
				for (SizeType j = 0; j < ins; ++j) {
					SizeType k = absoluteLegNumber(i,j,ntensors);
					buffer_ += "\\coordinate (I" + label;
					buffer_ += ttos(k) + ") at (" + ttos(x[i]) + "," + ttos(y[i]) + ");\n";
				}
			}

			buffer_ += "\\node at (A" + ttos(i) + ") {" + ttos(i) + ",";
			buffer_ += ttos(tensorSrep(i).id()) + "};\n";
		}

		drawConnections(tensorSrep);
	}

	void computeCoordinates(VectorRealType& x,
	                        VectorRealType& y,
	                        RealType dx,
	                        const TensorSrep& tensorSrep) const
	{
		SizeType ntensors = tensorSrep.size();
		RealType xwoffset = 1.5*dx;
		SizeType yoffset0 = dx;

		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType tensorX = 0;
			SizeType tensorY = 0;
			ProgramGlobals::unpackTimeAndSpace(tensorY,
			                                   tensorX,
			                                   tensorSrep(i).id(),
			                                   tauMax_);
			SizeType type = tensorSrep(i).type();
			RealType ysign = (tensorSrep(i).isConjugate()) ? -1.0 : 1.0;
			RealType xsep = 3.0*dx*(1+tensorY);
			RealType xoffset = 3.0*pow(2,tensorY);
			if (tensorX == 0 && tensorY > 0 && type == TensorStanza::TENSOR_TYPE_U) {
				SizeType id = tensorSrep(i).id();
				SizeType j = findTensor(tensorSrep,id,TensorStanza::TENSOR_TYPE_W);
				if (tensorSrep(i).legTag(0,TensorStanza::INDEX_DIR_OUT) ==
				        tensorSrep(j).legTag(1,TensorStanza::INDEX_DIR_IN)) {
					xwoffset = -1.5*dx;
				} else {
					xwoffset = 1.5*dx;
				}
			}

			if (type == TensorStanza::TENSOR_TYPE_U) {
				x[i] = xsep*dx*tensorX + xoffset;
				y[i] = 3.5*tensorY*ysign + yoffset0*ysign;
			} else if (type == TensorStanza::TENSOR_TYPE_W) {
				x[i] = xsep*dx*tensorX + xwoffset + xoffset;
				y[i] = ysign*(3.5*tensorY + 1.5) + yoffset0*ysign;
			} else if (type == TensorStanza::TENSOR_TYPE_ROOT) {
				x[i] = 1.5*tauMax_*xsep;
				y[i] = ysign*3.5*tauMax_ + yoffset0*ysign;
			} else if (type == TensorStanza::TENSOR_TYPE_H) {
				x[i] = xsep*dx*tensorX + xoffset;
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
			SizeType layerX = 0;
			SizeType layerY = 0;
			ProgramGlobals::unpackTimeAndSpace(layerY,layerX,tensorSrep(i).id(),tauMax_);
			if (layerY != layer) continue;
			counter++;
		}

		return counter;
	}

	SizeType findTensor(const TensorSrep& tensorSrep,
	                    SizeType id,
	                    TensorStanza::TensorTypeEnum type) const
	{
		SizeType ntensors = tensorSrep.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			TensorStanza::TensorTypeEnum t = tensorSrep(i).type();
			if (t != type) continue;
			if (tensorSrep(i).id() != id) continue;
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
				if (tensorSrep(i).legType(j,TensorStanza::INDEX_DIR_IN) !=
				        TensorStanza::INDEX_TYPE_SUMMED) continue;
				SizeType k1 = absoluteLegNumber(i,j,ntensors);
				PsimagLite::String label1 = "(I";
				label1 += tensorSrep(i).label();
				label1 += ttos(k1) + ")";

				SizeType what = tensorSrep(i).legTag(j,TensorStanza::INDEX_DIR_IN);
				PairSizeType k = findTarget(tensorSrep,i,what,TensorStanza::INDEX_DIR_OUT);
				if (k.first < ntensors) {
					PsimagLite::String label2 = "(O";
					label2 += tensorSrep(k.first).label();
					SizeType k2 = absoluteLegNumber(k.first,k.second,ntensors);
					label2 += ttos(k2) + ")";
					buffer_ += "\\draw " + label1 + " -- " + label2 + ";\n";
				}

				k = findTarget(tensorSrep,i,what,TensorStanza::INDEX_DIR_IN);
				if (k.first >= ntensors) continue;
				PsimagLite::String label2 = "(I";
				label2 += tensorSrep(k.first).label();
				SizeType k2 = absoluteLegNumber(k.first,k.second,ntensors);
				label2 += ttos(k2) + ")";
				buffer_ += "\\draw " + label1 + " -- " + label2 + ";\n";
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
			SizeType s = (dir == TensorStanza::INDEX_DIR_IN) ?
			            tensorSrep(i).ins() : tensorSrep(i).outs();
			for (SizeType j = 0; j < s; ++j) {
				if (tensorSrep(i).legType(j,dir) !=
				        TensorStanza::INDEX_TYPE_SUMMED) continue;
				if (tensorSrep(i).legTag(j,dir) == what)
					return PairSizeType(i,j);
			}
		}

		return PairSizeType(ntensors,0);
	}

	static void printHeader(std::ostream& os)
	{
		PsimagLite::String str = "\\documentclass{standalone}\n";
		str += "\\usepackage{tikz}\n";
		str += "\\usepackage{pgfplots}\n";
		str += "\\pgfplotsset{width=7cm,compat=1.8}\n";
		str += "\\usetikzlibrary{positioning}\n";
		str += "\\usetikzlibrary{shadows}\n";
		str += "\\usetikzlibrary{shapes}\n";
		str += "\\usetikzlibrary{arrows}\n";
		str += "\\usepackage{xcolor}\n";
		str += "\\definecolor{myfuchsia}{HTML}{FF12BE}\n";
		str += "\\definecolor{myblue}{HTML}{0074D9}\n";
		str += "\\definecolor{mygreen}{HTML}{2ECC40}\n";
		str += "\\definecolor{myyellow}{HTML}{FFDC00}\n";
		str += "\\begin{document}";
		str += "\\begin{tikzpicture}[\n";
		str += "disen/.style={fill=mygreen},\n";
		str += "isom/.style={fill=myfuchsia},\n";
		str += "myfreelink/.style={thick,mygreen},\n";
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
}; // class MeraToTikz

template<typename T>
PsimagLite::String MeraToTikz<T>::buffer_;
} // namespace Mera
#endif // MERATOTIKZ_H
