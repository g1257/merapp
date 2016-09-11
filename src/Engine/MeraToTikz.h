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

		buffer_ = "%" + srep_;
		buffer_ += "\n";

		RealType dx = 1.;
		RealType dy = 1.;

		SizeType ntensors = tensorSrep.size();
		VectorRealType x(ntensors);
		VectorRealType y(ntensors);
		computeCoordinates(x,y,dx,dy,tensorSrep);

		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType ins = tensorSrep(i).ins();
			SizeType outs = tensorSrep(i).outs();

			if (tensorSrep(i).type() == TensorStanza::TENSOR_TYPE_U) {
				buffer_ += "\\coordinate (A";
				buffer_ += ttos(i) + ") at (" + ttos(x[i]) + "," + ttos(y[i]) + ");\n";
				buffer_ += "\\coordinate (B";
				buffer_ += ttos(i) + ") at (" + ttos(x[i] + dx);
				buffer_ += "," + ttos(y[i] + dy) + ");\n";
				buffer_ += "\\draw[disen] (A" + ttos(i) + ") rectangle (B";
				buffer_ += ttos(i)+ ");\n";
				// ins for u
				assert(ins > 0);
				RealType a = dx/(ins - 1);
				for (SizeType j = 0; j < ins; ++j) {
					SizeType k = absoluteLegNumber(i,j,ntensors);
					RealType xtmp = a*j + x[i];
					buffer_ += "\\coordinate (IU";
					buffer_ += ttos(k) + ") at (" + ttos(xtmp) + "," + ttos(y[i]) + ");\n";
					if (tensorSrep(i).legType(j,TensorStanza::INDEX_DIR_IN) ==
					        TensorStanza::INDEX_TYPE_FREE) {
						buffer_ += ("\\coordinate (IUF");
						buffer_ += ttos(k) + ") at (" + ttos(xtmp) + ",";
						buffer_ += ttos(y[i]-0.5*dy) + ");\n";
						buffer_ += "\\draw (IU" + ttos(k) + ") -- (IUF" + ttos(k) + ");\n";
					}
				}

				// outs for u
				a = dx/(outs - 1);
				for (SizeType j = 0; j < outs; ++j) {
					SizeType k = absoluteLegNumber(i,j,ntensors);
					RealType xtmp = a*j + x[i];
					buffer_ += "\\coordinate (OU";
					buffer_ += ttos(k) + ") at (" + ttos(xtmp) + "," + ttos(y[i] + dy) + ");\n";
				}
			} else {
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
					buffer_ += "\\coordinate (IW";
					buffer_ += ttos(k) + ") at (" + ttos(xtmp) + "," + ttos(y[i]) + ");\n";
					if (tensorSrep(i).legType(j,TensorStanza::INDEX_DIR_IN) ==
					        TensorStanza::INDEX_TYPE_FREE) {
						buffer_ += ("\\coordinate (IWF");
						buffer_ += ttos(k) + ") at (" + ttos(xtmp) + ",";
						buffer_ += ttos(y[i]-0.5*dy-1.5) + ");\n";
						buffer_ += "\\draw (IW" + ttos(k) + ") -- (IWF" + ttos(k) + ");\n";
					}
				}

				// outs for w
				if (outs == 0) continue;
				if (outs > 1) {
					PsimagLite::String str("Dont' know how to draw a w with outs > 1\n");
					throw PsimagLite::RuntimeError(str);
				}

				RealType xtmp = x[i] + 0.5*dx;
				buffer_ += "\\coordinate (OW";
				SizeType k = absoluteLegNumber(i,0,ntensors);
				buffer_ +=  ttos(k) + ") at (" + ttos(xtmp) + "," + ttos(y[i]+dy) + ");\n";
			}

			buffer_ += "\\node at (A" + ttos(i) + ") {" + ttos(i) + ",";
			buffer_ += ttos(tensorSrep(i).id()) + "};\n";
		}

		drawConnections(tensorSrep);
	}

	void computeCoordinates(VectorRealType& x,
	                        VectorRealType& y,
	                        RealType dx,
	                        RealType dy,
	                        const TensorSrep& tensorSrep) const
	{
		const RealType slope = 1.0;
		SizeType ntensors = tensorSrep.size();
		RealType xoffset = 0.0;
		RealType xwoffset = 1.5*dx;
		SizeType firstWofLayer = ntensors;
		SizeType savedY = ntensors;
		for (SizeType i = 0; i < ntensors; ++i) {
			SizeType tensorX = 0;
			SizeType tensorY = 0;
			ProgramGlobals::unpackTimeAndSpace(tensorY,
			                                   tensorX,
			                                   tensorSrep(i).id(),
			                                   tauMax_);
			RealType xsep = 3.0*(dx+tensorY);
			SizeType type = tensorSrep(i).type();
			if (tensorX == 0 && tensorY > 0 && type == TensorStanza::TENSOR_TYPE_U) {
				// compute xshift
				std::cerr<<"firstWofLayer "<<firstWofLayer<<"\n";
				assert(firstWofLayer < ntensors);
				RealType yu = 3.5*tensorY;
				RealType xw = x[firstWofLayer];
				RealType yw = y[firstWofLayer];
				xoffset = xw + (yu - yw)/slope;
			}

			bool b = (tensorSrep(i).outs() == 0)
			        ? false : (tensorSrep(i).legType(0,TensorStanza::INDEX_DIR_OUT) ==
			                   TensorStanza::INDEX_TYPE_SUMMED);
			if (type == TensorStanza::TENSOR_TYPE_W && savedY != tensorY && b) {
				savedY = tensorY;
				firstWofLayer = i;
			}

			if (tensorX == 0 && tensorY > 0 && type == TensorStanza::TENSOR_TYPE_U) {
				SizeType id = tensorSrep(i).id();
				SizeType j = findTensor(tensorSrep,id,TensorStanza::TENSOR_TYPE_W);
				std::cerr<<"for i= "<<i<<" FOUND j="<<j<<"\n";
				if (tensorSrep(i).legTag(0,TensorStanza::INDEX_DIR_OUT) ==
				        tensorSrep(j).legTag(1,TensorStanza::INDEX_DIR_IN)) {
					xwoffset = -1.5*dx;
				} else {
					xwoffset = 1.5*dx;
				}
			}

			if (type == TensorStanza::TENSOR_TYPE_U) {
				x[i] = xsep*dx*tensorX + xoffset;
				y[i] = 3.5*tensorY;
			} else {
				assert(type == TensorStanza::TENSOR_TYPE_W);
				x[i] = xsep*dx*tensorX + xwoffset + xoffset;
				y[i] = 3.5*tensorY + 1.5;
			}
		}
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
			SizeType outs = tensorSrep(i).outs();
			TensorStanza::TensorTypeEnum type = tensorSrep(i).type();
			for (SizeType j = 0; j < outs; ++j) {
				if (tensorSrep(i).legType(j,TensorStanza::INDEX_DIR_OUT) !=
				        TensorStanza::INDEX_TYPE_SUMMED) continue;
				SizeType k1 = absoluteLegNumber(i,j,ntensors);
				PsimagLite::String label1 = "(O";
				label1 += (type == TensorStanza::TENSOR_TYPE_U) ? "U" : "W";
				label1 += ttos(k1) + ")";
				SizeType what = tensorSrep(i).legTag(j,TensorStanza::INDEX_DIR_OUT);
				PairSizeType k = findTarget(tensorSrep,what,type);
				if (k.first >= ntensors) continue;
				PsimagLite::String label2 = "(I";
				label2 += (tensorSrep(k.first).type() == TensorStanza::TENSOR_TYPE_U) ? "U" : "W";
				SizeType k2 = absoluteLegNumber(k.first,k.second,ntensors);
				label2 += ttos(k2) + ")";
				buffer_ += "\\draw " + label1 + " -- " + label2 + ";\n";
			}
		}
	}

	PairSizeType findTarget(const TensorSrep& tensorSrep,
	                        SizeType what,
	                        TensorStanza::TensorTypeEnum type) const
	{
		SizeType ntensors = tensorSrep.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			if (tensorSrep(i).type() == type) continue;
			SizeType ins = tensorSrep(i).ins();
			for (SizeType j = 0; j < ins; ++j) {
				if (tensorSrep(i).legType(j,TensorStanza::INDEX_DIR_IN) !=
				        TensorStanza::INDEX_TYPE_SUMMED) continue;
				if (tensorSrep(i).legTag(j,TensorStanza::INDEX_DIR_IN) == what)
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
		str += "\\begin{document}";
		str += "\\begin{tikzpicture}[\n";
		str += "disen/.style={fill=green},\n";
		str += "isom/.style={fill=myfuchsia},mylink2/.style={very thick}]";
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
