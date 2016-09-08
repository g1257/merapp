#ifndef MERATOTIKZ_H
#define MERATOTIKZ_H
#include "Vector.h"
#include "Tokenizer.h"
#include "TypeToString.h"
#include <cassert>
#include <iostream>
#include "TensorSrep.h"

namespace Mera {

template<typename RealType>
class MeraToTikz {

public:

	MeraToTikz(PsimagLite::String srep)
	    : srep_(srep)
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
		for (SizeType i = 0; i < tensorSrep.size(); ++i) {

			RealType xsep = 3.0*(dx+tensorSrep(i).y());
			RealType xdisen = xsep*dx*tensorSrep(i).x() + 1.5*tensorSrep(i).y();
			RealType ydisen = 3.5*tensorSrep(i).y();
			RealType xisom = xdisen + 1.5*dx;
			RealType yisom = ydisen + 1.5;
			SizeType ins = tensorSrep(i).ins();
			if (tensorSrep(i).type() == TensorStanza::TENSOR_TYPE_U) {
				buffer_ += "\\coordinate (A";
				buffer_ += ttos(i) + ") at (" + ttos(xdisen) + "," + ttos(ydisen) + ");\n";
				buffer_ += "\\coordinate (B";
				buffer_ += ttos(i) + ") at (" + ttos(xdisen+dx);
				buffer_ += "," + ttos(ydisen+dy) + ");\n";
				buffer_ += "\\draw[disen] (A" + ttos(i) + ") rectangle (B";
				buffer_ += ttos(i)+ ");\n";
				// ins
				assert(ins > 0);
				RealType a = dx/(ins - 1);
				for (SizeType j = 0; j < ins; ++j) {
					RealType xtmp = a*j + xdisen;
					buffer_ += "\\coordinate (IU";
					buffer_ += ttos(i) + ") at (" + ttos(xtmp) + "," + ttos(ydisen) + ");\n";
					if (tensorSrep(i).indexType(j,TensorStanza::INDEX_DIR_IN) ==
					        TensorStanza::INDEX_TYPE_FREE) {
						buffer_ += ("\\coordinate (IUF");
						buffer_ += ttos(i) + ") at (" + ttos(xtmp) + ",";
						buffer_ += ttos(ydisen-0.5*dy) + ");\n";
						buffer_ += "\\draw (IU" + ttos(i) + ") -- (IUF" + ttos(i) + ");\n";
					}
				}
			} else {
				buffer_ += "\\coordinate (A";
				buffer_ += ttos(i) + ") at (" + ttos(xisom) + "," + ttos(yisom) + ");\n";
				buffer_ += "\\coordinate (B";
				buffer_ += ttos(i) + ") at (" + ttos(xisom+dx);
				buffer_ += "," + ttos(yisom) + ");\n";
				buffer_ += "\\coordinate (C";
				buffer_ += ttos(i) + ") at (" + ttos(xisom+0.5*dx) + ",";
				buffer_ += ttos(yisom+dy) + ");\n";
				buffer_ += "\\draw[isom] (A" + ttos(i) + ") -- (B" + ttos(i) + ") -- ";
				buffer_ += "(C" + ttos(i) + ") -- cycle;\n";
				RealType a = dx/(ins - 1);
				for (SizeType j = 0; j < ins; ++j) {
					RealType xtmp = a*j + xisom;
					buffer_ += "\\coordinate (IW";
					buffer_ += ttos(i) + ") at (" + ttos(xtmp) + "," + ttos(yisom) + ");\n";
					if (tensorSrep(i).indexType(j,TensorStanza::INDEX_DIR_IN) ==
					        TensorStanza::INDEX_TYPE_FREE) {
						buffer_ += ("\\coordinate (IWF");
						buffer_ += ttos(i) + ") at (" + ttos(xtmp) + ",";
						buffer_ += ttos(ydisen-0.5*dy) + ");\n";
						buffer_ += "\\draw (IW" + ttos(i) + ") -- (IWF" + ttos(i) + ");\n";
					}
				}
			}
		}
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
}; // class MeraToTikz

template<typename T>
PsimagLite::String MeraToTikz<T>::buffer_;
} // namespace Mera
#endif // MERATOTIKZ_H
