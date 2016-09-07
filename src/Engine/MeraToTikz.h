#ifndef MERATOTIKZ_H
#define MERATOTIKZ_H
#include "Vector.h"
#include "Tokenizer.h"

namespace Mera {

template<typename RealType>
class MeraToTikz {

	enum TensorTypeEnum {TENSOR_TYPE_UNKNOWN,TENSOR_TYPE_W,TENSOR_TYPE_U};

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	class OpaqueTensor {

	public:

		OpaqueTensor(PsimagLite::String srep)
		    :  conjugate_(false),
		      srep_(srep),name_(""),
		      type_(TENSOR_TYPE_UNKNOWN)
		{
			VectorStringType tokens;
			PsimagLite::tokenizer(srep_,tokens,":");
			SizeType ts = tokens.size();
			if (ts < 6) {
				PsimagLite::String str("OpaqueTensor: malformed partial srep ");
				throw PsimagLite::RuntimeError(str + srep_ + "\n");
			}

			name_ = tokens[0];
			SizeType l = name_.length();
			if (l == 0) {
				PsimagLite::String str("OpaqueTensor: malformed partial srep, empty name ");
				throw PsimagLite::RuntimeError(str + srep_ + "\n");
			}

			if (name_.substr(l-1,l) == "*") conjugate_ = true;
			SizeType l2 = (conjugate_) ? l - 1 : l;
			if (name_.substr(0,l2) == "u") type_ = TENSOR_TYPE_U;
			if (name_.substr(0,l2) == "w") type_ = TENSOR_TYPE_W;

			if (type_ == TENSOR_TYPE_UNKNOWN) {
				PsimagLite::String str("OpaqueTensor: partial srep, tensor type");
				throw PsimagLite::RuntimeError(str + srep_ + "\n");
			}

			timeIndex_ = atoi(tokens[1].c_str());
			spaceIndex_ = atoi(tokens[2].c_str());
			ins_ = atoi(tokens[3].c_str());
			outs_ = atoi(tokens[4].c_str());

			if (ins_ + outs_ + 5 != ts) {
				PsimagLite::String str("OpaqueTensor: malformed partial srep ");
				throw PsimagLite::RuntimeError(str + srep_ + "\n");
			}
		}

		SizeType x() const
		{
			return spaceIndex_;
		}

		SizeType y() const
		{
			return timeIndex_;
		}

		TensorTypeEnum type() const { return type_; }

	private:

		SizeType timeIndex_;
		SizeType spaceIndex_;
		SizeType ins_;
		SizeType outs_;
		bool conjugate_;
		PsimagLite::String srep_;
		PsimagLite::String name_;
		TensorTypeEnum type_;
	};

public:

	MeraToTikz(PsimagLite::String srep)
	    : srep_(srep)
	{
		fillBuffer();
	}

	~MeraToTikz()
	{
		for (SizeType i = 0; i < tensor_.size(); ++i) {
			delete tensor_[i];
			tensor_[i] = 0;
		}
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
		VectorStringType tokens;
		PsimagLite::tokenizer(srep_,tokens,";");
		tensor_.resize(tokens.size(),0);
		for (SizeType i = 0; i < tokens.size(); ++i)
			tensor_[i] = new OpaqueTensor(tokens[i]);

		buffer_ = "%" + srep_;
		buffer_ += "\n";

		RealType dx = 1.;
		RealType dy = 1.;
		for (SizeType i = 0; i < tensor_.size(); ++i) {
			buffer_ += ("\\coordinate (A");
			RealType x = 1.5*tensor_[i]->x();
			RealType y = 3.0*tensor_[i]->y();
			buffer_ += ttos(i) + ") at (" + ttos(x) + "," + ttos(y) + ");\n";
			buffer_ += ("\\coordinate (B");
			buffer_ += ttos(i) + ") at (" + ttos(x+dx);
			buffer_ += "," + ttos(y+dy) + ");\n";
			if (tensor_[i]->type() == TENSOR_TYPE_U) {
				buffer_ += "\\draw[disen] (A" + ttos(i) + ") rectangle (B";
				buffer_ += ttos(i)+ ");\n";
			} else {
				//buffer_ += "\\draw (A" + ttos(i) + ") circle (0.5);\n";
			}
		}
	}

	static void printHeader(std::ostream& os)
	{
		PsimagLite::String str = "\\documentclass{article}\n";
		str += "\\usepackage{tikz}\n";
		str += "\\usepackage{pgfplots}\n";
		str += "\\pgfplotsset{width=7cm,compat=1.8}\n";
		str += "\\usetikzlibrary{positioning}\n";
		str += "\\usetikzlibrary{shadows}\n";
		str += "\\usetikzlibrary{shapes}\n";
		str += "\\usetikzlibrary{arrows}\n";
		str += "\\usepackage{xcolor}\n";

		str += "\\begin{document}";
		str += "\\begin{tikzpicture}[\n";
		str += "disen/.style={fill=green},mylink2/.style={very thick}]";
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
	typename PsimagLite::Vector<OpaqueTensor*>::Type tensor_;
}; // class MeraToTikz

template<typename T>
PsimagLite::String MeraToTikz<T>::buffer_;
} // namespace Mera
#endif // MERATOTIKZ_H
