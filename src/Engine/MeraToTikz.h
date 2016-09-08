#ifndef MERATOTIKZ_H
#define MERATOTIKZ_H
#include "Vector.h"
#include "Tokenizer.h"
#include "TypeToString.h"
#include <cassert>
#include <iostream>

namespace Mera {

template<typename RealType>
class MeraToTikz {

	enum TensorTypeEnum {TENSOR_TYPE_UNKNOWN,TENSOR_TYPE_W,TENSOR_TYPE_U};

	enum IndexDirectionEnum {INDEX_DIR_IN, INDEX_DIR_OUT};

	enum IndexTypeEnum {INDEX_TYPE_SUMMED, INDEX_TYPE_FREE};

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef std::pair<char,SizeType> PairCharSizeType;
	typedef PsimagLite::Vector<PairCharSizeType>::Type VectorPairCharSizeType;

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
			SizeType ins = atoi(tokens[3].c_str());
			SizeType outs = atoi(tokens[4].c_str());

			if (ins + outs + 5 != ts) {
				PsimagLite::String str("OpaqueTensor: malformed partial srep ");
				throw PsimagLite::RuntimeError(str + srep_ + "\n");
			}

			for (SizeType i = 0; i < ins; ++i) {
				PsimagLite::String item = tokens[5+i];
				l = item.length();
				if (l == 0) {
					PsimagLite::String str("OpaqueTensor: malformed srep, in s-index ");
					throw PsimagLite::RuntimeError(str + ttos(i) + ", " + srep_ + "\n");
				}

				if (l == 1 && item != "d") {
					PsimagLite::String str("OpaqueTensor: malformed srep, in s-index ");
					throw PsimagLite::RuntimeError(str + ttos(i) + ", " + srep_ + "\n");
				}

				PsimagLite::String tmp = item.substr(1,l);
				insSi_.push_back(PairCharSizeType(item[0],atoi(tmp.c_str())));
			}

			for (SizeType i = 0; i < outs; ++i) {
				PsimagLite::String item = tokens[5+ins+i];
				l = item.length();

				if (l == 0) {
					PsimagLite::String str("OpaqueTensor: malformed srep, out s-index ");
					throw PsimagLite::RuntimeError(str + ttos(i) + ", " + srep_ + "\n");
				}

				if (l == 1 && item != "d") {
					PsimagLite::String str("OpaqueTensor: malformed srep, out s-index ");
					throw PsimagLite::RuntimeError(str + ttos(i) + ", " + srep_ + "\n");
				}

				PsimagLite::String tmp = (l == 1) ? "0" : item.substr(1,l);
				outsSi_.push_back(PairCharSizeType(item[0],atoi(tmp.c_str())));
			}
		}

		SizeType x() const { return spaceIndex_; }

		SizeType y() const { return timeIndex_; }

		SizeType ins() const { return insSi_.size(); }

		SizeType outs() const { return outsSi_.size(); }

		IndexTypeEnum indexType(SizeType ind, IndexDirectionEnum dir) const
		{
			char c = (dir == INDEX_DIR_IN) ? insSi_[ind].first : outsSi_[ind].first;
			return (c == 's') ? INDEX_TYPE_SUMMED : INDEX_TYPE_FREE;
		}

		TensorTypeEnum type() const { return type_; }

	private:

		SizeType timeIndex_;
		SizeType spaceIndex_;
		bool conjugate_;
		PsimagLite::String srep_;
		PsimagLite::String name_;
		TensorTypeEnum type_;
		VectorPairCharSizeType insSi_;
		VectorPairCharSizeType outsSi_;
	};

public:

	MeraToTikz(PsimagLite::String srep)
	    : srep_(srep)
	{
		cleanWhiteSpace(srep_);
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

			RealType xsep = 1.5*(1.0+tensor_[i]->y());
			RealType xdisen = xsep*tensor_[i]->x() + 1.5*tensor_[i]->y();
			RealType ydisen = 3.5*tensor_[i]->y();
			RealType xisom = xdisen + 1;
			RealType yisom = ydisen + 1.5;
			SizeType ins = tensor_[i]->ins();
			if (tensor_[i]->type() == TENSOR_TYPE_U) {
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
					buffer_ += "\\coordinate (I";
					buffer_ += ttos(i) + ") at (" + ttos(xtmp) + "," + ttos(ydisen) + ");\n";
					if (tensor_[i]->indexType(j,INDEX_DIR_IN) == INDEX_TYPE_FREE) {
						buffer_ += ("\\coordinate (IF");
						buffer_ += ttos(i) + ") at (" + ttos(xtmp) + ",";
						buffer_ += ttos(ydisen-0.5*dy) + ");\n";
						buffer_ += "\\draw (I" + ttos(i) + ") -- (IF" + ttos(i) + ");\n";
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
	typename PsimagLite::Vector<OpaqueTensor*>::Type tensor_;
}; // class MeraToTikz

template<typename T>
PsimagLite::String MeraToTikz<T>::buffer_;
} // namespace Mera
#endif // MERATOTIKZ_H
