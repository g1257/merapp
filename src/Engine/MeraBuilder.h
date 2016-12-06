#ifndef MERABUILDER_H
#define MERABUILDER_H

namespace Mera {

class MeraBuilder {

public:

	MeraBuilder(SizeType sites,SizeType arity,SizeType dimension)
	    : srep_("")
	{
		if (dimension != 1)
			throw PsimagLite::RuntimeError("MeraBuilder: dimension must be 1 for now\n");

		if (arity != 2)
			throw PsimagLite::RuntimeError("MeraBuilder: arity must be 2 for now\n");

		SizeType ln = 0;
		if ((ln = logBase2Strict(sites)) == 0)
			throw PsimagLite::RuntimeError("MeraBuilder: sites must be a power of 2\n");

		SizeType tensors = sites/2;
		srep_ = "ln=" + ttos(ln) + "\n";
		SizeType counter = 0;
		while (tensors > 1) {
			createUlayer(tensors,counter);
			createWlayer(tensors,counter);
			tensors /= 2;
			counter++;
		}

		srep_ += "r(,)";
	}

	const PsimagLite::String& operator()() const
	{
		return srep_;
	}

private:

	void createUlayer(SizeType n, SizeType counter)
	{
		PsimagLite::String I0("");
		PsimagLite::String I1("");
		PsimagLite::String O0("");
		PsimagLite::String O1("");
		for (SizeType i = 0; i < n; ++i) {
			bool hasTwoOutputs = true;
			if (counter == 0) {
				I0 = "f" + ttos(2*i);
				I1 = "f" + ttos(2*i+1);
			}

			if (counter & 1) {
				if (i + 1 == n) hasTwoOutputs = false;
			} else {
				if (i == 0) hasTwoOutputs = false;
			}

			srep_ += "u("+I0+"," + I1 + "|" + O0;
			if (hasTwoOutputs) srep_ += "," + O1;
			srep_ += ")";
		}
	}

	void createWlayer(SizeType n, SizeType counter)
	{
		for (SizeType i = 0; i < n; ++i)
			srep_ += "w(,|)";
	}

	SizeType logBase2Strict(SizeType x) const
	{
		if (x == 0) return 0;
		SizeType counter = 0;
		SizeType ones = 0;
		while (x > 0) {
			if (x & 1) ones++;
			if (ones > 1) return 0;
			x >>= 1;
			counter++;
		}

		assert(counter > 0);
		return counter - 1;
	}

	PsimagLite::String srep_;
}; // class MeraBuilder
} // namespace Mera
#endif // MERABUILDER_H
