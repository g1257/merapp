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
		SizeType counter = 0;
		SizeType summed = 0;
		SizeType savedSummedForU = 0;
		SizeType idsU = 0;
		SizeType idsW = 0;
		while (tensors > 1) {
			SizeType savedSummedForW = summed;
			createUlayer(summed,idsU,savedSummedForU,tensors,counter);
			savedSummedForU = summed;
			createWlayer(summed,idsW,savedSummedForW,tensors,counter);
			tensors /= 2;
			counter++;
		}

		assert(summed > 1);
		srep_ += "r(s" + ttos(summed-2) + ",s" + ttos(summed-1) + ")";
	}

	const PsimagLite::String& operator()() const
	{
		return srep_;
	}

private:

	void createUlayer(SizeType& summed,
	                  SizeType& idsU,
	                  SizeType savedSummed,
	                  SizeType n,
	                  SizeType counter)
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
			} else {
				I0 = "s" + ttos(savedSummed++);
				I1 = "s" + ttos(savedSummed++);
			}

			if (counter & 1) {
				if (i + 1 == n) hasTwoOutputs = false;
			} else {
				if (i == 0) hasTwoOutputs = false;
			}

			O0 = "s" + ttos(summed++);
			srep_ += "u" + ttos(idsU++) + "("+I0+"," + I1 + "|" + O0;
			if (hasTwoOutputs) {
				O1 = "s" + ttos(summed++);
				srep_ += "," + O1;
			}

			srep_ += ")";
		}
	}

	void createWlayer(SizeType& summed,
	                  SizeType& idsW,
	                  SizeType savedSummed,
	                  SizeType n,
	                  SizeType counter)
	{
		PsimagLite::String I0("");
		PsimagLite::String I1("");
		PsimagLite::String O0("");
		for (SizeType i = 0; i < n; ++i) {
			bool hasTwoInputs = true;

			if (counter & 1) {
				if (i == 0) hasTwoInputs = false;
			} else {
				if (i + 1 == n) hasTwoInputs = false;
			}

			I0 = "s" + ttos(savedSummed++);
			srep_ += "w" + ttos(idsW++) +"(" + I0;
			if (hasTwoInputs) {
				I1 = "s" + ttos(savedSummed++);
				srep_ += "," + I1;
			}

			O0 = "s" + ttos(summed++);
			srep_ += "|" + O0 + ")";
		}
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
