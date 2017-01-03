#ifndef IRREDUCIBLEIDENTITY_H
#define IRREDUCIBLEIDENTITY_H

namespace Mera {

class IrreducibleIdentity {

public:

	IrreducibleIdentity()
	    : maxIndex_(0)
	{}

	void increase()
	{
		++maxIndex_;
	}

	SizeType maxIndex() const
	{
		return maxIndex_;
	}

	PsimagLite::String dimensionSrep() const
	{
		PsimagLite::String str("");
		for (SizeType i = 0; i < maxIndex_; ++i) {
			PsimagLite::String args = "D2,D2"; // FIXME
			str += "i" + ttos(i) + "(" + args + ")";
		}

		return str;
	}

private:

	SizeType maxIndex_;

}; // class IrreducibleIdentity

} // namespace Mera
#endif // IRREDUCIBLEIDENTITY_H
