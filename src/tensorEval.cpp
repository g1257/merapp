#include "TensorEval.h"
#include "Vector.h"

int main()
{
	PsimagLite::String str = "u0(s0)u1(s0)";

	SizeType dim0 = 5;
	typedef Mera::TensorEval<double> TensorEvalType;
	typedef TensorEvalType::TensorType TensorType;
	TensorEvalType::VectorTensorType vt(2,0);

	vt[0] = new TensorType(dim0);
	vt[1] = new TensorType(dim0);

	TensorEvalType tensorEval(str,vt);
	std::cout<<tensorEval.eval();

	for (SizeType i = 0; i < vt.size(); ++i) {
		delete vt[i];
		vt[i] = 0;
	}
}
