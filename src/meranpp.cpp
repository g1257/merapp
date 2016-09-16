#include "Vector.h"
#include "TensorOptimizer.h"

int main(int argc, char** argv)
{
	typedef Mera::TensorOptimizer<double> TensorOptimizerType;
	typedef TensorOptimizerType::IoInType IoInType;
	PsimagLite::String file = "meraEnviron.txt";
	if (argc >= 2) file = argv[1];

	IoInType io(file);
	SizeType tauMax = 0;
	io.readline(tauMax,"TauMax=");

	PsimagLite::String str("");
	io.readline(str,"Y_FOR_TENSOR=");
	PsimagLite::Vector<PsimagLite::String>::Type tokens;
	PsimagLite::tokenizer(str,tokens,",");
	if (tokens.size() != 2)
		throw PsimagLite::RuntimeError(PsimagLite::String(argv[0])
	        + ": Error reading Y_FOR_TENSOR\n");

	TensorOptimizerType to(io,tokens[0],std::atoi(tokens[1].c_str()));

	to.optimize(10);
}
