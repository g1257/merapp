#include "Vector.h"
#include "IoSimple.h"

int main(int argc, char** argv)
{
	PsimagLite::String file = "meraEnviron.txt";
	if (argc >= 2) file = argv[1];

	PsimagLite::IoSimple::In io(file);
	SizeType tauMax = 0;
	io.readline(tauMax,"TauMax=");

}
