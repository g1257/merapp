#include "MeraToTikz.h"
#include <iostream>
#include "Vector.h"
#include <cstdlib>

int main(int argc, char** argv)
{
	char c = '\0';
	PsimagLite::String srep;
	while (std::cin>>c) {
		srep += c;
	}

	if (argc == 1) {
		std::cerr<<"USAGE: "<<argv[0]<<" tauMax\n";
		return 1;
	}

	Mera::MeraToTikz<double> obj(srep,atoi(argv[1]));
	std::cout<<obj;
}
