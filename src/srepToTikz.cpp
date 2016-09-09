#include "MeraToTikz.h"
#include <iostream>
#include "Vector.h"
#include <cstdlib>
#include "Version.h"

int main(int argc, char** argv)
{
	if (argc == 1) {
		std::cerr<<"USAGE: "<<argv[0]<<" tauMax\n";
		return 1;
	}

	char c = '\0';
	PsimagLite::String srep;
	while (std::cin>>c) {
		srep += c;
	}

	Mera::MeraToTikz<double> obj(srep,atoi(argv[1]));
	std::cout<<"%Created by "<<argv[0]<<" version "<<MERA_VERSION<<"\n";
	std::cout<<obj;
}
