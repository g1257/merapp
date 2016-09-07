# DO NOT EDIT!!! Changes will be lost. Modify Config.make instead
# This Makefile was written by configure.pl
# MPS++ by G.A.

include Config.make
CPPFLAGS += -I../../PsimagLite -I../../PsimagLite/src -IEngine
all: merapp

merapp: merapp.o 
	$(CXX) merapp.o -o merapp $(LDFLAGS)

merapp.o: merapp.cpp  Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c merapp.cpp

../../PsimagLite/lib/libpsimaglite.a:
	$(MAKE) -f Makefile -C ../../PsimagLite/lib/

Makefile.dep: merapp.cpp   
	$(CXX) $(CPPFLAGS) -MM merapp.cpp   > Makefile.dep

clean: Makefile.dep
	rm -f core* merapp  *.o *.dep  

include Makefile.dep 
