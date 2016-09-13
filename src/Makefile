# DO NOT EDIT!!! Changes will be lost. Modify Config.make instead
# This Makefile was written by configure.pl
# MPS++ by G.A.

include Config.make
CPPFLAGS += -I../../PsimagLite -I../../PsimagLite/src -IEngine
all: merapp srepToTikz tensorEval

merapp: merapp.o 
	$(CXX) merapp.o -o merapp $(LDFLAGS)

merapp.o: merapp.cpp  Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c merapp.cpp

srepToTikz: srepToTikz.o
	$(CXX) srepToTikz.o -o srepToTikz $(LDFLAGS)

srepToTikz.o: srepToTikz.cpp  Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c srepToTikz.cpp

tensorEval: tensorEval.o
	$(CXX) tensorEval.o -o tensorEval $(LDFLAGS)

tensorEval.o: tensorEval.cpp Makefile   Config.make
	$(CXX) $(CPPFLAGS) -c  tensorEval.cpp

../../PsimagLite/lib/libpsimaglite.a:
	$(MAKE) -f Makefile -C ../../PsimagLite/lib/

Makefile.dep: merapp.cpp   
	$(CXX) $(CPPFLAGS) -MM merapp.cpp  srepToTikz.cpp tensorEval.cpp > Makefile.dep

clean: Makefile.dep
	rm -f core* merapp  *.o *.dep  

include Makefile.dep 
