# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl
# DMRG++ (v2.0) by G.A.
# Platform: Linux
# MPI: 0

LDFLAGS =      -llapack    -lm  -lpthread 
CPPFLAGS = -Werror -Wall -I../PartialPsimag  -I../src
CXX = g++ -g3 -DNDEBUG


all: clean ninj 

ninj.o: ninj.cpp 
	$(CXX) $(CPPFLAGS) -c ninj.cpp  

ninj: clean  ninj.o
	$(CXX) -o  ninj ninj.o $(LDFLAGS)

clean:
	rm -f core* ninj *.o *.ii *.tt

######## End of Makefile ########

