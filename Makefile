# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl
# DMRG++ (v2.0) by G.A.
# Platform: Linux
# MPI: 0

LDFLAGS =      -llapack    -lm  -lpthread 
CPPFLAGS = -Werror -Wall -IPartialPsimag  -Isrc
CXX = g++ -g3 -DNDEBUG
all: $(EXENAME)
OBJECTSCOMMON = 
OBJECTSDMRG    = main.o  

all: clean ./freeFermions 

freeFermions: $(OBJECTSDMRG) $(OBJECTSCOMMON) 
	$(CXX) -o freeFermions  $(OBJECTSCOMMON) $(OBJECTSDMRG)  $(LDFLAGS) 

freeSystem: clean  $(OBJECTSCOMMON)
	$(CXX) -o  $(OBJECTSCOMMON) $(LDFLAGS)

clean:
	rm -f core* $(EXENAME) *.o *.ii *.tt

######## End of Makefile ########

