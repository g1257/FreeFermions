#!/usr/bin/perl

use strict;
use warnings;

my @drivers = ("cicj","deltaIdeltaJ","EasyExciton","HolonDoublon",
	       "reducedDensityMatrix","cicjBetaGrand",
	       "dynamics","hd2","ninj","niVsBetaGrand","splusSminus","szsz");

backupMakefile();
writeMakefile();
make();

sub make
{
	system("make");
}

sub backupMakefile
{
	system("cp Makefile Makefile.bak") if (-r "Makefile");
	print "Backup of Makefile in Makefile.bak\n";
}

sub writeMakefile
{
	my $allExecutables = combineAllDrivers("");
	my $allCpps = combineAllDrivers(".cpp");

open(FILE,">Makefile") or die "Cannot open Makefile for writing: $!\n";
print FILE<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl
# DMRG++ (v2.0) by G.A.
# Platform: Linux
# MPI: 0

LDFLAGS =      -llapack    -lm  -lpthread 
CPPFLAGS = -Werror -Wall -I../PartialPsimag  -I../../PsimagLite/src -I../src
CXX = g++ -O3 -DNDEBUG

all: $allExecutables
EOF

foreach my $what (@drivers) {
print FILE<<EOF;
$what.o: $what.cpp  Makefile
	\$(CXX) \$(CPPFLAGS) -c $what.cpp  

$what: $what.o
	\$(CXX) -o  $what $what.o \$(LDFLAGS)

EOF
}

print FILE<<EOF;
Makefile.dep: $allCpps
	\$(CXX) \$(CPPFLAGS) -MM $allCpps  > Makefile.dep

clean:
	rm -f core* $allExecutables *.o *.ii *.tt

include Makefile.dep
######## End of Makefile ########

EOF

close(FILE);
print "Done writing Makefile\n";
}

sub combineAllDrivers
{
	my ($extension) = @_;
	my $buffer = "";
	foreach my $what (@drivers) {
		my $tmp = $what.$extension." ";
		$buffer .= $tmp;
	}
	return $buffer;
}

