#!/usr/bin/perl -w

my ($what)=@ARGV;
my $possibilities = getPossibilities();
defined($what) or die "perl make.pl whatever\n \twhere whatever is one of $possibilities\n";
$what=~s/\.cpp$//;

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
open(FILE,">Makefile") or die "Cannot open Makefile for writing: $!\n";
my $headers = join(' ',glob("../Engine/*.h"));
print FILE<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl
# DMRG++ (v2.0) by G.A.
# Platform: Linux
# MPI: 0

LDFLAGS =      -llapack    -lm  -lpthread 
CPPFLAGS = -Werror -Wall -I../PartialPsimag  -I../src
CXX = g++ -O2 -pg -DNDEBUG


all: clean $what 

$what.o: $what.cpp $headers Makefile
	\$(CXX) \$(CPPFLAGS) -c $what.cpp  

$what: clean  $what.o
	\$(CXX) -o  $what $what.o \$(LDFLAGS)

clean:
	rm -f core* $what *.o *.ii *.tt

######## End of Makefile ########

EOF

close(FILE);
print "Done writing Makefile\n";
}

sub getPossibilities
{
	open(PIPE,"ls *.cpp |") or return "";
	my $a = "";
	while(<PIPE>) {
		chomp;
		s/\.cpp$//;
		$a .= $_." ";
	}
	close(PIPE);
	return $a;
}

