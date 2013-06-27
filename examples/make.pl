#!/usr/bin/perl

use strict;
use warnings;

use lib '../../PsimagLite/scripts';
use Make;

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
	open(my $fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";

	my $cppflags = "-Werror -Wall -I../../PsimagLite  -I../../PsimagLite/src -I../src";
	my $cxx="g++ -O3 -DNDEBUG";
	my $lapack = Make::findLapack();

	Make::make($fh,\@drivers,"FreeFermions","Linux",0,"$lapack    -lm  -lpthread",$cxx,$cppflags,"true"," "," ");

	close($fh);
	print "Done writing Makefile\n";
}
