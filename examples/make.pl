#!/usr/bin/perl

use strict;
use warnings;

use lib '../../PsimagLite/scripts';
use Make;

my @drivers = ("cicj","deltaIdeltaJ","EasyExciton","HolonDoublon","decay",
	       "reducedDensityMatrix","cicjBetaGrand",
	       "dynamics","hd2","ninj","niVsBetaGrand","splusSminus","szsz",
               "SzSzTime","etd","sqOmega","WavePacket");

Make::backupMakefile();
writeMakefile();
#make();

sub make
{
	system("make");
}

sub writeMakefile
{
	open(my $fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";

	my $cppflags = " -I../../PsimagLite ";
	$cppflags .= " -I../../PsimagLite/src -I../src";
	my $cxx="g++ -ansi -pedantic -O3 -DNDEBUG";
	my $lapack = Make::findLapack();

	Make::make($fh,\@drivers,"FreeFermions","Linux",0,
	"-lpsimaglite $lapack    -lm  -lpthread ",$cxx,$cppflags,"true"," "," ");

	close($fh);
	print "$0: Done writing Makefile\n";
}
