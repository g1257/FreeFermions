#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use lib '../../PsimagLite/scripts';
use NewMake;
use PsiTag;
my ($flavor, $lto) = (NewMake::noFlavor() , 0);
my $usage = "USAGE: $0 [-f flavor] [-lto] [-c config]\n";
my $config;

GetOptions('f=s' => \$flavor,
           'lto' => \$lto,
           'c=s' => \$config) or die "$usage\n";

my $gccdash = "";
if ($lto == 1) {
	$gccdash = "gcc-";
	$lto = "-flto";
} else {
	$lto = "";
}

my @configFiles = ("../TestSuite/inputs/ConfigBase.psiTag");
push @configFiles, $config if (defined($config));

my %args;
$args{"CPPFLAGS"} = $lto;
$args{"LDFLAGS"} = $lto;
$args{"flavor"} = $flavor;
$args{"code"} = "FreeFermions";
$args{"configFiles"} = \@configFiles;

my @drivers = ("cicj","deltaIdeltaJ","EasyExciton","HolonDoublon","decay",
	       "reducedDensityMatrix","cicjBetaGrand",
	       "dynamics","hd2","ninj","niVsBetaGrand","splusSminus","szsz",
               "SzSzTime","etd","sqOmega","WavePacket", "WavePacket2");

createMakefile(\@drivers, \%args);

sub createMakefile
{
	my ($drivers, $args) = @_;
	NewMake::backupMakefile();

	my $fh;
	open($fh, ">", "Makefile") or die "Cannot open Makefile for writing: $!\n";

	NewMake::main($fh, $args, \@drivers);

	close($fh);
	print STDERR "$0: File Makefile has been written\n";
}

