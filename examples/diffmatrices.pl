#!/usr/bin/perl -w

use strict;
my $M_PI=3.1415927;

my @m;
my $n = loadMatrix("matrix.txt",\@m);
my @m2;
my $n2 = loadMatrix("matrix2.txt",\@m2);
my @mdiff;
diffMatrices(\@mdiff,\@m,\@m2,$n*$n);

printMatrix(\@mdiff,$n);

sub loadMatrix
{
	my ($file,$ma)=@_;
	my $n;
	my $j = 0;
	open(FILE,$file) or die "Cannot open file $file: $!\n";
	while(<FILE>) {
		chomp;
		next if ($_ eq "" || /^#/);
		my @temp=split;
		$n = $#temp+1 if (!defined($n));
		for (my $i=0;$i<$n;$i++) {
			$ma->[$i+$j*$n] = $temp[$i];
		}
		$j++;
	}
	close(FILE);
	return $n;
}



sub printMatrix
{
	my ($ma,$n)=@_;
	for (my $i=0;$i<$n;$i++) {
		for (my $j=0;$j<$n;$j++) {
			print "$ma->[$i+$j*$n] ";
		}
		print "\n";
	}
}

sub diffMatrices
{
	my ($ma,$m1a,$m2a,$n)=@_;
	for (my $i=0;$i<$n;$i++) {
		$ma->[$i] = $m1a->[$i] - $m2a->[$i];;
	}
}
