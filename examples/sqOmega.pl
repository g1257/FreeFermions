#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use OmegaUtils;
use Getopt::Long qw(:config no_ignore_case);

my $usage = "-f input -o omegaBegin -i omegaStep -t omegaTotal ";
$usage .= " [-w observable] [-M mMax] [-p] [-r]\n";

my ($input1,$geometry,$GlobalNumberOfSites);
my ($isPeriodic,$mMax,$wantsRealPart);
my ($omega0,$omegaStep,$centralSite,$total);
my $observable = "sz";

GetOptions('f=s' => \$input1,
           'p' => \$isPeriodic,
           'M:i' => \$mMax,
           'r' => \$wantsRealPart,
           't:i' => \$total,
           'i:f' => \$omegaStep,
           'o:f' => \$omega0,
		   'w:s' => \$observable) or die "$usage\n";

(defined($input1) && defined($total) && defined($omegaStep)) or die "$0: USAGE: $usage\n";
defined($isPeriodic) or $isPeriodic = 0;
defined($omega0) or $omega0 = 0;

my $input = $input1;
$input =~ s/\.(.*$)//;
$input .= ".dat";
my $cmd = "./sqOmega -f $input1 -t $total -i $omegaStep -o $omega0 -w $observable";
$cmd .= " > $input 2> /dev/null";
print STDERR "$0: Trying to exec $cmd\n";
my $ret = system($cmd);

die "$0: Command $cmd failed\n" if ($ret != 0);

my $hptr = {"TSPSites 1" => \$centralSite,
            "GeometryKind" =>\$geometry,
            "TotalNumberOfSites" => \$GlobalNumberOfSites};

OmegaUtils::getLabels($hptr,$input);

my @spaceValues;
readSpace(\@spaceValues,$input);

my @qvalues;
my @omegas;
for (my $i = 0; $i < $total; ++$i) {
	$omegas[$i] = $omega0 + $omegaStep*$i;
	$qvalues[$i] = procThisOmega($omegas[$i],$spaceValues[$i]);
}

printGnuplot(\@qvalues,\@omegas,$geometry);

sub procThisOmega
{
	my ($omega,$spaceForThisOmega) = @_;

	my @qValues;
	fourier(\@qValues,$spaceForThisOmega,$geometry);
	return \@qValues;
}

sub readSpace
{
	my ($space,$inFile) = @_;
	my $counter = 0;

	open(FIN,"$inFile") or die "$0: Cannot open $inFile : $!\n";
	while(<FIN>) {
		if (/^#/) {
		        next;
		}

		my @temp=split;
		my $n = scalar(@temp);
		next unless ($n == 2*$GlobalNumberOfSites + 1);
		my @temp2;
		for (my $i = 1; $i < $n; ++$i) {
			$temp2[$i-1] = $temp[$i];
		}

		$space->[$counter++] = \@temp2;
	}

	close(FIN);
	print STDERR "$0: Read $counter omegas\n";
}

sub printGnuplot
{
	my ($array,$omegas,$geometry) = @_;

	my $numberOfOmegas = scalar(@$omegas);

	my $factor = 0;
	my @fileIndices=(0);
	if ($geometry eq "chain") {
		$factor = 1.0;
	} elsif ($geometry eq "ladder") {
		$factor = 0.5;
		@fileIndices=(0,1);
	} else {
		die "$0: Unknown geometry $geometry\n";
	}

	foreach my $fileIndex (@fileIndices) {
		my $outFile = "outSpectrum$fileIndex.gnuplot";
		open(FOUT,"> $outFile") or die "$0: Cannot write to $outFile : $!\n";

		for (my $i = 0; $i < $numberOfOmegas; ++$i) {
			my $omega = $omegas->[$i];
			my $a = $array->[$i];
			my $numberOfQs = int($factor*scalar(@$a));
			for (my $m = 0; $m < $numberOfQs; ++$m) {
				my $q = getQ($m,$numberOfQs);
				my $realAndImag = $a->[$m + $fileIndex*$numberOfQs];
				my $n2 = scalar(@$realAndImag);
				$n2 == 2 or die "$0: Error $n2 != 2\n";
				my $realPart = $realAndImag->[0];
				my $imagPart = $realAndImag->[1];
				print FOUT "$q $omega $realPart $imagPart\n";
			}
		}

		close(FOUT);
		print "$0: Written $outFile\n";
	}
}
sub fourier
{
	my ($f,$v,$geometry) = @_;

	if ($geometry eq "chain") {
		return fourierChain($f,$v);
	}

	if ($geometry eq "ladder") {
		return fourierLadder($f,$v);
	}

	die "$0: ft: undefined geometry $geometry\n";
}

sub fourierChain
{
	my ($f,$v) = @_;
	my $n = int(0.5*scalar(@$v));
	my $numberOfQs = (defined($mMax)) ? $mMax : $n;
	for (my $m = 0; $m < $numberOfQs; ++$m) {
		my @sum = (0,0);
		my $q = getQ($m,$numberOfQs);
		for (my $i = 0; $i < $n; $i++) {
			my @temp = ($v->[2*$i],$v->[2*$i+1]);
			my $arg = $q*($i-$centralSite);
			my $carg = cos($arg);
			$sum[0] += $temp[0]*$carg;
			$sum[1] += $temp[1]*$carg;
		}

		$f->[$m] = \@sum;
	}
}

sub fourierLadder
{
	my ($f,$v) = @_;
	my $n = int(0.5*scalar(@$v));
	my $numberOfQs = (defined($mMax)) ? $mMax : int(0.5*$n);
	for (my $m = 0; $m < $numberOfQs; ++$m) {
		my $q = getQ($m,$numberOfQs);
		my @f0 = fourierF0($v,$q);
		my @f1 = fourierF1($v,$q);
		for (my $x = 0; $x < 2; ++$x) {
			my $sign = 1-2*$x;
			my $realPart = $f0[0] + $sign*$f1[0];
			my $imagPart = $f0[1] + $sign*$f1[1];
			my @sum = ($realPart,$imagPart);
			$f->[$m+$numberOfQs*$x] = \@sum;
		}
	}
}

sub fourierF0
{
	my ($v,$q) = @_;
	my $n = int(0.5*scalar(@$v));
	my @sum;
	for (my $i = 0; $i < $n; $i+=2) {
		my @temp = ($v->[2*$i],$v->[2*$i+1]);
		my $arg = $q*($i-$centralSite)*0.5;
		my $carg = cos($arg);
		$sum[0] += $temp[0]*$carg;
		$sum[1] += $temp[1]*$carg;
	}

	return @sum;
}

sub fourierF1
{
	my ($v,$q) = @_;
	my $n = int(0.5*scalar(@$v));
	my @sum;
	for (my $i = 1; $i < $n; $i+=2) {
		my @temp = ($v->[2*$i],$v->[2*$i+1]);
		my $arg = $q*distanceLadder($i,$centralSite);
		my $carg = cos($arg);
		$sum[0] += $temp[0]*$carg;
		$sum[1] += $temp[1]*$carg;
	}

	return @sum;
}

sub distanceLadder
{
	my ($ind, $jnd) = @_;
	my $first = ($ind-1)/2;
	my $second = $jnd/2;
	return $first - $second;
}

sub getQ
{
	my ($m,$n) = @_;
	return ($isPeriodic) ? 2.0*pi*$m/$n : pi*$m/($n+1.0);
}



