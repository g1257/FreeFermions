#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use OmegaUtils;
use Getopt::Long qw(:config no_ignore_case);

my $usage = "-f input -o omegaBegin -i omegaStep -t omegaTotal ";
$usage .= " [-w observable] [-M mMax] [-p] [-r]\n";

my ($input1,$orbitals,$GlobalNumberOfSites);
my ($isPeriodic,$mMax,$wantsRealPart);
my ($omega0,$omegaStep,$centralSite,$total);
my $observable = "sz";
my $zeroAtCenter = 0;

GetOptions('f=s' => \$input1,
           'p' => \$isPeriodic,
           'M:i' => \$mMax,
           'r' => \$wantsRealPart,
           'z' => \$zeroAtCenter,
           't:i' => \$total,
           'i:f' => \$omegaStep,
           'o:f' => \$omega0,
	   'w:s' => \$observable) or die "$usage\n";

(defined($input1) && defined($total) && defined($omegaStep)) or die "$0: USAGE: $usage\n";
defined($isPeriodic) or $isPeriodic = 0;
defined($omega0) or $omega0 = 0;

my $geometryName;
my $geometryLeg = 1;
my $subgeometry = "UNDEFINED";
my $hptr = {"TSPSites 1" => \$centralSite,
            "GeometrySubKind" =>\$subgeometry,
            "GeometryKind" => \$geometryName,
            "LadderLeg" => \$geometryLeg,
            "Orbitals" =>\$orbitals,
            "TotalNumberOfSites" => \$GlobalNumberOfSites};

OmegaUtils::getLabels($hptr,$input1);
$hptr->{"isPeriodic"} = $isPeriodic;
$hptr->{"mMax"} = $mMax;
$hptr->{"centralSite"} = $centralSite;

my $geometry = {"name" => $geometryName, "leg" => $geometryLeg, "subname" => $subgeometry};

my $input = $input1;
$input =~ s/\.(.*$)//;
$input .= ".dat";

#my $cmd = "./sqOmega -f $input1 -t $total -i $omegaStep -o $omega0 -c $centralSite -w $observable";
#$cmd .= " > $input 2> /dev/null";
#print STDERR "$0: Trying to exec $cmd\n";
#my $ret = system($cmd);
#die "$0: Command $cmd failed\n" if ($ret != 0);

my @spaceValues;
readSpace(\@spaceValues, $input);

my @freqK;
my %h;
for (my $i = 0; $i < $total; ++$i) {
	my $thisOmega = $spaceValues[$i];
	my $omega = shift @$thisOmega;
	defined($omega) or die "$0: Undefined for $i\n";
	my @qValues;
	OmegaUtils::fourier(\@qValues, $thisOmega, $geometry, $hptr);
	my $nf = scalar(@qValues); # n. of m values
	my @final = ($omega);
	print STDERR "$0: Found $nf m values\n";
	for (my $j = 0; $j < $nf; ++$j) {
		my $temp = $qValues[$j];
		my $ntemp = scalar(@$temp); # == 2
		die "$0: Expected array of 2 values (real, imag) but found $ntemp instead\n" if ($ntemp != 2);
		push @final, @$temp;
	}
		
	$h{$omega} = \@final;
}

OmegaUtils::printGnuplot(\%h, $geometry, $isPeriodic, $zeroAtCenter);

sub readSpace
{
	my ($space,$inFile) = @_;
	my $counter = 0;

	open(FIN, "<", "$inFile") or die "$0: Cannot open $inFile : $!\n";
	while(<FIN>) {
		if (/^#/) {
		        next;
		}

		my @temp=split;
		my $n = scalar(@temp);
		next unless ($n == 2*$GlobalNumberOfSites + 1);
		my @rAi;
		my @oneOmega = ($temp[0]);
		for (my $i = 0; $i < $GlobalNumberOfSites; ++$i) {
			my $offset = 2*$i + 1;
			my @rAi = ($temp[$offset], $temp[$offset + 1]);
			$oneOmega[$i + 1] = \@rAi;
		}

		$space->[$counter++] = \@oneOmega;
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
        } elsif ($geometry eq "LongRange") {

                defined($subgeometry) or die "$0 LongeRange geometry: need to specify geometry (chain or ladder) for the Fourier Transform\n";

		if ($subgeometry eq "chain") {
			$factor = 1.0;
		} elsif ($subgeometry eq "ladder") {
			$factor = 0.5;
			@fileIndices=(0,1);
		}
	} else {
		die "$0: Unknown geometry $geometry\n";
	}

	foreach my $fileIndex (@fileIndices) {
		my $outFile = "outSpectrum$fileIndex.gnuplot";
		open(FOUT, ">", "$outFile") or die "$0: Cannot write to $outFile : $!\n";

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






