#!/usr/bin/perl -w

use strict;
my $M_PI=3.1415927;

my ($option)=@ARGV; # what order the matrix is in?
my @m;
my $n = loadMatrix(\@m);
my $GlobalDim = 2; # square lattice
checkMatrix2(\@m,$n);

my @ms;
my $sites=sumOverFlavors(\@ms,\@m,$n,$option);

my @v;
sumOverDirs(\@v,\@ms,$sites);
checkMatrix(\@ms,$sites);
printMatrix(\@ms,$sites);
printVector(\@v,$sites);
fourierTransform(\@v,$sites);

sub sumOverDirs
{
	my ($va,$msa,$n)=@_;
	for (my $x=0;$x<$n;$x++) {
		$va->[$x] = calcDirSum($msa,$x,$n)/$n;
	}
}

sub calcDirSum
{
	my ($msa,$x,$n)=@_;
	my $sum = 0;
	for (my $i=0;$i<$n;$i++) {
		my $j = sumDirection($i,$x,$n); # j = i+x
		$sum += $msa->[$i+$j*$n];
	}
	return $sum;
}

sub sumDirection
{
	my ($i,$x,$n) = @_;
	my $length = int(sqrt($n));
	my @v1 = vectorialDirection($i,$length);
	my @v2 = vectorialDirection($x,$length);
	
	my @v;
	vectorSum(\@v,\@v1,\@v2,$length);
	print STDERR "@v1 <--> @v\n" if ($x==4);
	return vectorialIndex(\@v,$length);
}

sub vectorSum
{
	my ($va,$v1,$v2,$length)=@_;
	for (my $i=0;$i<$GlobalDim;$i++) {
		$va->[$i] = $v1->[$i] + $v2->[$i];
		$va->[$i] -= $length if ($va->[$i]>=$length);
	}
}

sub vectorialDirection
{
	my ($ind,$length)=@_;
	my $x = $ind % $length;
	my $y = int($ind/$length);
	return ($x,$y);
}

sub vectorialIndex
{
	my ($va,$length)=@_;
	return $va->[0]+$va->[1]*$length;
}

sub sumOverFlavors
{
	my ($msa,$ma,$n,$option)=@_;
	my $sites = int($n/2);
	for (my $i=0;$i<$sites*$sites;$i++) {
		my $r1 = $i % $sites;
		my $r2 = int($i / $sites);
		$msa->[$i] = 0;
		for (my $orb1=0;$orb1<2;$orb1++) {
			my $j1 = $orb1 + $r1 * 2;
			$j1 = $r1 + $orb1*$sites if ($option);
			for (my $orb2=0;$orb2<2;$orb2++) {
				my $j2 = $orb2 + $r2*2;
				$j2 = $r2 + $orb2*$sites if ($option);
				my $j = $j1 + $j2*$n; 
				defined ($ma->[$j]) or die "Undefined for $j $j1 $j2 $sites $n\n";
				$msa->[$i] += $ma->[$j];
			}
		}
	}
	return $sites;
}

sub loadMatrix
{
	my ($ma)=@_;
	my $n;
	my $j = 0;
	while(<STDIN>) {
		chomp;
		next if ($_ eq "" || /^#/);
		my @temp=split;
		$n = $#temp+1 if (!defined($n));
		for (my $i=0;$i<$n;$i++) {
			$ma->[$i+$j*$n] = $temp[$i];
		}
		$j++;
	}
	return $n;
}

sub printVector
{
	my ($va,$n)=@_;
	for (my $i=0;$i<$n;$i++) {
		$_ = $va->[$i];
		print "$i $_\n";
	}
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

sub checkMatrix
{
	my ($ma,$n)=@_;
	my $length = int(sqrt($n));
	for (my $i=0;$i<$n;$i++) {
		my $ir = rotate($i,$length);
		my $j = $i;
		#for (my $j=0;$j<$n;$j++) {
			my $jr = rotate($j,$length);
			my $a = $ma->[$i+$j*$n];
			my $b = $ma->[$ir+$jr*$n];
			die "Different for $i $j $a $b\n"
				if ($a!=$b);
		#}
	}
}

sub checkMatrix2
{
	my ($ma,$n)=@_;
	my $sites = $n/2;
	my $length = int(sqrt($sites));
	for (my $i=0;$i<$n;$i++) {
		my $ri = $i%$sites;
		my $orb = int($i/$sites);
		my $ir = rotate($ri,$length)+$orb*$sites;
		my $j = $i;
		#for (my $j=0;$j<$n;$j++) {
			my $rj = $j%$sites;
			my $orb2 = int($j/$sites);
			next if ($orb!=$orb2);
			my $jr = rotate($rj,$length)+$orb2*$sites;
			my $a = $ma->[$i+$j*$n];
			my $b = $ma->[$ir+$jr*$n];
			die "check2: Different for $i $j $ir $jr $a $b\n"
				if ($a!=$b);
		#}
	}
}

sub rotate
{
	my ($i,$length)=@_;
	my @v = vectorialDirection($i,$length);
	return $v[1] + $v[0]*$length;
}

			
sub fourierTransform
{
	my ($va,$sites)=@_;
	my $length = int(sqrt($sites));
	my $sum = 0;
	for (my $k=0;$k<$sites;$k++) {
		my ($tmpr,$tmpi)=(0,0);
		my $kx = $k%$length;
		my $ky = int($k/$length);
			
		for (my $i=0;$i<$sites;$i++) {
			my $x = $i%$length;
			my $y = int($i/$length);
			my $tmpe = $x*$kx + $y*$ky;
			$tmpe *= 2.0*$M_PI/$length;
			$tmpr += $va->[$i]*cos($tmpe);
			$tmpi += $va->[$i]*sin($tmpe);
			$sum += $va->[$i];
		}
		print "$k $tmpr $tmpi\n";
	}
	$sum /= $sites;
	print "#sum of all cx = $sum\n";
}
