#!/usr/bin/perl6

use v6;

class Matrix {

	has $.rows;
	has $.cols;
	has $!data;

	method new($rows, $cols)
	{
		self.bless(:$rows, :$cols);
	}

	method addRow($r, $row)
	{
		loop (my $col = 0; $col < $.cols; ++$col) {
			$!data.[$col + $row*$.cols] = $r.[$col];
		}
	}

	method value($i, $j)
	{
		return $!data.[$i + $j*$.cols];
	}
}

sub MAIN($file, $lx) 
{
	my $input = open($file, :r);
	my $line = $input.get;
	my ($rows, $cols) = split(/\s/, $line);
	die "$*PROGRAM-NAME: Expected rows cols as first line of $file\n" unless
	($rows and $cols and $rows == $cols);

	my $row = 0;
	my $matrix = Matrix.new($rows, $cols);
	while ($line = $input.get) {
		my @temp = split(/\s/, $line);
		$matrix.addRow(@temp, $row++);
		last if ($row == $rows);
	}

	$input.close;

	loop (my $i = 0; $i < $rows; ++$i) {
		my $val = $matrix.value($i, $i);
		print "$val";
		print ($i > 0 && $i % $lx == 0) ?? "\n" :: " ";
	}
}




