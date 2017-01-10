#!/usr/bin/perl
#
use strict;
use warnings;
use Math::MatrixReal;

my $pixelPerElement = 8;

my ($label) = @ARGV;
defined($label) or die "USAGE: cat matrix.txt | $0 label\n";

while (<STDIN>) {
	last if (/^$label/);
}

$_ = <STDIN>;
chomp;
my @temp = split;
if (scalar(@temp) != 2) {
	die "$0: Line $. expected rows cols\n";
}

my ($rows, $cols) = @temp;

print STDERR "$0: Matrix $rows times $cols\n";

my @m;
my $counter = 0;
while (<STDIN>) {
	chomp;
	my $row = $_;
	loadThisRow(\@m,$counter,$row,$cols);
	$counter++;
	last if ($counter >= $rows);
}

if ($counter != $rows) {
	die "$0: Error: Too few rows, expecting $rows found $counter\n";
}

my $matrix = Math::MatrixReal->new_from_rows(\@m);

my %usedColors;
my @colors = ("0 0 0","255 0 0", "0 255 0", "0 0 255","255 255 0", "255 0 255", "0 255 255",
"255 128 0", "128 255 0", "128 0 255", "255 255 128", "255 128 255", "128 255 255");
printArtistic($matrix,$pixelPerElement, \%usedColors, \@colors);

hermiticity($matrix);

isBlockDiagonal($matrix);

sub loadThisRow
{
	my ($m,$row,$data,$cols) = @_;
	my @temp = split(/ /,$data);
	my $n = scalar(@temp);
	if ($cols != $n) {
		die "$0: Row $row has $n columns, $cols expected\n";
	}

	$m->[$counter] = \@temp;
}

sub printArtistic
{
	my ($m, $ppe, $usedColors, $colors) = @_;
	my ($rows, $cols) = $m->dim();
	my ($drows, $dcols) = ($ppe*$rows, $ppe*$cols);
	print "P3\n$drows $dcols\n255\n";
	for (my $i = 0; $i < $rows; ++$i) {
		# compute row
		my $thisRow = "";
		for (my $j = 0; $j < $cols; ++$j) {
			my $e = $m->element($i+1,$j+1);
			my $val = getArtisticRow($e, $ppe, $usedColors, $colors);
			$thisRow .= "$val ";
		}

		#print row
		for (my $j = 0; $j < $ppe; ++$j) {
			print "$thisRow\n";
		}
	}
}

sub getArtisticRow
{
	my ($e, $ppe, $usedColors, $colors) = @_;
	my $val = getArtisticValue($e, $usedColors, $colors);
	my $result = "";
	for (my $i = 0; $i < $ppe; ++$i) {
		$result .= "$val ";
	}

	return $result;
}

sub getArtisticValue
{
	my ($e,$usedColors,$colors) = @_;
	return "255 255 255" if ($e == 0);
	my $max = keys %{$usedColors};

	my $color = $usedColors->{"$e"};
	if (defined($color)) {
		$color =~ s/,/ /g;
		return $color;
	}

	die "$0: We run out of colors\n" if ($max >= scalar(@$colors));
	$color = $colors->[$max];
	print STDERR "New color $color for $e\n";
	my $colorWithCommas = $color;
	$colorWithCommas =~ s/ /,/g;
	$usedColors->{"$e"} = $colorWithCommas;
	return $color;
}

sub hermiticity
{
	my ($m) = @_;
	my ($rows, $cols) = $m->dim();
	for (my $i = 0; $i < $rows; ++$i) {
		for (my $j = $i + 1; $j < $cols; ++$j) {
			my $e = $m->element($i+1,$j+1);
			my $et = $m->element($j+1,$i+1);
			next if ($e == $et);
			print STDERR "hermiticity: $i $j     $e != $et\n";
		}
	}
}

sub isBlockDiagonal
{
	my ($m) = @_;
	my ($rows, $cols) = $m->dim();
	my $offset = 0;
	my $counter = 0;
	while ($offset < $rows) {
		my $c = findBlockStartingAt($m,$offset);
		my $s = $c + 1 - $offset;
		print STDERR "$0: Found block $counter starting at $offset of size $s\n";
		$offset = $c + 1;
		$counter++;
	}
}

sub findBlockStartingAt
{
	my ($m, $start) = @_;
	my ($rows, $cols) = $m->dim();
	my $c = $start;
	while ($c < $rows) {
		my $tmp = findLastNzColumn($m,$c);
		last if ($tmp <= $c);
		$c = $tmp + 1;
	}

	return $c;
}

sub findLastNzColumn
{
	my ($m, $row) = @_;
	my ($rows, $cols) = $m->dim();
	my $col = $row;
	for (my $j = $row; $j < $cols; ++$j) {
		my $e = $m->element($row+1,$j+1);
		$col = $j if ($e != 0);
	}

	return $col;
}

