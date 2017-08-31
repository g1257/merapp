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
defined($_) or die "$0: Could not find $label in STDIN\n";

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

findConnected($matrix);

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

	# print summary of values
	foreach my $key (keys %$usedColors) {
		my $ptr = $usedColors{"$key"};
		my $counter = $ptr->[0];
		print STDERR "Values $key appeared $counter times\n";
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

	my $ptr = $usedColors->{"$e"};
	if (defined($ptr)) {
		my ($counter, $color) = @$ptr;
		$ptr->[0] = $counter + 1;
		$color =~ s/,/ /g;
		return $color;
	}

	die "$0: We run out of colors\n" if ($max >= scalar(@$colors));
	my $color = $colors->[$max];
	print STDERR "New color $color for $e\n";
	my $colorWithCommas = $color;
	$colorWithCommas =~ s/ /,/g;
	$usedColors->{"$e"} = [1, $colorWithCommas];
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

sub depthFirstSearch
{
    my ($visited, $m, $ind) = @_;
    my $cols = scalar(@$visited);
    for (my $j = 0; $j < $cols; ++$j) {
        my $e = $m->element($ind+1,$j+1);
        next if ($e == 0);
        next if ($visited->[$j]);
        $visited->[$j] = 1;
        depthFirstSearch($visited, $m, $j);
    }
}

sub findConnected
{
    my ($m) = @_;
    my ($rows, $cols) = $m->dim();
    my @visited;

    for (my $i = 0; $i < $rows; ++$i) {
        $visited[$i] = 0;
    }

    my @permutation;
    my $counter = 0;
    my $x = 0;
    for (my $i = 0; $i < $rows; ++$i) {
        next if ($visited[$i]);
        my @copy = @visited;
        $visited[$i] = 1;
        depthFirstSearch(\@visited, $m, $i);
        $x = evalDiffs(\@copy, \@visited, $x, \@permutation);
        ++$counter;
    }

    print STDERR "$0: Found $counter connected sub-graphs\n";
	$counter = scalar(@permutation);
	my $sum = 0;
	open(FOUT, ">", "permutation.txt") or die "$0: writing to permutation.txt : $!\n";
	print FOUT "PERMUTATION\n";
	print FOUT "$counter $counter\n";
	for (my $i = 0; $i < $counter; ++$i) {
		$sum += $permutation[$i];
		for (my $j = 0; $j < $counter; ++$j) {
			my $val = ($permutation[$i] == $j) ? 1.0 : 0.0;
			print FOUT "$val ";
		}

		print FOUT "\n";
	}

	close(FOUT);

	my $expectedSum = $counter*($counter - 1)/2;
	($sum == $expectedSum) or die "$0: Permutation check failed\n";


}

sub evalDiffs
{
    my ($a,$b, $x, $permutation) = @_;
    my $n = scalar(@$a);
    ($n == scalar(@$b)) or die "$0: printDiffs: vectors of different sizes\n";
    #print STDERR "$0: Subgraph: ";
    my $c = 0;
    for (my $i = 0; $i < $n; ++$i) {
        next if ($a->[$i] == $b->[$i]);
	#print STDERR "$i ";
	$permutation->[$x++] = $i;
        ++$c;
    }

    #print STDERR " [$c]\n";
    return $x;
}



