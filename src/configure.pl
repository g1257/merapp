#!/usr/bin/perl
=pod
Copyright (c) 2009-2018-2019, UT-Battelle, LLC
All rights reserved

[MERA++, Version 0.]

=cut
use warnings;
use strict;

use Getopt::Long qw(:config no_ignore_case);
use lib "../../PsimagLite/scripts";
use NewMake;
use lib ".";
use PsiTag;

my ($flavor, $lto) = (NewMake::noFlavor(), 0);
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

my $basicConfig = "../TestSuite/inputs/ConfigBase.psiTag";
my @configFiles = NewMake::configFilesList($basicConfig, $config);
my %driver1 = (name => 'meranpp');
my %driver2 = (name => 'merapp');
my %driver3 = (name => 'srepToTikz');
my %driver4 = (name => 'tensorBreakup');
my %driver5 = (name => 'tensorEval');


my @drivers = (\%driver1, \%driver2, \%driver3, \%driver4, \%driver5);


my $aInc = " -I/home/gonza/.exatn/include/exatn -I/home/gonza/.exatn/include -I/home/gonza/.exatn/include/cppmicroservices4 -std=gnu++11 -fPIC  -DPATH_MAX=4096 -Wno-attributes -DNO_GPU -DEXATN_SERVICE ";
my $aLibs = " -rdynamic -Wl,-rpath,/home/gonza/.exatn/lib -L /home/gonza/.exatn/lib -lCppMicroServices -ltalsh -lexatn -lexatn-numerics -lexatn-runtime -lexatn-runtime-graph -ldl -lpthread ";

my %args;
$args{"CPPFLAGS"} = $lto." $aInc";
$args{"LDFLAGS"} = $lto." $aLibs";
$args{"flavor"} = $flavor;
$args{"code"} = "DMRG++";
$args{"configFiles"} = \@configFiles;
$args{"additional3"} = "GitRevision.h";
$args{"additional4"} = $args{"additional3"};

#system("./createGitRevision.pl GitRevision.h");

createMakefile(\@drivers, \%args);

sub createMakefile
{
	my ($drivers, $args) = @_;
	unlink("Makefile.dep");
	NewMake::backupMakefile();

	my $fh;
	open($fh, ">", "Makefile") or die "Cannot open Makefile for writing: $!\n";

	NewMake::main($fh, $args, $drivers);
	local *FH = $fh;
print FH<<EOF;

.PHONY: GitRevision.h

GitRevision.h:
	./createGitRevision.pl GitRevision.h

EOF

	close($fh);
	print STDERR "$0: File Makefile has been written\n";
}
