#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Test::More tests=>4;
use FindBin qw/$RealBin/;
use lib "$RealBin/../lib";

use_ok 'Bio::Minimizer';

srand(42);
my @nt = qw(A T C G);
my $alphabetSize = scalar(@nt);
my $sequence = "";
for(1..240){
  $sequence .= $nt[int(rand($alphabetSize))]
}

my $minimizer = Bio::Minimizer->new($sequence,{numcpus=>2});

is($$minimizer{k}, 31, "Expected default kmer length");
is($$minimizer{l}, 21, "Expected default lmer length");

subtest 'Kmer => minimizer' => sub{
  plan tests=>6;
  note "Forward kmer";
  is($$minimizer{minimizers}{TCAGTCACAAGAGGCGCTCAGACCGACCTGC}, "AAGAGGCGCTCAGACCGACCT", "TCAGTCACAAGAGGCGCTCAGACCGACCTGC");
  is($$minimizer{minimizers}{TTGCTTCACCGCTACGCAGGCCTCTATTCCA}, "ACCGCTACGCAGGCCTCTATT", "TTGCTTCACCGCTACGCAGGCCTCTATTCCA");
  is($$minimizer{minimizers}{GTCCAGCGTCTTTGAGGGTAATCATTCGAGG}, "AGCGTCTTTGAGGGTAATCAT", "GTCCAGCGTCTTTGAGGGTAATCATTCGAGG");

  # revcom minimizers
  note "Revcom kmer";
  is($$minimizer{minimizers}{TGGGTCTGTCCCTCACGTACCGACTAAAGTA}, "CCCTCACGTACCGACTAAAGT", "TGGGTCTGTCCCTCACGTACCGACTAAAGTA");
  is($$minimizer{minimizers}{CCGACTAAAGTATTAAAAGTGATTCTGGGGC}, "AAAGTATTAAAAGTGATTCTG", "CCGACTAAAGTATTAAAAGTGATTCTGGGGC");
  is($$minimizer{minimizers}{CCGACTAAAGTATTAAAAGTGATTCTGGGGC}, "AAAGTATTAAAAGTGATTCTG", "CCGACTAAAGTATTAAAAGTGATTCTGGGGC");
};

