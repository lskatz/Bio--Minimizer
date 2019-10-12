#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Test::More tests=>5;
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

my $minimizer = Bio::Minimizer->new($sequence);

is($$minimizer{k}, 31, "Expected default kmer length");
is($$minimizer{l}, 21, "Expected default lmer length");

subtest 'Kmer => minimizer' => sub{
  plan tests=>6;
  is($$minimizer{minimizers}{TCAGTCACAAGAGGCGCTCAGACCGACCTGC}, "AAGAGGCGCTCAGACCGACCT", "forward kmer");
  is($$minimizer{minimizers}{TTGCTTCACCGCTACGCAGGCCTCTATTCCA}, "ACCGCTACGCAGGCCTCTATT", "forward kmer");
  is($$minimizer{minimizers}{GTCCAGCGTCTTTGAGGGTAATCATTCGAGG}, "AGCGTCTTTGAGGGTAATCAT", "forward kmer");

  # revcom minimizers
  is($$minimizer{minimizers}{TGGGTCTGTCCCTCACGTACCGACTAAAGTA}, "CCCTCACGTACCGACTAAAGT", "revcom kmer");
  is($$minimizer{minimizers}{CCGACTAAAGTATTAAAAGTGATTCTGGGGC}, "AAAGTATTAAAAGTGATTCTG", "revcom kmer");
  is($$minimizer{minimizers}{CCGACTAAAGTATTAAAAGTGATTCTGGGGC}, "AAAGTATTAAAAGTGATTCTG", "revcom kmer");
};

subtest 'Minimizer => kmer' => sub{
  plan tests=>12;

  # Six different minimizers at different array elements, leading to different kmers
  is($$minimizer{kmers}{TTATGAAGGGTTGCTTCACCG}[0], "GGGCTGATTGTTATGAAGGGTTGCTTCACCG", "forward kmer");
  is($$minimizer{kmers}{AAGAGGCGCTCAGACCGACCT}[1], "TTCAGTCACAAGAGGCGCTCAGACCGACCTG", "forward kmer");
  is($$minimizer{kmers}{AGGAACCGGACCTTTAATCAC}[2], "ATCATTCGAGGAACCGGACCTTTAATCACGG", "forward kmer");
  is($$minimizer{kmers}{ACCGCTACGCAGGCCTCTATT}[3], "TTGCTTCACCGCTACGCAGGCCTCTATTCCA", "forward kmer");
  is($$minimizer{kmers}{CCAGAATCACTTTTAATACTT}[5], "GGGCCCCAGAATCACTTTTAATACTTTAGTC", "forward kmer");
  is($$minimizer{kmers}{AGCGTCTTTGAGGGTAATCAT}[6], "GTCCAGCGTCTTTGAGGGTAATCATTCGAGG", "forward kmer");

  # revcom minimizers
  is($$minimizer{kmers}{CGGTTCCTCGAATGATTACCC}[0], "ATTAAAGGTCCGGTTCCTCGAATGATTACCC", "revcom kmer");
  is($$minimizer{kmers}{CGGTTCCTCGAATGATTACCC}[1], "TTAAAGGTCCGGTTCCTCGAATGATTACCCT", "revcom kmer");
  is($$minimizer{kmers}{ATAGAGGCCTGCGTAGCGGTG}[2], "GGTCTGGAATAGAGGCCTGCGTAGCGGTGAA", "revcom kmer");
  is($$minimizer{kmers}{ATAGAGGCCTGCGTAGCGGTG}[3], "GTCTGGAATAGAGGCCTGCGTAGCGGTGAAG", "revcom kmer");
  is($$minimizer{kmers}{CCCTCACGTACCGACTAAAGT}[4], "GTCTGTCCCTCACGTACCGACTAAAGTATTA", "revcom kmer");
  is($$minimizer{kmers}{CCCTCACGTACCGACTAAAGT}[5], "TCTGTCCCTCACGTACCGACTAAAGTATTAA", "revcom kmer");
};

