#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Test::More tests=>6;
use FindBin qw/$RealBin/;
use lib "$RealBin/../lib";
use List::MoreUtils qw/uniq/;

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

subtest 'Minimizer => kmer' => sub{
  my %expected = (
    AATCAGCCCCGGTACCTTTGG =>
      [qw(ACAATCAGCCCCGGTACCTTTGGGTCTGTCC CAATCAGCCCCGGTACCTTTGGGTCTGTCCC AATCAGCCCCGGTACCTTTGGGTCTGTCCCT)], 
    AATCACTTTTAATACTTTAGT =>
      [qw(AGGGCCCCAGAATCACTTTTAATACTTTAGT GGGCCCCAGAATCACTTTTAATACTTTAGTC GGCCCCAGAATCACTTTTAATACTTTAGTCG GCCCCAGAATCACTTTTAATACTTTAGTCGG CCCCAGAATCACTTTTAATACTTTAGTCGGT CCCAGAATCACTTTTAATACTTTAGTCGGTA CCAGAATCACTTTTAATACTTTAGTCGGTAC CAGAATCACTTTTAATACTTTAGTCGGTACG AGAATCACTTTTAATACTTTAGTCGGTACGT GAATCACTTTTAATACTTTAGTCGGTACGTG)],
    AGCCCCGGTACCTTTGGGTCT =>
      [qw(ATCAGCCCCGGTACCTTTGGGTCTGTCCCTC TCAGCCCCGGTACCTTTGGGTCTGTCCCTCA)],
    CCTTTGGGTCTGTCCCTCACG =>
      [qw(CCTTTGGGTCTGTCCCTCACGTACCGACTAA)],
    AATACTTTAGTCGGTACGTGA =>
      [qw(AATCACTTTTAATACTTTAGTCGGTACGTGA ATCACTTTTAATACTTTAGTCGGTACGTGAG TCACTTTTAATACTTTAGTCGGTACGTGAGG CACTTTTAATACTTTAGTCGGTACGTGAGGG ACTTTTAATACTTTAGTCGGTACGTGAGGGA CTTTTAATACTTTAGTCGGTACGTGAGGGAC TTTTAATACTTTAGTCGGTACGTGAGGGACA TTTAATACTTTAGTCGGTACGTGAGGGACAG TTAATACTTTAGTCGGTACGTGAGGGACAGA TAATACTTTAGTCGGTACGTGAGGGACAGAC AATACTTTAGTCGGTACGTGAGGGACAGACC)],
    );

  plan tests=>scalar(keys(%expected));

  while(my($m,$kmers) = each(%expected)){
    my @got = sort {$a cmp $b} uniq(@$kmers);
    my @exp = sort {$a cmp $b} uniq(@{ $$minimizer{kmers}{$m} });
    is_deeply(\@got, \@exp, $m);
  }

};

my $numKmerTargets = 0;
my $numMinimizers  = 0;
while(my($min, $kmers) = each(%{$$minimizer{kmers}})){
  $numKmerTargets += scalar(@$kmers);
  $numMinimizers++;
}
is($numKmerTargets, 420, "Number of kmer targets");
is($numMinimizers, 78, "Number of minimizers");

