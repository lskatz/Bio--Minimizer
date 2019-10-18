#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Minimizer;

my $sequence = "";
while(<STDIN>){
  s/^\s+|\s+$//g;  # whitespace trim
  next if(/^>/);

  $sequence .= $_;
}
close STDIN;

$sequence =~ s/\s+//g; # remove any whitespace in the sequence

my $minimizer = Bio::Minimizer->new($sequence);
my $starts = $minimizer->starts;
for my $m(keys(%$starts)){
  print "$m\t";
  print join(",", @{$$starts{$m}})."\n";
}

