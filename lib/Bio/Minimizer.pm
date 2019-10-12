#!/usr/bin/env perl

# Minimizer.pm: a minimizer package for kmers
# Author: Lee Katz <lkatz@cdc.gov>

package Bio::Minimizer;
require 5.12.0;
our $VERSION=0.4;

use strict;
use warnings;
use Data::Dumper;
use Carp qw/carp croak/;

# TODO if/when threads, multithread based on large substrings
# of DNA b/c future heuristic 'saving lmers' between adjacent
# kmers.
our $iThreads; # boolean for whether threads are loaded
BEGIN{
  eval{
    require threads;
    require threads::shared;
    $iThreads = 1;
  };
  if($@){
    $iThreads = 0;
  }
}

my $startTime = time();
sub logmsg{
  local $0 = basename $0; 
  my $tid = 0;
  if($iThreads){
    $tid = threads->tid;
  }
  my $elapsedTime = time() - $startTime;
  print STDERR "$0.$tid $elapsedTime @_\n";
}

=pod

=head1 NAME

Bio::Minimizer - minimizer package

Based on the ideas put forth by Roberts et al 2004:
https://academic.oup.com/bioinformatics/article/20/18/3363/202143

=head1 SYNOPSIS

    my $minimizer = Bio::Minimizer->new($sequenceString);
    my $kmers     = $minimizer->kmers;     # hash of minimizer => kmer
    my $minimizers= $minimizer->minimizers;# hash of minimizer => [kmer1,kmer2,...]

    # With more options
    my $lmer      = Bio::Minimizer->new($sequenceString,{k=>31,l=>21});

=head1 DESCRIPTION

Creates a set of minimizers from sequence

=head1 VARIABLES

=over

=item $Bio::Minimizer::iThreads

Boolean describing whether the module instance is using threads

=back

=head1 METHODS

=over

=item Bio::Minimizer->new()

    Arguments:
      sequence     A string of ACGT
      settings     A hash
        k          Kmer length
        l          Minimizer length (some might call it lmer)
        numcpus    Number of threads to use. Does not currently do anything.

=back

=cut

sub new{
  my($class,$sequence,$settings) = @_;

  # Extract from $settings or set defaults
  my $k = $$settings{k} || 31;
  my $l = $$settings{l} || 21;
  my $numcpus = $$settings{numcpus} || 1;

  # Alter the sequence a bit
  $sequence = uc($sequence); # work in uppercase only
  $sequence =~ s/\s+//g;     # Remove all whitespace

  my $self={
    sequence   => $sequence,
    k          => $k,        # kmer length
    l          => $l,        # minimizer length
    numcpus    => $numcpus,  
    
    # Filled in by _minimizers()
    minimizers => {},        # kmer      => minimizer
    kmers      => {},        # minimizer => [kmer1,kmer2,...]
  };

  bless($self,$class);

  # Set $$self{minimizers} right away
  $self->_minimizers($sequence);

  return $self;
}

# Argument: string of nucleotides
sub _minimizers{
  my($self,$seq) = @_;
  my %MINIMIZER; 
  my %KMER;

  my ($k,$l)=($$self{k}, $$self{l}); 
  my $defaultSmallestLmer = 'Z' x $l;

  # Create a small array of lmers along the way
  # so that they don't have to be recalculated
  # all the time between kmers.
  my @lmer;

  # How many minimizers we'll get per kmer: the difference in lengths, plus 1
  my $minimizersPerKmer = $k-$l+1;

  # Also reverse-complement the sequence
  my $revcom = reverse($seq);
  $revcom =~ tr/ATCGatcg/TAGCtagc/;

  for my $sequence($seq, $revcom){
    # Number of kmers in the seq is the length of the seq, minus $k, plus 1
    my $numKmers = length($sequence) - $k + 1;
    for(my $i=0; $i<$numKmers; $i++){

      # The kmer is the subsequence starting at $i, length $k
      my $kmer=substr($sequence,$i,$k);
      
      # Get lmers along the length of the sequence into the @lmer buffer.
      # The start counter $j how many lmers are already in the buffer.
      for(my $j=scalar(@lmer); $j < $minimizersPerKmer; $j++){
        # The lmer will start at $i plus how many lmers are already
        # in the buffer @lmer, for length $l.
        my $lmer = substr($sequence, $i+$j, $l);
        push(@lmer, $lmer);
      }

      # The minimizer is the lowest lmer lexicographically sorted.
      my $minimizer = (sort {$a cmp $b} @lmer)[0];

      # Remove one lmer to reflect the step size of one
      # for the next iteration of the loop.
      shift(@lmer);

      $MINIMIZER{$kmer} = $minimizer;
      push(@{ $KMER{$minimizer} }, $kmer);
    } 
  }

  $$self{minimizers} = \%MINIMIZER;
  $$self{kmers}      = \%KMER;

  # Go ahead and return kmer=>minimizer
  return \%MINIMIZER;
}
 
1;
