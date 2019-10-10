#!/usr/bin/env perl

# Minimizer.pm: a minimizer package for kmers
# Author: Lee Katz <lkatz@cdc.gov>

package Bio::Minimizer;
require 5.10.0;
our $VERSION=0.1;

use strict;
use warnings;

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

# TODO if 'die' is imported by a script, redefine
# sig die in that script as this function.
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

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

=head1 SYNOPSIS

    my $minimizer = Bio::Minimizer->new($sequenceString);
    my $kmers     = $minimizer->kmers;     # hash of minimizer => kmer
    my $minimizers= $minimizer->minimizers;# hash of minimizer => [kmer1,kmer2,...]

=head1 DESCRIPTION

Creates a set of minimizers from sequence

=head1 VARIABLES

=over

=item $Bio::Minimizer::iThreads

Boolean describing whether the module instance is using threads

=back

=cut

sub new{
  my($class,$sequence) = @_;

  my $self={
    sequence   => $sequence,
    k          => 31,        # kmer length
    l          => 21,        # minimizer length
    
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
  my %LMER; 
  my %KMER;

  my ($k,$l)=($$self{k}, $$self{l}); 
  my $minimizersPerKmer = $k-$l+1;
  $seq=~s/\s+//g;
  my $length=length($seq);
  my $numKmers = $length - $k + 1;
  for(my $i=0; $i<$numKmers; $i++){
    my $kmer=substr($seq,$i,$k);
    my $smallestMinimizer = ~0;
    for(my $j=0; $j<$minimizersPerKmer; $j++){
      my $minimizer = substr($kmer, $j, $l);
      if($minimizer lt $smallestMinimizer){
        $smallestMinimizer = $minimizer;
      } 
      $LMER{$kmer} = $minimizer;
      push(@{ $KMER{$minimizer} }, $kmer);
    } 
  } 

  $$self{minimizers} = \%LMER;
  $$self{kmers}      = \%KMER;

  return \%LMER;
}
 
