#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Minimizer;
use Data::Dumper;
use Getopt::Long qw/GetOptions/;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings, qw(help l updatefrequency=i)) or die $!;
  if($$settings{help}){
    print usage();
    return 0;
  }

  my $l = $$settings{l} || 21;
  my $updateFrequency  = $$settings{updatefrequency} || 100000;
  $updateFrequency *= 4; # reflect the number of lines per entry

  my @entry;
  while(my $id = <STDIN>){
    my $seq = <STDIN>;
    my $plus= <STDIN>;
    my $qual= <STDIN>;
    chomp($id, $seq, $plus, $qual);

    # Minimizer object
    my $minimizer = Bio::Minimizer->new($seq,{k=>length($seq)-1,l=>$l});
    my $theOnlyMinimizer = (sort{$a cmp $b} values(%{$$minimizer{minimizers}}))[0] || $seq || "";

    my $GC = $seq =~ tr/GC//;    # GC count
       $GC = $GC / length($seq); # GC frequency

    # combine the minimum minimizer with the entry, for
    # sorting later.
    # Save the entry as a string so that we don't have to
    # parse it later.
    my $entry = [$theOnlyMinimizer, $seq, $GC, "$id\n$seq\n$plus\n$qual\n"];
    push(@entry,$entry);

    if($. % $updateFrequency == 0){
      my $numEntries = $. / 4;
      print STDERR "Finished with $numEntries\n";
    }
  }
  close STDIN;

  print STDERR "Finished getting minimizers and fastq entries. Now sorting/printing.\n";

  # Sort entries by minimizer
  my @sorted = sort{
       $$a[0] cmp $$b[0]                   # lexicographical of the minimizer
    || length($$a[1]) <=> length($$b[1])   # length
    || $$a[2] <=> $$b[2]                   # GC content
    || $$a[1] cmp $$b[1]                   # lexicographical of the sequence itself
  } @entry;

  for my $e(@sorted){
    # Print entries in order of minimizer
    print $$e[-1];
  }

  return 0;
}
sub usage{
  "$0 [options] < file.fastq > sorted.fastq

  -l                21      lmer length
  --updatefrequency 100000  See an update message after X fastq entries
  "
}

