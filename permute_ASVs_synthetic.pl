#!/usr/bin/perl -w

# example usage:
# perl permute_ASVs_synthetic.pl 84 10 1000 > permute.ASVs.Phyla84.ratio10.seq1000.txt

use strict;

my $N_genera  = $ARGV[0]; # the number of genera present in the sample. Here we will use genera, but this could be any taxonomic level (the higher of the two)
my $ASV_per_g = $ARGV[1]; # the number of ASVs per genus. Here we will use ASVs, but this could be any taxonomic level nested within a higher level
my $N_reads = $ARGV[2]; # the number of reads sampled

my %taxonomy; #keys=ASV_ID, values=taxonomy
my @ASVs; # list of all ASVs

# Assign ASVs to genera
my $ASV_ID = 0;
for (my $g=1; $g<=$N_genera; $g++) {
  for (my $n=1; $n<=$ASV_per_g; $n++) {
    $taxonomy{$ASV_ID} = $g;
    $ASV_ID++;
  }
}

# Print headers
print join "\t", "true_ASVs", "true_G", "obs_ASVs", "obs_G";
print "\n";

# Print the true number of ASVs and genera in each sample, in a range of samples containing from 1 to $N_genera
for (my $g=1; $g<=$N_genera; $g++) {
  my @ASVs;
  for (my $genus=1; $genus<=$g; $genus++) {
    for (my $n=0; $n<$ASV_ID; $n++) {
      next unless ($taxonomy{$n} == $genus);
      push @ASVs, $n;
    }
  }
  print join "\t", (@ASVs+0), $g;
  print "\t";
  
  # in each sample, count the total number of observd unique ASVs and observed unique genera at a given sampling effort ($N_reads)
  my $N_ASV = 0;
  my $N_genus = 0;
  my %seen_ASV;
  my %seen_genus;
  
  for (my $j=0; $j<$N_reads; $j++) { # sample a set number of reads at random from the ASVs in the sample
      my $j = int rand ( @ASVs+0 );
      my $ASV = $ASVs[$j];
      my $genus = $taxonomy{$ASV};
      unless (exists $seen_ASV{$ASV}) {
	  $seen_ASV{$ASV} = 1;
	  $N_ASV++;
      }
      unless (exists $seen_genus{$genus}) {
	  $seen_genus{$genus} = 1;
	  $N_genus++
      }
  }
  print join "\t", $N_ASV, $N_genus;
  print "\n";
}
