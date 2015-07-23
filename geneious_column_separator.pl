#!/usr/bin/perl

# Simple Perl script to select desired columns from full Geneious TSV export

use File::Basename;
#use strict;
use warnings;

# Read provided file name
my $GENEIOUSFILE = $ARGV[0];
my $GENEIOUSFILENAME = basename($GENEIOUSFILE, ".tsv");
# Check if it is readable
open(FH, "< $GENEIOUSFILE") || die "Cannot open $GENEIOUSFILE for reading: $!";
my @array = <FH>;
close FH || die "Could not open file: $!";

# Set name of output file
my $OUTPUTFILE = "$GENEIOUSFILENAME.columns.tsv";

open(MYFILE, "< $GENEIOUSFILE") || die "Cannot open $GENEIOUSFILE for read access: $!";

# Create new file for writing of output
open($OUTPUTFILE, "> $OUTPUTFILE") || die "Cannot create file for output: $!";

# List of desired columns
my @wanted_fields = ("# Sequences", "% Pairwise Identity", "Description", "Mean Coverage", "Name", "Sequence Length");

# Retrieve wanted fields from first line, as it is the header
my @fields = split /\t/, <MYFILE>;
chomp @fields;

# Write output file
while(<MYFILE>)
{
  chomp;
  # Read each line from the array and turn it into a hash and assign each column to matching header
  my %row;
  @row{@fields} = split /\t/;
  # Use map to find wanted column in the hash for matching data on row
  my @wanted_data = map{$row{$_}} @wanted_fields;
  print $OUTPUTFILE join("\t", @wanted_data), "\n";
  }

close $OUTPUTFILE || die "Error closing $OUTPUTFILE: $!";

exit;
