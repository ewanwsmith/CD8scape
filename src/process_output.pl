#!/usr/bin/perl
use strict;
use warnings;

my $input_file = shift or die "Usage: $0 input_file.tsv > output.csv\n";
open my $fh, '<', $input_file or die "Could not open '$input_file': $!";

# Read the first line (allele line), trim whitespace, split on tabs, remove empty fields
my $allele_line = <$fh>;
chomp $allele_line;
$allele_line =~ s/^\s+|\s+$//g;   # Trim leading and trailing whitespace
my @alleles = grep { /\S/ } split /\t/, $allele_line;

my $num_haplotypes = scalar(@alleles);

# Read the second line for headers
my $header_line = <$fh>;
chomp $header_line;
$header_line =~ s/^\s+|\s+$//g;   # Also trim just in case
my @header_cols = split /\t/, $header_line;

my @final_header;
push @final_header, $header_cols[0]; # Pos
push @final_header, $header_cols[1]; # Peptide
push @final_header, $header_cols[2]; # ID

my $cols_per_haplotype = 6;
my $expected_allele_cols = $num_haplotypes * $cols_per_haplotype;
my $total_cols = scalar(@header_cols);

# Check if we have enough columns for all haplotypes (if we have haplotypes at all)
if ($num_haplotypes > 0 && ($total_cols - 3) < $expected_allele_cols) {
    die "Not enough columns in the header line for $num_haplotypes haplotype(s) with 6 columns each.\n";
}

my $allele_start = 3;
my $allele_end   = $allele_start + $expected_allele_cols - 1;

# Attach allele names to their columns
my $allele_index = 0;
for (my $j = $allele_start; $j <= $allele_end && $num_haplotypes > 0; $j += $cols_per_haplotype) {
    my @block = @header_cols[$j..($j+$cols_per_haplotype-1)];
    my $allele = $alleles[$allele_index++];
    foreach my $field (@block) {
        push @final_header, $allele . "_" . $field;
    }
}

# Any leftover columns after haplotypes
if ($allele_end < ($total_cols - 1)) {
    my @leftover = @header_cols[($allele_end+1) .. ($total_cols - 1)];
    push @final_header, @leftover;
}

# Print the new header line as CSV
print join(",", @final_header), "\n";

# Print the rest of the lines as CSV, converting tabs to commas
while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/^\s+|\s+$//g; # Trim whitespace
    my @cols = split /\t/, $line;
    print join(",", @cols), "\n";
}

close $fh;