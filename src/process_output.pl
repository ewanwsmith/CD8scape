#!/usr/bin/perl
use strict;
use warnings;

# Usage: perl process_output.pl input.tsv > output.csv

my $input_file = shift or die "Usage: $0 input_file.tsv > output.csv\n";
open my $fh, '<', $input_file or die "Could not open '$input_file': $!";

# Read and process the allele line
my $allele_line = <$fh>;
chomp $allele_line;
$allele_line =~ s/^\s+|\s+$//g;  # Trim leading/trailing whitespace
my @alleles = grep { /\S/ } split /\t/, $allele_line;
my $num_haplotypes = scalar(@alleles);

# Read and process the header line
my $header_line = <$fh>;
chomp $header_line;
$header_line =~ s/^\s+|\s+$//g;
my @header_cols = split /\t/, $header_line;

die "Not enough columns for Pos, Peptide, ID\n" if @header_cols < 3;

# Fixed columns
my @fixed_header = @header_cols[0..2];

my $total_cols = scalar(@header_cols);
my $haplotype_start = 3;
my $available_haplotype_cols = $total_cols - 3;

# Determine per-haplotype column width robustly and extract the per-haplotype column names
my $cols_per_haplotype = 4;
if ($num_haplotypes > 0) {
    # Prefer exact division if possible
    if ($available_haplotype_cols % $num_haplotypes == 0) {
        $cols_per_haplotype = int($available_haplotype_cols / $num_haplotypes);
    } else {
        # Fallback heuristic: try 6 then 4
        if ($available_haplotype_cols >= $num_haplotypes * 6) {
            $cols_per_haplotype = 6;
        } else {
            $cols_per_haplotype = 4;
        }
    }
}

# Verify enough columns for haplotype data
my $expected_allele_cols = $num_haplotypes * $cols_per_haplotype;
die "Not enough columns in the header line for $num_haplotypes haplotype(s) using $cols_per_haplotype columns each.\n"
    if ($available_haplotype_cols < $expected_allele_cols);

# Extract per-haplotype column names from the header (use the first block)
my @per_haplotype_cols = @header_cols[$haplotype_start .. ($haplotype_start + $cols_per_haplotype - 1)];

# Construct new header using the actual per-haplotype column names
my @final_header = ('Pos', 'Peptide', 'ID', 'HLA', @per_haplotype_cols);
print join(",", @final_header), "\n";

# Process data rows
while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/^\s+|\s+$//g;
    my @cols = split /\t/, $line;
    my @fixed_cols = @cols[0..2];
    
    for (my $i = 0; $i < $num_haplotypes; $i++) {
        my $start_idx = $haplotype_start + ($i * $cols_per_haplotype);
        my @haplotype_data = @cols[$start_idx .. ($start_idx + $cols_per_haplotype - 1)];
        print join(",", @fixed_cols, $alleles[$i], @haplotype_data), "\n";
    }
}

close $fh;
