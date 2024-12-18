#!/usr/bin/perl
use strict;
use warnings;

# Usage: perl process_output.pl input.tsv > output.csv
#
# This script can handle two patterns for each haplotype:
# - 6 columns per haplotype: core, icore, EL-score, EL_Rank, BA-score, BA_Rank
# - 4 columns per haplotype: core, icore, EL-score, EL_Rank
#
# It determines which pattern to use by checking the first haplotype’s columns.
# If "BA-score" is present in the first haplotype block, it uses 6 columns.
# Otherwise, it uses 4 columns.
#
# The script also:
# - Trims whitespace and removes empty fields from the allele line.
# - Dynamically handles leftover columns after haplotypes.

my $input_file = shift or die "Usage: $0 input_file.tsv > output.csv\n";
open my $fh, '<', $input_file or die "Could not open '$input_file': $!";

# Read and process the allele line
my $allele_line = <$fh>;
chomp $allele_line;
$allele_line =~ s/^\s+|\s+$//g;   # Trim leading/trailing whitespace
my @alleles = grep { /\S/ } split /\t/, $allele_line;
my $num_haplotypes = scalar(@alleles);

# Read and process the header line
my $header_line = <$fh>;
chomp $header_line;
$header_line =~ s/^\s+|\s+$//g;
my @header_cols = split /\t/, $header_line;

# Basic columns: Pos, Peptide, ID must be present
die "Not enough columns for Pos, Peptide, ID\n" if @header_cols < 3;
my @final_header = @header_cols[0..2];

# If there are no haplotypes, just print the existing header and rows as CSV
if ($num_haplotypes == 0) {
    # Just print header and the rest as is
    print join(",", @final_header), "\n";
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/^\s+|\s+$//g;
        my @cols = split /\t/, $line;
        print join(",", @cols), "\n";
    }
    close $fh;
    exit;
}

my $total_cols = scalar(@header_cols);
my $haplotype_start = 3; # where haplotype columns start
my $available_haplotype_cols = $total_cols - 3;

# We need to determine if we have a 4-column or 6-column pattern.
# Let's check the first haplotype's columns:
# We'll attempt to detect BA-score. If present, we assume 6 columns.
# If not present, we assume 4 columns.

# Check we have enough columns to at least test a block
if ($available_haplotype_cols < 4) {
    die "Not enough columns to even form one haplotype block of 4 columns.\n";
}

# Extract a tentative first haplotype block of 4 columns
my @first_block = @header_cols[$haplotype_start .. ($haplotype_start + 3)];

# Check if we have BA-score after the first 4 columns.
# If we do and it's in the pattern, we will assume 6 columns per haplotype.
my $cols_per_haplotype;

if ($available_haplotype_cols >= $num_haplotypes * 6) {
    # We can potentially afford 6 columns per haplotype. Check if BA-score is in the first 6 columns.
    my @first_6_block = @header_cols[$haplotype_start .. ($haplotype_start + 5)];
    if (grep { $_ eq 'BA-score' } @first_6_block) {
        $cols_per_haplotype = 6;
    } else {
        # If BA-score not found in the first 6 columns, fall back to 4 columns
        $cols_per_haplotype = 4;
    }
} else {
    # Not enough columns for 6 columns per haplotype even if we wanted to.
    # Must be 4 columns if possible.
    $cols_per_haplotype = 4;
}

# Now verify we have enough columns for the chosen pattern
my $expected_allele_cols = $num_haplotypes * $cols_per_haplotype;
if ($available_haplotype_cols < $expected_allele_cols) {
    die "Not enough columns in the header line for $num_haplotypes haplotype(s) using $cols_per_haplotype columns each.\n";
}

# Construct headers
my $allele_end = $haplotype_start + $expected_allele_cols - 1;

my $allele_index = 0;
for (my $j = $haplotype_start; $j <= $allele_end; $j += $cols_per_haplotype) {
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

# Print the rest of the lines as CSV (convert tabs to commas)
while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/^\s+|\s+$//g;
    my @cols = split /\t/, $line;
    print join(",", @cols), "\n";
}

close $fh;