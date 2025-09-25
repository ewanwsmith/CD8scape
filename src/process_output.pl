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

# Determine whether haplotypes use 4 or 6 columns
my $cols_per_haplotype;
if ($available_haplotype_cols >= $num_haplotypes * 6) {
    my @first_6_block = @header_cols[$haplotype_start .. ($haplotype_start + 5)];
    if (grep { $_ eq 'BA-score' } @first_6_block) {
        $cols_per_haplotype = 6;
    } else {
        $cols_per_haplotype = 4;
    }
} else {
    $cols_per_haplotype = 4;
}

# Verify enough columns for haplotype data
my $expected_allele_cols = $num_haplotypes * $cols_per_haplotype;
die "Not enough columns in the header line for $num_haplotypes haplotype(s) using $cols_per_haplotype columns each.\n"
    if ($available_haplotype_cols < $expected_allele_cols);

# Construct new header
my @final_header = ('Pos', 'Peptide', 'ID', 'HLA', 'core', 'icore', 'EL-score', 'EL_Rank');
push @final_header, ('BA-score', 'BA_Rank') if $cols_per_haplotype == 6;
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
        # Replace any undefined values with empty string
        for (my $j = 0; $j < @haplotype_data; $j++) {
            $haplotype_data[$j] = defined $haplotype_data[$j] ? $haplotype_data[$j] : "";
        }
        my $allele = defined $alleles[$i] ? $alleles[$i] : "";
        print join(",", @fixed_cols, $allele, @haplotype_data), "\n";
    }
}

close $fh;
