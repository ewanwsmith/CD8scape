#!/usr/bin/env perl
use strict;
use warnings;

my $input = $ARGV[0];
my $output = $input;
$output =~ s/netMHCpan_output\.tsv/context_processed_netMHCpan_output\.csv/;

open(my $in, '<', $input) or die "Cannot open $input: $!\n";
open(my $out, '>', $output) or die "Cannot open $output: $!\n";

# Read and process the allele line
my $allele_line = <$in>;
chomp $allele_line;
$allele_line =~ s/^\s+|\s+$//g;  # Trim leading/trailing whitespace
my @alleles = grep { /\S/ } split /\t/, $allele_line;
my $num_haplotypes = scalar(@alleles);

# Read and process the header line
my $header_line = <$in>;
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

print $out join(",", @final_header), "\n";

# Process data rows
while (my $line = <$in>) {
    chomp $line;
    $line =~ s/^\s+|\s+$//g;
    my @cols = split /\t/, $line;
    my @fixed_cols = @cols[0..2];

    for (my $i = 0; $i < $num_haplotypes; $i++) {
        my $start_idx = $haplotype_start + ($i * $cols_per_haplotype);
        my @haplotype_data;
        for (my $j = 0; $j < $cols_per_haplotype; $j++) {
            my $idx = $start_idx + $j;
            $haplotype_data[$j] = (defined $cols[$idx]) ? $cols[$idx] : 'NA';
        }
        print $out join(",", @fixed_cols, $alleles[$i], @haplotype_data), "\n";
    }
}

close $in;
close $out;
print "Processed output written to $output\n";
