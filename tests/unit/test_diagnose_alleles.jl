using Test, CSV, DataFrames

@testset "diagnose_alleles parity (self-contained)" begin
    # Create synthetic best_ranks-like DataFrame and a squished panel DataFrame
    br = DataFrame(Allele = ["HLA-A01:01", "HLA-A*02:01", " HLA-A*03:01 "])
    spl = DataFrame(allele = ["HLA-A01:01", "HLA-A02:01", "HLA-A03:01"])  # squished form

    function normalize_allele(a)
        s = String(a)
        s = replace(s, '\u00A0' => ' ')
        s = uppercase(strip(s))
        s = replace(s, r"\s+" => "")
        # normalize by removing '*' so HLA-A*02:01 and HLA-A02:01 match
        s = replace(s, "*" => "")
        return s
    end

    best = Set(normalize_allele.(unique(string.(br.Allele))))
    panel_alleles = Set(normalize_allele.(String.(spl[!, 1])))
    @test issubset(best, panel_alleles)
end
