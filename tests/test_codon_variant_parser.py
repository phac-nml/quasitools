"""
Copyright Government of Canada 2017

Written by: Camy Tran, National Microbiology Laboratory,
            Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import pytest
import os
from quasitools.parsers.codon_variant_file_parser import parse_codon_variants
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.codon_variant import CodonVariant, CodonVariantCollection


TEST_PATH = os.path.dirname(os.path.abspath(__file__))

class TestCodonVariantParser:
    def test_valid_csv_file(self):
        """Tests to make sure that a valid codon variant csv file is properly
        parsed into a CodonVariantCollection object.
        """

        reference = TEST_PATH + "/data/hxb2_pol.fas"
        rs = parse_references_from_fasta(reference)

        var_obj = CodonVariantCollection(rs)

        for i in range(0, 30):
            variant = CodonVariant(
                    chrom="hxb2_pol",
                    pos=i,
                    gene="gag",
                    nt_start_gene=1309+i,
                    nt_end_gene=2841+i,
                    nt_start=2077+i,
                    nt_end=2079+i,
                    ref_codon="ata",
                    mutant_codon="aAa",
                    ref_aa="I",
                    mutant_aa="K",
                    coverage=563+i,
                    mutant_freq=1.60+i,
                    mutant_type="S",
                    ns_count=1.0000,
                    s_count=1.5000)

            pos = int(variant.nt_start)-int(variant.nt_start_gene)
            var_obj.variants["gag"][pos]["aAa"] = variant

        valid_csv = TEST_PATH + "/data/valid_csv.csv"

        with open(valid_csv, "w+") as f:
            f.write("#gene,nt position (gene),nt start position,"
                "nt end position,ref codon,mutant codon,ref AA,mutant AA,"
                "coverage,mutant frequency,mutant type,NS count,S count")
            
            for gene in var_obj.variants:
                for pos in var_obj.variants[gene]:
                    for codon in var_obj.variants[gene][pos]:
                        variant = var_obj.variants[gene][pos][codon]

                        f.write("%s,%i-%i,%i,%i,%s,%s,%s,%s,%i,%.2f,%s,%0.4f,%0.4f\n" % (
                            variant.gene, variant.nt_start_gene, variant.nt_end_gene,
                            variant.nt_start, variant.nt_end, variant.ref_codon,
                            variant.mutant_codon, variant.ref_aa, variant.mutant_aa,
                            variant.coverage, variant.mutant_freq, variant.mutant_type,
                            variant.ns_count, variant.s_count))


        parsed_codon_variants = parse_codon_variants(valid_csv, rs)

        for gene in parsed_codon_variants.variants:
            for pos in parsed_codon_variants.variants[gene]:
                for codon in parsed_codon_variants.variants[gene][pos]:
                    parsed_variant = parsed_codon_variants.variants[gene][pos][codon]
                    variant = var_obj.variants[gene][pos][codon]

                    assert parsed_variant.chrom == variant.chrom
                    assert parsed_variant.nt_start_gene == variant.nt_start_gene
                    assert parsed_variant.nt_end_gene == variant.nt_end_gene
                    assert parsed_variant.nt_start == variant.nt_start
                    assert parsed_variant.nt_end == variant.nt_end
                    assert parsed_variant.ref_codon == variant.ref_codon
                    assert parsed_variant.mutant_codon == variant.mutant_codon
                    assert parsed_variant.ref_aa == variant.ref_aa
                    assert parsed_variant.mutant_aa == variant.mutant_aa
                    assert parsed_variant.coverage == variant.coverage
                    assert parsed_variant.mutant_freq == variant.mutant_freq
                    assert parsed_variant.mutant_type == variant.mutant_type
                    assert parsed_variant.ns_count == variant.ns_count
                    assert parsed_variant.s_count == variant.s_count

        os.remove(valid_csv)

