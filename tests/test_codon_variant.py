"""
Copyright Government of Canada 2015

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

import pdb
import pytest
import os
from quasitools.codon_variant import CodonVariant, CodonVariantCollection
from quasitools.nt_variant import NTVariant, NTVariantCollection
from quasitools.aa_census import AACensus
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.codon_variant_file_parser import parse_codon_variants
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.genes_file_parser import parse_genes_file


TEST_PATH = os.path.dirname(os.path.abspath(__file__))

class TestCodonVariant:
    @classmethod
    def setup_class(self):
        self.offset = 1269
        self.variant = CodonVariant(
                        chrom="hxb2_pol",
                        pos=1,
                        gene="gag",
                        nt_start_gene=1309,
                        nt_end_gene=2841,
                        nt_start=2077,
                        nt_end=2079,
                        ref_codon="ata",
                        mutant_codon="aAa",
                        ref_aa="I",
                        mutant_aa="K",
                        coverage=563,
                        mutant_freq=1.60,
                        mutant_type="S",
                        ns_count=1.0000,
                        s_count=1.5000)
    
    def test_to_csv_entry(self):
        assert self.variant.to_csv_entry(self.offset) == ( 
            "gag,%i-%i,%i,%i,ata,aAa,I,K,563,1.60,S,1.0000,1.5000\n" % (
            1309+self.offset, 2841+self.offset, 2077+self.offset, 2079+self.offset) )

class TestCodonVariantCollection:
    @classmethod
    def setup_class(self):
        self.reference = TEST_PATH + "/data/hxb2_pol.fas"
        self.references = parse_references_from_fasta(self.reference)
        self.variant_collection = CodonVariantCollection(self.references)
        self.offset = 1269

        self.variant_collection.variants['gag']['3']['aTa'] = CodonVariant(
                                                                chrom="hxb2_pol",
                                                                pos=1,
                                                                gene="gag",
                                                                nt_start_gene=1309,
                                                                nt_end_gene=2841,
                                                                nt_start=2077,
                                                                nt_end=2079,
                                                                ref_codon="ata",
                                                                mutant_codon="aTa",
                                                                ref_aa="I",
                                                                mutant_aa="K",
                                                                coverage=563,
                                                                mutant_freq=1.60,
                                                                mutant_type="S",
                                                                ns_count=1.0000,
                                                                s_count=1.5000)
        self.variant_collection.variants['tat']['10']['aAa'] = CodonVariant(
                                                                chrom="hxb2_pol",
                                                                pos=2,
                                                                gene="tat",
                                                                nt_start_gene=3309,
                                                                nt_end_gene=4841,
                                                                nt_start=4000,
                                                                nt_end=4002,
                                                                ref_codon="ata",
                                                                mutant_codon="aAa",
                                                                ref_aa="I",
                                                                mutant_aa="K",
                                                                coverage=563,
                                                                mutant_freq=1.60,
                                                                mutant_type="S",
                                                                ns_count=1.0000,
                                                                s_count=1.5000)

    def test_from_aacensus(self):
        bam = TEST_PATH + "/data/align.bam"
        genes_file = TEST_PATH + "/data/hxb2_pol.bed"
        mapped_read_collection_arr = []
        error_rate = 0.0038

        # Create a MappedReadCollection object
        for r in self.references:
            mapped_read_collection_arr.append(parse_mapped_reads_from_bam(r, bam))

            variants = NTVariantCollection.from_mapped_read_collections(
                    error_rate, self.references, *mapped_read_collection_arr)
            variants.filter('q30', 'QUAL<30', True)
            variants.filter('ac5', 'AC<5', True)
            variants.filter('dp100', 'DP<100', True)

        # Mask the unconfident differences
        for mrc in mapped_read_collection_arr:
            mrc.mask_unconfident_differences(variants)
        
        # Parse the genes from the gene file
        genes = parse_genes_file(genes_file, self.references[0].name)

        # Determine which frames our genes are in
        frames = set()

        for gene in genes:
            frames.add(genes[gene]['frame'])

        aa_census = AACensus(self.reference, mapped_read_collection_arr, genes, frames)

        test_variants = CodonVariantCollection.from_aacensus(
                            aa_census, next(iter(frames)))
        ref_seq = self.references[0].seq

        for gene in test_variants.variants:
            assert gene in genes
            for pos in test_variants.variants[gene]:
                for frame in frames:
                    nt_pos = pos/3 - frame
                    assert nt_pos >= genes[gene]['start'] or nt_pos <= genes[gene]['end'] 
                for codon in test_variants.variants[gene][pos]:
                    ref_codon = ref_seq[(pos):(pos) + 3].lower()   
                    assert codon != ref_codon
                    
    def test_to_csv_file(self):
        csv_file_string = self.variant_collection.to_csv_file(self.offset)
        # pdb.set_trace()
        assert ( "#gene,nt position (gene),nt start position,nt end position,"
                 "ref codon,mutant codon,ref AA,mutant AA,coverage,mutant frequency,"
                 "mutant type,NS count,S count" ) in csv_file_string
        assert ( "gag,%i-%i,%i,%i,ata,aTa,I,K,563,1.60,S,1.0000,1.5000" % (
            1309+self.offset, 2841+self.offset, 2077+self.offset, 2079+self.offset) ) in csv_file_string
        assert ( "tat,%i-%i,%i,%i,ata,aAa,I,K,563,1.60,S,1.0000,1.5000" % (
            3309+self.offset, 4841+self.offset, 4000+self.offset, 4002+self.offset) ) in csv_file_string

    def test_report_dnds_values(self):
        valid_dnds_report = TEST_PATH + "/data/output/dnds_report.csv"
        csv = TEST_PATH + "/data/output/mutant_types.csv"
        ref_seq = self.references[0].seq
        codon_variants = parse_codon_variants(csv, self.references)

        # Read from file and make sure there are no empty lines
        with open(valid_dnds_report, "r") as input:
            valid_report = input.read()

        # Sort and filter for comparison
        valid_dnds_values = sorted(filter(None, 
            valid_report.split("\n")))

        # Create the report string
        test_report = codon_variants.report_dnds_values(ref_seq, self.offset)

        # Split into lines and sort
        test_values = sorted(test_report.split("\n"))

        assert len(valid_dnds_values) == len(test_values)
        
        # Compare each line in the test report to the valid report
        for pos in range(0, len(valid_dnds_values)):
            if valid_dnds_values[pos][0:1] != "#":
                assert valid_dnds_values[pos] == \
                    test_values[pos]
    