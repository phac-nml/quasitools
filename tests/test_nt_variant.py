"""
Copyright Government of Canada 2015

Written by: Eric Enns, National Microbiology Laboratory, Public Health Agency of Canada

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
import copy
from quasitools.nt_variant import NTVariant, NTVariantCollection
from quasitools.parsers.reference_parser import parse_references_from_fasta


class TestNTVariant:
    @classmethod
    def setup_class(self):
        self.variant = NTVariant('hxb2_pol', 1, ref='g', alt='t', qual='30', filter='PASS', info={
                                 'DP': 400, 'AC': 12, 'AF': 0.03})

    def test_to_vcf_entry(self):
        assert self.variant.to_vcf_entry(
        ) == 'hxb2_pol\t1\t.\tg\tt\t30\tPASS\tDP=400;AC=12;AF=0.0300'


class TestNTVariantCollection:
    @classmethod
    def setup_class(self):
        self.references = parse_references_from_fasta('tests/data/ref1.fasta')
        self.variant_collection = NTVariantCollection(self.references)

        self.variant_collection.variants['ref1']['3']['t'] = NTVariant(
            chrom='ref1', pos=3, ref='c', alt='t', qual=30, info={'DP': 400, 'AC': 12, 'AF': 0.03})
        self.variant_collection.variants['ref1']['10']['a'] = NTVariant(
            chrom='ref1', pos=10, ref='a', alt='t', qual=23, info={'DP': 200, 'AC': 7, 'AF': 0.035})

    def test_from_mapped_read_collections(self):
        # TODO: add actual test for Variants.from_mapped_reads method
        assert True

    def test_to_vcf_file(self):
        vcf_file_string = self.variant_collection.to_vcf_file()

        assert "##fileformat=VCFv4.2" in vcf_file_string
        assert "##fileDate=" in vcf_file_string
        assert "##source=quasitools" in vcf_file_string
        assert "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" in vcf_file_string
        assert "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">" in vcf_file_string
        assert "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">" in vcf_file_string
        assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" in vcf_file_string
        assert "ref1\t3\t.\tc\tt\t30\t.\tDP=400;AC=12;AF=0.0300" in vcf_file_string
        assert "ref1\t10\t.\ta\tt\t23\t.\tDP=200;AC=7;AF=0.0350" in vcf_file_string

    def test_calculate_variant_qual(self):
        qual1 = self.variant_collection._NTVariantCollection__calculate_variant_qual(
            0.01, 12, 400)
        qual2 = self.variant_collection._NTVariantCollection__calculate_variant_qual(
            0.01, 7, 200)

        assert qual1 == 30
        assert qual2 == 23

    def test_filter(self):
        variant_collection_copy = copy.deepcopy(self.variant_collection)

        variant_collection_copy.filter('q30', 'QUAL<30', True)
        assert variant_collection_copy.variants['ref1']['3']['t'].filter == 'PASS'
        assert variant_collection_copy.variants['ref1']['10']['a'].filter == 'q30'

        variant_collection_copy.filter('ac5', 'AC<5', True)
        assert variant_collection_copy.variants['ref1']['3']['t'].filter == 'PASS'
        assert variant_collection_copy.variants['ref1']['10']['a'].filter == 'q30'

        variant_collection_copy.filter('dp100', 'DP<100', True)
        assert variant_collection_copy.variants['ref1']['3']['t'].filter == 'PASS'
        assert variant_collection_copy.variants['ref1']['10']['a'].filter == 'q30'
