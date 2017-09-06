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
from quasitools.variants import Variants
from quasitools.variant import Variant
from quasitools.parsers.reference_parser import parse_reference_fasta

class TestVariants:
    @classmethod
    def setup_class(self):
        self.references = parse_reference_fasta('tests/data/ref1.fasta')
        self.variants_obj = Variants(0.01, self.references)

        self.variants_obj.variants['ref1']['3']['t'] = Variant(chrom='ref1', pos=3, ref='c', alt='t', qual=30, info={'DP':400,'AC':12,'AF':0.03})
        self.variants_obj.variants['ref1']['10']['a'] = Variant(chrom='ref1', pos=10, ref='a', alt='t', qual=23, info={'DP':200,'AC':7,'AF':0.035})

    def test_from_bam(self):
        #TODO: add actual test for Variants.from_bam method
        assert True

    def test_from_mapped_reads(self):
        #TODO: add actual test for Variants.from_mapped_reads method
        assert True

    def test_str(self):
        vcf_string = str(self.variants_obj)

        assert "##fileformat=VCFv4.2" in vcf_string
        assert "##fileDate=" in vcf_string
        assert "##source=quasitools" in vcf_string
        assert "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" in vcf_string
        assert "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">" in vcf_string
        assert "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">" in vcf_string
        assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" in vcf_string
        assert "ref1\t3\t.\tc\tt\t30\t.\tDP=400;AC=12;AF=0.0300" in vcf_string
        assert "ref1\t10\t.\ta\tt\t23\t.\tDP=200;AC=7;AF=0.0350" in vcf_string

    def test_calculate_variant_qual(self):
        qual1 = self.variants_obj._Variants__calculate_variant_qual(12, 400)
        qual2 = self.variants_obj._Variants__calculate_variant_qual(7, 200)

        assert qual1 == 30
        assert qual2 == 23

    def test_filter(self):
        variants_obj_copy = copy.deepcopy(self.variants_obj)

        variants_obj_copy.filter('q30', 'QUAL<30', True)
        assert variants_obj_copy.variants['ref1']['3']['t'].filter == 'PASS'
        assert variants_obj_copy.variants['ref1']['10']['a'].filter == 'q30'

        variants_obj_copy.filter('ac5', 'AC<5', True)
        assert variants_obj_copy.variants['ref1']['3']['t'].filter == 'PASS'
        assert variants_obj_copy.variants['ref1']['10']['a'].filter == 'q30'

        variants_obj_copy.filter('dp100', 'DP<100', True)
        assert variants_obj_copy.variants['ref1']['3']['t'].filter == 'PASS'
        assert variants_obj_copy.variants['ref1']['10']['a'].filter == 'q30'
