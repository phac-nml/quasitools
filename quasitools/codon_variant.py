"""
Copyright Government of Canada 2017

Written by: Camy Tran, National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

from collections import defaultdict
from quasitools.aa_census import CONFIDENT, UNCONFIDENT
from quasitools.variant import Variant, VariantCollection

class CodonVariant(Variant):

    def __init__(self, gene=".", gene_nt_pos=0,
                 nt_start_pos=0, nt_end_pos=0,
                 ref_codon=".", mutant_codon=".",
                 ref_aa=".", mutant_aa=".", coverage=0,
                 mutant_freq="0", mutant_type=".",
                 ns_count=0, s_count=0)
        super(CodonVariant, self).__init__(*kwargs)

        self.gene = gene
        self.gene_nt_pos = gene_nt_pos
        self.nt_start_pos = nt_start_pos
        self.nt_end_pos = nt_end_pos
        self.ref_codon = ref_codon
        self.mutant_codon = mutant_codon
        self.ref_aa = ref_aa
        self.mutant_aa = mutant_aa
        self.coverage = coverage
        self.mutant = mutant_freq
        self.mutant_type = mutant_type
        self.ns_count = ns_count
        self.s_count = s_count
