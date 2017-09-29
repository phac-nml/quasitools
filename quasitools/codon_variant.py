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
from numpy import array as np_array

class CodonVariant(Variant):

    def __init__(self, gene=".", gene_start=0, gene_end=0,,
                 codon_start=0, codon_end=0,
                 ref_codon=".", mutant_codon=".",
                 ref_aa=".", mutant_aa=".", coverage=0,
                 mutant_freq="0", mutant_type=".",
                 ns_count=0, s_count=0)
        super(CodonVariant, self).__init__(*kwargs)

        self.gene = gene
        self.gene_start = gene_start
        self.gene_end = gene_end
        self.codon_start = codon_start
        self.codon_end = codon_end
        self.ref_codon = ref_codon
        self.mutant_codon = mutant_codon
        self.ref_aa = ref_aa
        self.mutant_aa = mutant_aa
        self.coverage = coverage
        self.mutant_freq = mutant_freq
        self.mutant_type = mutant_type
        self.ns_count = ns_count
        self.s_count = s_count

    def to_csv_entry(self):
        return "%s,%i-%i,%i,%i,%s,%s,%s,%s,%i,%.2f,%s,%0.4f,%0.4f\n" % (
            self.gene, self.gene_start, self.gene_end,
            self.codon_start, self.codon_end, self.ref_codon,
            self.mutant_codon, self.ref_aa, self.mutant_aa,
            self.coverage, self.mutant_freq, self.mutant_type,
            self.ns_count, self.s_count
        )

class CodonVariantCollection(VariantCollection):

    def __init__(self, references):
        super(CodonVariantCollection, self).__init__(references)
        self.variants = defaultdict(
            lambda: defaultdict(lambda: defaultdict(dict)))

    @classmethod
    def from_aa_census(cls, aa_census, *args):
        """Build the CodonVariantCollection from any number of
        AACensus objects"""
        
        # Handle aa_census being an array or single element nicely
        aa_census = np_array([aa_census]).flatten()

        var_collect = cls(aa_census)

        # Build up the collection of CodonVariants from many census
        for census_ind, census in enumerate(aa_census):

            # For each gene in this census
            for gene key in census.genes:
                gene = census.genes[gene_key]
                frame = gene['frame']


















