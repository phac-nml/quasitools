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

import pdb
from collections import defaultdict
from quasitools.aa_census import CONFIDENT, UNCONFIDENT
from quasitools.variant import Variant, VariantCollection
from numpy import array as np_array
from Bio.Seq import Seq


class CodonVariant(Variant):

    def __init__(self, gene=".", nt_start_gene=0, nt_end_gene=0,
                 nt_start=0, nt_end=0,
                 ref_codon=".", mutant_codon=".",
                 ref_aa=".", mutant_aa=".", coverage=0,
                 mutant_freq="0", mutant_type=".",
                 ns_count=0, s_count=0, **kwargs):
        super(CodonVariant, self).__init__(**kwargs)

        self.gene = gene
        self.nt_start_gene = nt_start_gene
        self.nt_end_gene = nt_end_gene
        self.nt_start = nt_start
        self.nt_end = nt_end
        self.ref_codon = ref_codon
        self.mutant_codon = mutant_codon
        self.ref_aa = ref_aa
        self.mutant_aa = mutant_aa
        self.coverage = coverage
        self.mutant_freq = mutant_freq
        self.mutant_type = mutant_type
        self.ns_count = ns_count
        self.s_count = s_count

    def to_csv_entry(self, offset):
        return "%s,%i-%i,%i,%i,%s,%s,%s,%s,%i,%.2f,%s,%0.4f,%0.4f\n" % (
            self.gene, self.nt_start_gene+offset, self.nt_end_gene+offset,
            self.nt_start+offset, self.nt_end+offset, self.ref_codon,
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
    def from_aacensus(cls, aa_census, *args):
        """Build the CodonVariantCollection from any number of
        AACensus objects"""
        
        # Handle aa_census being an array or single element nicely
        aa_census = np_array([aa_census]).flatten()

        var_collect = cls(aa_census)

        codon_permutations = [[[0]], [[0,1],[1,0]], [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]]]


        # Build up the collection of CodonVariants from many census
        for census_ind, census in enumerate(aa_census):

            # For each gene in this census
            for gene_key in census.genes:
                gene = census.genes[gene_key]
                frame = gene['frame']
                codon_start = int((gene['start'] - frame) / 3) # don't understand the diff starts
                codon_end = int((gene['end'] - frame - 2) / 3)
                chrom = gene['chrom']

                for ref_codon_pos in range(codon_start, codon_end):
                    coverage = census.coverage_at(frame, ref_codon_pos)
                    ref_seq = census.mapped_read_collections[0].reference.seq
                    ref_codon = ref_seq[(ref_codon_pos*3+frame): (ref_codon_pos*3+frame) + 3].lower()
                    ref_aa = Seq(ref_codon).translate()[0]

                    for aa in census.aminos_at(frame, ref_codon_pos, CONFIDENT):
                        frequency = census.amino_frequency_at(frame, ref_codon_pos, aa, CONFIDENT)
                        if frequency >= 0.01:
                            # if (ref_codon_pos*3 + frame) == 7145:
                            #     pdb.set_trace()
                            for codon in census.amino_to_codons_at(frame, ref_codon_pos, aa, CONFIDENT):
                                # pdb.set_trace()
                                if codon != ref_codon:
                                    if aa == ref_aa:
                                        mutation_type = "S"
                                    else:
                                        mutation_type = "NS"

                                    base_change_count = sum(1 for c in codon if not c.islower())
                                    base_change_pos = []

                                    for codon_pos in range(0, 3):
                                        nucleotide = codon[codon_pos:codon_pos+1]
                                        if nucleotide.upper() == nucleotide:
                                            base_change_pos.append(codon_pos)

                                    ns_count = 0
                                    s_count = 0

                                    for codon_permutation in codon_permutations[base_change_count-1]:
                                        codon_pathway = ref_codon
                                        for base_pos in codon_permutation:
                                            mutant_pos = base_change_pos[base_pos]
                                            mutant_nt = codon[mutant_pos:mutant_pos+1]
                                            codon_pathway = codon_pathway[:mutant_pos] + mutant_nt + \
                                                    codon_pathway[mutant_pos+1:]

                                            if Seq(codon_pathway).translate()[0] == ref_aa:
                                                s_count += 1
                                            else:
                                                ns_count += 1

                                    # Create CodonVariant and add to collection
                                    mutation = CodonVariant(chrom=chrom,
                                                            gene=gene_key,
                                                            id="mutation",
                                                            ref=ref_codon,
                                                            alt=codon,
                                                            mutant_freq=census.codon_frequency_for_amino_at(frame, ref_codon_pos, aa, CONFIDENT, codon) / coverage*100,
                                                            coverage=coverage,
                                                            pos=(ref_codon_pos - (
                                                                gene['start'] // 3
                                                            ) + 1),
                                                            nt_start_gene=gene['start'],
                                                            nt_end_gene=gene['end'],
                                                            nt_start=(ref_codon_pos*3 + frame),
                                                            nt_end=(ref_codon_pos*3 + frame+2),
                                                            ref_codon=ref_codon,
                                                            mutant_codon=codon,
                                                            ref_aa=ref_aa,
                                                            mutant_aa=aa,
                                                            mutant_type=mutation_type,
                                                            ns_count=ns_count,
                                                            s_count=s_count
                                                            )
                                    var_collect.variants[gene_key][(ref_codon_pos*3 + frame)][codon] = mutation
        return var_collect

    def to_csv_file(self, offset):
        """"Build a string representation of our CodonVariant objects
        (i.e. a csv file)."""

        report = ("#gene,nt position (gene),nt start position,nt end position,ref codon,mutant codon,"
        "ref AA,mutant AA,coverage,mutant frequency,mutant type,NS count,S count\n")

        for gene in self.variants:
            for pos in self.variants[gene]:
                for codon in self.variants[gene][pos]:
                    report += self.variants[gene][pos][codon].to_csv_entry(offset)

        return report[:-1]

                            
                                                        
                                                        












