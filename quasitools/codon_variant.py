"""
Copyright Government of Canada 2017

Written by: Camy Tran, Eric Enns, National Microbiology Laboratory,
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

from collections import defaultdict
from quasitools.aa_census import CONFIDENT
from quasitools.variant import Variant, VariantCollection
from numpy import array as np_array
from Bio.Seq import Seq
from numpy import log


class CodonVariant(Variant):

    def __init__(self, gene, nt_start_gene, nt_end_gene,
                 nt_start, nt_end, ref_codon, mutant_codon,
                 ref_aa, mutant_aa, coverage, mutant_freq,
                 mutant_type, ns_count, s_count, **kwargs):
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

    @classmethod
    def from_aacensus(cls, gene_key, aa, codon, census, frame, nt_pos):
        codon_permutations = [
            [[0]], [[0, 1], [1, 0]],
            [
                [0, 1, 2], [0, 2, 1],
                [1, 0, 2], [1, 2, 0],
                [2, 0, 1], [2, 1, 0]
            ]
        ]

        gene = census.genes[gene_key]
        frame = gene['frame']
        chrom = gene['chrom']

        coverage = census.coverage_at(frame, nt_pos)
        ref_seq = census.mapped_read_collections[0].reference.seq
        ref_codon = ref_seq[(nt_pos*3+frame):
                            (nt_pos*3+frame) + 3].lower()
        ref_aa = Seq(ref_codon).translate()[0]

        if aa == ref_aa:
            mutation_type = "S"
        else:
            mutation_type = "NS"

        nt_change_count = sum(1 for c in codon if not c.islower())
        base_change_pos = []

        for codon_pos in range(0, 3):
            nucleotide = codon[
                codon_pos:codon_pos+1]
            if nucleotide.upper() == nucleotide:
                base_change_pos.append(codon_pos)

        ns_count = 0
        s_count = 0

        for codon_permutation in codon_permutations[nt_change_count-1]:
            codon_pathway = ref_codon
            for base_pos in codon_permutation:
                mutant_pos = base_change_pos[base_pos]
                mutant_nt = codon[mutant_pos:mutant_pos+1]
                codon_pathway = codon_pathway[:mutant_pos] + \
                    mutant_nt + codon_pathway[mutant_pos+1:]

                if Seq(codon_pathway).translate()[0] == ref_aa:
                    s_count += 1
                else:
                    ns_count += 1

        return cls(
            chrom=chrom,
            gene=gene_key,
            id="mutation",
            coverage=coverage,
            mutant_freq=census.codon_frequency_for_amino_at(
                frame,
                nt_pos,
                aa,
                CONFIDENT,
                codon)/float(coverage)*100.0,
            pos=(nt_pos - (gene['start'] // 3) + 1),
            nt_start_gene=gene['start'],
            nt_end_gene=gene['end'],
            nt_start=nt_pos*3 + frame,
            nt_end=nt_pos*3 + frame+2,
            ref_codon=ref_codon,
            mutant_codon=codon,
            ref_aa=ref_aa,
            mutant_aa=aa,
            mutant_type=mutation_type,
            ns_count=ns_count,
            s_count=s_count)


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

        # Build up the collection of CodonVariants from many census
        for census_ind, census in enumerate(aa_census):

            # For each gene in this census
            for gene_key in census.genes:
                gene = census.genes[gene_key]
                frame = gene['frame']
                gene_start = int((gene['start'] - frame) / 3)
                gene_end = int((gene['end'] - frame - 2) / 3)

                for nt_pos in range(gene_start, gene_end):
                    ref_seq = census.mapped_read_collections[0].reference.seq
                    ref_codon = ref_seq[(nt_pos*3+frame):
                                        (nt_pos*3+frame) + 3].lower()

                    for aa in census.aminos_at(frame, nt_pos, CONFIDENT):
                        frequency = census.amino_frequency_at(
                            frame, nt_pos, aa, CONFIDENT)
                        if frequency >= 0.01:
                            for codon in census.amino_to_codons_at(
                                                          frame, nt_pos,
                                                          aa, CONFIDENT):
                                    if codon != ref_codon:
                                        mutation = CodonVariant.from_aacensus(
                                                gene_key,
                                                aa, codon,
                                                census,
                                                frame,
                                                nt_pos)

                                        var_collect.variants[gene_key][
                                            (nt_pos*3 + frame)][
                                                codon] = mutation
        return var_collect

    def to_csv_file(self, offset):
        """"Build a string representation of our CodonVariant objects
        (i.e. a csv file)."""

        report = ("#gene,nt position (gene),nt start position,"
                  "nt end position,ref codon,mutant codon,"
                  "ref AA,mutant AA,coverage,mutant frequency,"
                  "mutant type,NS count,S count\n")

        for gene in self.variants:
            for pos in self.variants[gene]:
                for codon in self.variants[gene][pos]:
                    report += self.variants[gene][pos][
                        codon].to_csv_entry(offset)

        return report[:-1]

    def report_dnds_values(self, ref_seq, offset):
        report = "#gene,pn,ps,pn_sites,ps_sites,dn/ds\n"

        # Iterate through the variants to
        # create a gene map and report on each gene
        genes = defaultdict(lambda: defaultdict(lambda:
                            defaultdict(lambda: defaultdict(int))))

        for gene in self.variants:
            for pos in self.variants[gene]:
                for codon in self.variants[gene][pos]:
                    variant = self.variants[gene][pos][codon]

                    genes[gene]['start'] = variant.nt_start_gene
                    genes[gene]['end'] = variant.nt_end_gene

                    if variant.ns_count > 0:
                        genes[gene][pos]['NS'][variant.ns_count] += (
                            variant.mutant_freq/100.0)
                    if variant.s_count > 0:
                        genes[gene][pos]['S'][variant.s_count] += (
                            variant.mutant_freq/100.0)

        for gene in genes:
            s_sites = 0
            ns_sites = 0
            gene_seq = ref_seq[(genes[gene]['start'] - offset):
                               (genes[gene]['end'] - offset + 1)]

            pn = 0
            ps = 0

            ns_ncod = 0
            s_ncod = 0

            pn_ncod = 0
            ps_ncod = 0

            for i in range(0, len(gene_seq)-1, 3):
                codon = gene_seq[i:i+3]
                aa = Seq(codon).translate()[0]
                non_syn = 0

                # synonymous sites only occur at 1st and 3rd pos in a codon
                for j in range(0, 3):
                    for nt in ('a', 'c', 'g', 't'):
                        if nt.lower() != codon[j:j+1].lower():
                            mod_codon = codon[:j] + nt + codon[j+1:]
                            mod_aa = Seq(mod_codon).translate()[0]
                            if mod_aa.upper() != aa.upper():
                                non_syn += 1/3.0

                ns_sites += non_syn
                s_sites += 3-non_syn

                if non_syn > 0:
                    ns_ncod += 1
                if 3-non_syn > 0:
                    s_ncod += 1

                pni = 0
                psi = 0

                if 'NS' in genes[gene][i]:
                    for count in genes[gene][i]['NS']:
                        pni += genes[gene][i]['NS'][int(count)] * (
                            count/non_syn)
                    pn += pni
                    pn_ncod += 1

                if 'S' in genes[gene][i]:
                    for count in genes[gene][i]['S']:
                        psi += genes[gene][i]['S'][int(count)] * (
                            count/(3-non_syn))
                    ps += psi
                    ps_ncod += 1

            if pn_ncod > 0 and ps_ncod > 0:
                pn = pn/pn_ncod
                ps = ps/ps_ncod

                if (1-(4*pn/3.0) > 0) and (1-(4*ps/3.0) > 0):
                    dn = -(3/4.0)*log(1-(4*pn/3.0))
                    ds = -(3/4.0)*log(1-(4*ps/3.0))
                    report += "%s,%0.4f,%0.4f,%i,%i,%0.4f\n" % \
                        (gene, pn, ps, pn_ncod, ps_ncod, dn/ds)
                else:
                    report += "%s,%0.4f,%0.4f,%i,%i,N/A\n" % \
                            (gene, pn, ps, pn_ncod, ps_ncod)

        return report[:-1]
