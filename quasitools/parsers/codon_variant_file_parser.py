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

from quasitools.codon_variant import CodonVariant, CodonVariantCollection


def parse_codon_variants(csv, references):
    """Parse a codon variants csv and build a codon variants object"""

    variant_collect = CodonVariantCollection(references)

    with open(csv, "r") as f:
        for line in f:
            if line[0] != "#":
                (
                    gene, gene_start_end, nt_start,
                    nt_end, ref_codon, mutant_codon,
                    ref_aa, mutant_aa, coverage,
                    mutant_freq, mutant_type,
                    ns_count, s_count
                ) = line.rstrip().split(",")

                gene_start, gene_end = gene_start_end.split('-')

                pos = int(nt_start)-int(gene_start)

                variant = CodonVariant(
                    chrom=gene,
                    pos=pos,
                    gene=gene,
                    nt_start_gene=int(gene_start),
                    nt_end_gene=int(gene_end),
                    nt_start=int(nt_start),
                    nt_end=int(nt_end),
                    ref_codon=ref_codon,
                    mutant_codon=mutant_codon,
                    ref_aa=ref_aa,
                    mutant_aa=mutant_aa,
                    coverage=int(coverage),
                    mutant_freq=float(mutant_freq),
                    mutant_type=mutant_type,
                    ns_count=float(ns_count),
                    s_count=float(s_count))

                variant_collect.variants[gene][
                    pos][mutant_codon] = variant

    f.close()
    return variant_collect
