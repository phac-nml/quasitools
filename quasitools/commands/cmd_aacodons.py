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

import click
from Bio.Seq import Seq
from quasitools.cli import pass_context
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.genes_file_parser import parse_genes_file
from quasitools.nt_variant import NTVariantCollection
from quasitools.aa_census import AACensus, CONFIDENT

@click.command('aacodons', short_help='Generate a summary of codon to amino acid mutations.')
@click.argument('bam', required=True, type=click.Path(exists=True))
@click.argument('reference', required=True, type=click.Path(exists=True))
@click.argument('offset', required=True, type=float)
@click.argument('genes_file', required=True, type=click.Path(exists=True))
@click.argument('output', required=True)
@click.option('-e', '--error_rate', default=0.01, help='estimated sequencing error rate.')

@pass_context
def cli(ctx, bam, reference, offset, genes_file, output, error_rate):
    click.echo("Running cmd_aacodons command...")

    rs = parse_references_from_fasta(reference)
    mapped_read_collection_arr = []

    for r in rs:
        # Create a MappedReadCollection object
        mapped_read_collection_arr.append(parse_mapped_reads_from_bam(r, bam))

    variants = NTVariantCollection.from_mapped_read_collections(error_rate, rs, *mapped_read_collection_arr)
    variants.filter('q30', 'QUAL<30', True)
    variants.filter('ac5', 'AC<5', True)
    variants.filter('dp100', 'DP<100', True)

    click.echo("Creating variants vcf file...")
    vcf_filename = "{}.vcf".format(output.split('.')[0])
    vcf_file = open(vcf_filename, "w")
    vcf_file.write(variants.to_vcf_file())
    vcf_file.close()

    # Mask the unconfident differences
    for mrc in mapped_read_collection_arr:
        mrc.mask_unconfident_differences(vcf_filename)

    # Parse the genes from the gene file
    genes = parse_genes_file(genes_file, rs[0].name)

    # Determine which frames our genes are in
    frames = set()

    for gene in genes:
        frames.add(genes[gene]['frame'])

    # Create an AACensus object
    click.echo("Creating AACensus object...")
    aa_census = AACensus(reference, mapped_read_collection_arr, genes, frames)

    codon_permutations = [[[0]], [[0,1],[1,0]], [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]]]

    report = ("#gene,nt position (gene),nt start position,nt end position,ref codon,mutant codon,"
              "ref AA,mutant AA,coverage,mutant frequency,mutant type,NS count,S count\n")

    click.echo("Building report...")

    for gene in genes:
        frame = genes[gene]['frame']
        codon_start = int((genes[gene]['start'] - frame) / 3)
        codon_end = int((genes[gene]['end'] - frame - 2) / 3)

        for i in range(codon_start, codon_end):
            coverage = aa_census.coverage_at(frame, i)
            ref_codon = rs[0].seq[(i*3+frame): (i*3+frame) + 3].lower()
            ref_aa = Seq(ref_codon).translate()[0]

            for aa in aa_census.aminos_at(frame, i, CONFIDENT):
                frequency = aa_census.amino_frequency_at(frame, i, aa, CONFIDENT)
                if frequency >= 0.01:
                    for codon in aa_census.amino_to_codons_at(frame, i, aa, CONFIDENT):
                        if codon != ref_codon:
                            report += "%s,%i-%i,%i,%i,%s,%s,%s,%s,%i,%.2f" % \
                                      (gene, genes[gene]['start'] + offset,
                                      genes[gene]['end'] + offset, (i*3 + frame + offset),
                                      (i*3 + frame+2 + offset), ref_codon, codon, ref_aa, aa, coverage,
                                      aa_census.codon_frequency_for_amino_at(frame, i, aa, CONFIDENT, codon) / coverage*100)

                            if aa == ref_aa:
                                report += ",S"
                            else:
                                report += ",NS"

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
                                    codon_pathway = codon_pathway[:mutant_pos] + mutant_nt + codon_pathway[mutant_pos+1:]

                                    if Seq(codon_pathway).translate()[0] == ref_aa:
                                        s_count += 1
                                    else:
                                        ns_count += 1

                            report += ",%0.4f,%0.4f\n" % ( ns_count / len(codon_permutations[base_change_count-1]),
                                      s_count / len(codon_permutations[base_change_count-1]) )

    csv_file = open(output, "w")
    csv_file.write(report)
    csv_file.close()





