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

import click
from quasitools.cli import pass_context
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.genes_file_parser import parse_genes_file
from quasitools.nt_variant import NTVariantCollection
from quasitools.aa_census import AACensus
from quasitools.codon_variant import CodonVariantCollection

@click.command('mutanttypes', short_help='Identify the number of non-synonymous and synonymous mutations.')
@click.argument('bam', required=True, type=click.Path(exists=True))
@click.argument('reference', required=True, type=click.Path(exists=True))
@click.argument('offset', required=True, type=float)
@click.argument('genes_file', required=True, type=click.Path(exists=True))
@click.argument('output', required=True)
@click.option('-e', '--error_rate', default=0.01, help='estimated sequencing error rate.')

@pass_context
def cli(ctx, bam, reference, offset, genes_file, output, error_rate):
    click.echo("Running mutanttypes command...")

    rs = parse_references_from_fasta(reference)
    mapped_read_collection_arr = []

    # Create a MappedReadCollection object
    for r in rs:
        mapped_read_collection_arr.append(parse_mapped_reads_from_bam(r, bam))

    variants_obj = NTVariantCollection.from_mapped_read_collections(error_rate, rs, 
                    *mapped_read_collection_arr)
    variants_obj.filter('q30', 'QUAL<30', True)
    variants_obj.filter('ac5', 'AC<5', True)
    variants_obj.filter('dp100', 'DP<100', True)

    click.echo("Creating variants vcf file...")
    vcf_filename = "{}.vcf".format(output.split('.')[0])
    vcf_file = open(vcf_filename, "w")
    vcf_file.write(variants_obj.to_vcf_file())
    vcf_file.close()

    # Mask the unconfident differences
    for mrc in mapped_read_collection_arr:
        mrc.mask_unconfident_differences_from_file(vcf_filename)

    # Parse the genes from the gene file
    genes = parse_genes_file(genes_file, rs[0].name)

    # Determine which frames our genes are in
    frames = set()

    for gene in genes:
        frames.add(genes[gene]['frame'])

    click.echo("Creating AACensus object...")
    aa_census = AACensus(reference, mapped_read_collection_arr, genes, frames)

    click.echo("Creating CodonVariantCollection")
    codon_variants = CodonVariantCollection.from_aacensus(aa_census, next(iter(frames)))

    report = codon_variants.to_csv_file(offset)

    csv_file = open(output, "w")
    csv_file.write(report)
    csv_file.close()

    click.echo(report)

