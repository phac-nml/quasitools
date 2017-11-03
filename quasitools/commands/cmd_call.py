"""
Copyright Government of Canada 2017

Written by: Eric Enns, Cole Peters, Eric Chubaty, Camy Tran,
            National Microbiology Laboratory, Public Health Agency of Canada

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
from quasitools.nt_variant import NTVariantCollection
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.aa_census import AACensus, CONFIDENT
from quasitools.aa_variant import AAVariantCollection
from quasitools.mutations import MutationDB
from quasitools.parsers.genes_file_parser import parse_genes_file
from quasitools.parsers.nt_variant_file_parser \
    import parse_nt_variants_from_vcf
from quasitools.codon_variant import CodonVariantCollection


@click.group(invoke_without_command=False)
@click.pass_context
def cli(ctx):
    pass


@cli.command('ntvar', short_help='Call nucleotide variants from a BAM file.')
@click.argument('bam', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('reference', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-e', '--error_rate', default=0.01,
              help='estimated sequencing error rate.')
@click.option('-o', '--output', type=click.File('w'))
def ntvar(bam, reference, error_rate, output):
    rs = parse_references_from_fasta(reference)

    mapped_read_collection_arr = []
    for r in rs:
        # create MappedReadCollection object
        mapped_read_collection_arr.append(parse_mapped_reads_from_bam(r, bam))

    mapped_read_collection_arr = []
    for r in rs:
        # create MappedReadCollection object
        mapped_read_collection_arr.append(parse_mapped_reads_from_bam(r, bam))

    variants = NTVariantCollection.from_mapped_read_collections(
        error_rate, rs, *mapped_read_collection_arr)

    variants.filter('q30', 'QUAL<30', True)
    variants.filter('ac5', 'AC<5', True)
    variants.filter('dp100', 'DP<100', True)

    if output:
        output.write(variants.to_vcf_file())
        output.close()
    else:
        click.echo(variants.to_vcf_file())


@cli.command('aavar', short_help='Identifies amino acid mutations.')
@click.argument('bam', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('reference', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('variants', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('genes_file', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('mutation_db', required=False,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-f', '--min_freq', default=0.01,
              help='the minimum required frequency.')
@click.option('-o', '--output', type=click.File('w'))
def aavar(bam, reference, variants, genes_file, min_freq,
          mutation_db, output):
    rs = parse_references_from_fasta(reference)

    mapped_read_collection_arr = []
    for r in rs:
        # Create a MappedReadCollection object
        mapped_read_collection_arr.append(parse_mapped_reads_from_bam(r, bam))

    variants_obj = parse_nt_variants_from_vcf(variants, rs)

    # Mask the unconfident differences
    for mrc in mapped_read_collection_arr:
        mrc.mask_unconfident_differences(variants_obj)

    # Parse the genes from the gene file
    genes = parse_genes_file(genes_file, rs[0].name)

    # Determine which frames our genes are in
    frames = set()

    for gene in genes:
        frames.add(genes[gene]['frame'])

    # Create an AACensus object
    aa_census = AACensus(reference, mapped_read_collection_arr, genes, frames)

    # Create AAVar collection and print the hmcf file
    aa_vars = AAVariantCollection.from_aacensus(aa_census)

    # Filter for mutant frequency
    aa_vars.filter('mf0.01', 'freq<0.01', True)

    # Build the mutation database and update collection
    if mutation_db is not None:
        mutation_db = MutationDB(mutation_db, genes)
        aa_vars.apply_mutation_db(mutation_db)

    if output:
        output.write(aa_vars.to_hmcf_file(CONFIDENT))
    else:
        click.echo(aa_vars.to_hmcf_file(CONFIDENT))


@cli.command('codonvar', short_help='Identify the number of '
             'non-synonymous and synonymous mutations.')
@click.argument('bam', required=True, type=click.Path(exists=True,
                file_okay=True, dir_okay=False))
@click.argument('reference', required=True, type=click.Path(exists=True,
                file_okay=True, dir_okay=False))
@click.argument('offset', required=True, type=float)
@click.argument('genes_file', required=True, type=click.Path(exists=True,
                file_okay=True, dir_okay=False))
@click.option('-e', '--error_rate', default=0.01,
              help='estimated sequencing error rate.')
@click.option('-o', '--output', type=click.File('w'))
def codonvar(bam, reference, offset, genes_file, error_rate, output):
    rs = parse_references_from_fasta(reference)
    mapped_read_collection_arr = []

    # Create a MappedReadCollection object
    for r in rs:
        mapped_read_collection_arr.append(parse_mapped_reads_from_bam(r, bam))

    variants = NTVariantCollection.from_mapped_read_collections(
        error_rate, rs, *mapped_read_collection_arr)
    variants.filter('q30', 'QUAL<30', True)
    variants.filter('ac5', 'AC<5', True)
    variants.filter('dp100', 'DP<100', True)

    # Mask the unconfident differences
    for mrc in mapped_read_collection_arr:
        mrc.mask_unconfident_differences(variants)

    # Parse the genes from the gene file
    genes = parse_genes_file(genes_file, rs[0].name)

    # Determine which frames our genes are in
    frames = set()

    for gene in genes:
        frames.add(genes[gene]['frame'])

    aa_census = AACensus(reference, mapped_read_collection_arr, genes, frames)

    codon_variants = CodonVariantCollection.from_aacensus(aa_census)

    click.echo(codon_variants.to_csv_file(offset))

    variants.filter('dp100', 'DP<100', True)

    if output:
        output.write(variants.to_vcf_file())
        output.close()
    else:
        click.echo(variants.to_vcf_file())
