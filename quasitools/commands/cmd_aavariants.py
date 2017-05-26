"""
Copyright Government of Canada 2017

Written by: Cole Peters, Eric Chubaty, National Microbiology Laboratory, Public Health Agency of Canada

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
from quasitools.aa_census import AACensus, CONFIDENT
from quasitools.aa_variant import AAVariantCollection
from quasitools.mutations import MutationDB
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.genes_file_parser import parse_genes_file


@click.command('aavariants',
               short_help='Identifies amino acid mutations.')
@click.argument('bam', required=True, type=click.Path(exists=True))
@click.argument('reference', required=True, type=click.Path(exists=True))
@click.argument('variants', required=True, type=click.Path(exists=True))
@click.argument('genes_file', required=True, type=click.Path(exists=True))
@click.argument('mutation_db', required=False, type=click.Path(exists=True))
@click.option('-f', '--min_freq', default=0.01,
              help='the minimum required frequency.')
@pass_context
def cli(ctx, bam, reference, variants, genes_file, min_freq, mutation_db):
    rs = parse_references_from_fasta(reference)

    mapped_read_collection_arr = []
    for r in rs:
        # Create a MappedReadCollection object
        mapped_read_collection_arr.append(parse_mapped_reads_from_bam(r, bam))

    # Mask the unconfident differences
    for mrc in mapped_read_collection_arr:
        mrc.mask_unconfident_differences(variants)

    # Parse the genes from the gene file
    genes = parse_genes_file(genes_file, rs[0].name)

    # Determine which frames our genes are in
    frames = set()

    for gene in genes:
        frames.add(genes[gene]['frame'])

    # Create an AACensus object
    aa_census = AACensus(reference, mapped_read_collection_arr, genes, frames)

    # Build the mutation database
    if mutation_db is not None:
        mutation_db = MutationDB(mutation_db, genes)

    # Create AAVar collection and print the hmcf file
    aa_vars = AAVariantCollection.from_aacensus(
        aa_census, next(iter(frames)))

    # Filter for mutant frequency
    aa_vars.filter('mf0.01', 'freq<0.01', True)
    click.echo(aa_vars.to_hmcf_file(CONFIDENT, mutation_db))
