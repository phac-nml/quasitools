"""
Copyright Government of Canada 2017

Written by: Cole Peters, National Microbiology Laboratory,
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
from quasitools.aa_census import AACensus
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.genes_file_parser import parse_genes_file


@click.command('aa_census',
               short_help='Builds an amino acid census and returns its '
               'coverage.')
@click.argument('bam', required=True, type=click.Path(exists=True))
@click.argument('reference', required=True, type=click.Path(exists=True))
@click.argument('genes_file', required=True, type=click.Path(exists=True))
@pass_context
def cli(ctx, bam, reference, genes_file):
    rs = parse_references_from_fasta(reference)

    mapped_read_collection_arr = []
    for r in rs:
        # Create a MappedReadCollection object
        mapped_read_collection_arr.append(parse_mapped_reads_from_bam(r, bam))

    # Parse the genes from the gene file
    genes = parse_genes_file(genes_file, rs[0].name)

    # Determine which frames our genes are in
    frames = set()

    for gene in genes:
        frames.add(genes[gene]["frame"])

    # Create an AACensus object
    aa_census = AACensus(reference, mapped_read_collection_arr, genes, frames)

    click.echo(aa_census.coverage(next(iter(frames))))
