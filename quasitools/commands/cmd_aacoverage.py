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
from quasitools.parsers.genes_file_parser import parse_BED4_file


@click.command('aa_coverage',
               short_help='Builds an amino acid census and returns its '
               'coverage.')
@click.argument('bam', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('reference', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('bed4_file', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-o', '--output', type=click.File('w'))
@pass_context
def cli(ctx, bam, reference, bed4_file, output):
    """This script builds an amino acid census and returns its coverage.
    The BAM alignment file corresponds to a pileup of sequences aligned to
    the REFERENCE. A BAM index file (.bai) must also be present and, except
    for the extension, have the same name as the BAM file. The REFERENCE must
    be in FASTA format. The BED4_FILE must be a BED file with at least 4
    columns and specify the gene locations within the REFERENCE.

    The output is in CSV format."""

    rs = parse_references_from_fasta(reference)

    mapped_read_collection_arr = []
    for r in rs:
        # Create a MappedReadCollection object
        mapped_read_collection_arr.append(parse_mapped_reads_from_bam(r, bam))

    # Parse the genes from the gene file
    genes = parse_BED4_file(bed4_file, rs[0].name)

    # Determine which frames our genes are in
    frames = set()

    for gene in genes:
        frames.add(genes[gene]["frame"])

    # Create an AACensus object
    aa_census = AACensus(reference, mapped_read_collection_arr, genes, frames)

    if output:
        output.write(aa_census.coverage(frames))
        output.close()
    else:
        click.echo(aa_census.coverage(frames))
