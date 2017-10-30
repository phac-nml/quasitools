"""
Copyright Government of Canada 2015-2017

Written by: Eric Enns, National Microbiology Laboratory,
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
from quasitools.nt_variant import NTVariantCollection
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.reference_parser import parse_references_from_fasta


@click.command('call_callntvar',
               short_help='Call nucleotide variants from a BAM file.')
@click.argument('bam', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('reference', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-e', '--error_rate', default=0.01,
              help='estimated sequencing error rate.')
@pass_context
def cli(ctx, bam, reference, error_rate):
    rs = parse_references_from_fasta(reference)

    mapped_read_collection_arr = []
    for r in rs:
        # create MappedReadCollection object
        mapped_read_collection_arr.append(parse_mapped_reads_from_bam(r, bam))

    variants = NTVariantCollection.from_mapped_read_collections(
        error_rate, rs, *mapped_read_collection_arr)

    variants.filter('q30', 'QUAL<30', True)
    variants.filter('ac5', 'AC<5', True)
    variants.filter('dp100', 'DP<100', True)

    click.echo(variants.to_vcf_file())
