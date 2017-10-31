"""
Copyright Government of Canada 2015

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
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam


@click.command('consensus',
               short_help='Generate a consensus sequence from a BAM file.')
@click.argument('bam', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('reference', required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-p', '--percentage', default=100,
              help='percentage to include base in mixture.')
@pass_context
def cli(ctx, bam, reference, percentage):
    rs = parse_references_from_fasta(reference)

    for r in rs:
        mrc = parse_mapped_reads_from_bam(r, bam)

        conseq = mrc.to_consensus(percentage)

        click.echo('>{0}_{1}_{2}\n{3}'.format('blah', percentage, r.name,
                                              conseq))
