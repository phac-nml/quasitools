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
import pysam
import os
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
@click.option('-i', '--id',
              help='specify default FASTA sequence identifier to be used ' 
              'for sequences without an RG tag.')
@click.option('-o', '--output', type=click.File('w'))
@pass_context
def cli(ctx, bam, reference, percentage, id, output):
    rs = parse_references_from_fasta(reference)
    bam_header = pysam.Samfile(bam, "rb").header

    if id:
        fasta_id = id
    else:
        fasta_id = os.path.basename(bam).split('.')[0]

    for r in rs:
        mrc = parse_mapped_reads_from_bam(r, bam)

        conseq = mrc.to_consensus(percentage)

        if hasattr(bam_header, 'RG'):
            fasta_id = bam_header['RG']

        if output:
            output.write('>{0}_{1}_{2}\n{3}'.format(fasta_id, percentage, r.name,
                                                    conseq))
            output.close()
        else:
            click.echo('>{0}_{1}_{2}\n{3}'.format(fasta_id, percentage, r.name,
                                                  conseq))
