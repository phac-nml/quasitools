"""
Copyright Government of Canada 2015

Written by: Eric Enns, National Microbiology Laboratory, Public Health Agency of Canada

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
from quasitools.variants import Variants
from quasitools.reference import Reference

@click.command('call_callntvar', short_help='Call nucleotide variants from a BAM file.')
@click.argument('bam', required=True, type=click.Path(exists=True))
@click.argument('reference', required=True, type=click.Path(exists=True))
@pass_context
def cli(ctx, bam, reference):
    r = Reference.from_fasta(reference)
    variants = Variants.from_bam(r, 65, 75, bam)

    variants.filter('q30', 'QUAL<30', True)
    variants.filter('ac5', 'AC<5', True)
    variants.filter('dp100', 'DP<100', True)

    click.echo(str(variants))
