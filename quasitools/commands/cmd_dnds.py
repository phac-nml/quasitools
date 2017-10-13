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
from quasitools.parsers.reference_parser import parse_references_from_fasta
from quasitools.parsers.codon_variant_file_parser \
    import parse_codon_variants


@click.command('dnds', short_help='Calculate the dn/ds '
               'value for each region in a bed file.')
@click.argument('csv', required=True, type=click.Path(exists=True))
@click.argument('reference', required=True, type=click.Path(exists=True))
@click.argument('offset', required=True, type=int)
@click.pass_context
def cli(ctx, csv, reference, offset):
    rs = parse_references_from_fasta(reference)
    ref_seq = rs[0].seq

    codon_variants = parse_codon_variants(csv, rs)

    click.echo(codon_variants.report_dnds_values(ref_seq, offset))
