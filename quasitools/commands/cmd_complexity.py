"""
Copyright Government of Canada 2018

Written by: Eric Marinier, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

__version__ = '0.1.1'

import click

from quasitools.parsers.reference_parser import parse_references_from_fasta

from quasitools.parsers.mapped_read_parser \
    import parse_pileup_from_bam

BASES = ['A', 'C', 'T', 'G']
GAP = '-'


@click.command(
    'complexity', short_help='Calculates various quasispecies complexity \
    measures.')
@click.argument('reference', nargs=1, required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('bam', nargs=1,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.pass_context
def cli(ctx, reference, bam):
    """

    Reports the complexity of a quasispecies using several measures outlined
    in the following work:

    Gregori, Josep, et al. "Viral quasispecies complexity measures."
    Virology 493 (2016): 227-237.

    """

    click.echo("\nStarting...")

    click.echo("Using file %s as reference" % reference)
    click.echo("Reading input from file(s)  %s" % bam)

    complexity(ctx, reference, bam)

    click.echo("\nComplete!")


def complexity(ctx, reference, bam):
    """
    # ========================================================================

    COMPLEXITY


    PURPOSE
    -------

    Computes various complexity measures for the quasispecies.


    INPUT
    -----

    [(FASTA) FILE LOCATION] [fasta]
        The file location of an aligned FASTA file for which to calculate the
        complexity measures.


    RETURN
    ------

    [NONE]


    POST
    ----

    The complexity computation and reporting will be complete.

    # ========================================================================
    """

    references = parse_references_from_fasta(reference)
    # for x in bam:
    pileup = parse_pileup_from_bam(references, bam)
    # This will be removed just couldn't have an
    # unused variiable floating around or flake8 would get mad at me.
    pileup.count_polymorphic_sites()
