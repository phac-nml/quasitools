"""
Copyright Government of Canada 2018

Written by: Eric Marinier and Ahmed Kidwai, Public Health Agency of Canada

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

import quasitools.calculate as calculate
import quasitools.haplotype as haplotype
from quasitools.parsers.reference_parser import parse_references_from_fasta

from quasitools.parsers.mapped_read_parser \
    import parse_pileup_from_bam, parse_haplotypes_from_fasta

BASES = ['A', 'C', 'T', 'G']
GAP = '-'


@click.command(
    'complexity', short_help='Calculates various quasispecies complexity \
    measures.')
@click.argument('reference', nargs=1, required=True,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('bam', nargs=1,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument('fasta', nargs=1,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.pass_context
def cli(ctx, reference, bam, fasta):
    """

    Reports the complexity of a quasispecies using several measures outlined
    in the following work:

    Gregori, Josep, et al. "Viral quasispecies complexity measures."
    Virology 493 (2016): 227-237.

    """

    click.echo("\nStarting...")

    click.echo("Using file %s as reference" % reference)
    click.echo("Reading input from file(s)  %s" % bam)
    click.echo("Reading converted bam to fasta as %s" % fasta)

    complexity(ctx, reference, bam, fasta)

    click.echo("\nComplete!")


def complexity(ctx, reference, bam, fasta):
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
    pileup = parse_pileup_from_bam(references, bam)


    # build a consensus
    consensus_from_the_pileup = pileup.build_consensus
    haplotypes = parse_haplotypes_from_fasta(fasta, consensus_from_the_pileup)

    # Will create an array of dictionaries to hold each measurment
    measurements = [] 
    
    '''
    Set the Incidence - Entity Level
    '''

    number_of_haplotypes = get_number_of_haplotypes()
    measurements.append({'Number of haplotypes': number_of_haplotypes})
    
    number_of_polymorphic_sites = get_number_of_polymorphic_sites(pileup)
    measurements.append({'Number of polymorphic sites': number_of_polymorphic_sites})
    
    number_of_mutations = get_number_of_mutations(pileup)
    measurements.append({'Number of mutation': number_of_mutations})

    '''
    Set the Abundance - Molecular Level
    '''
    

    for x in range(len(measurements)):
       print(measurements[x])

        
def get_number_of_haplotypes(): 
   
    """""
    #========================================================================
    
    GET NUMBER OF HAPLOTYPES


    PURPOSE
    -------

    Reports the number of unique haplotypes


    INPUT
    -------

    [NONE]


    RETURN
    -------

    [UNIQUE NUMBER OF HAPLOTYPES]


    # ========================================================================
    """
    # The number of unique haplotypes is theoretically 4^1 i.e 4
    # in the future this method will prob be modified.
    
    haplotype_num = 4

    return haplotype_num


def get_number_of_polymorphic_sites(pileup):
    """""
    #========================================================================
    
    GET NUMBER OF POLYMORPHIC SITES
    

    PURPOSE
    -------
    
    Reports the number of polymorphic sites
    
    INPUT
    -------
    
    [PILEUP] [pileup]
        A pileup object that represents the pileup of aligned reads

    RETURN
    -------
    
    [INT] [number_of_polymorphic_sites]

    COMMENTS
    -------
    [N/A]

    #========================================================================
    """
    number_of_polymorphic_sites = pileup.count_polymorphic_sites()
    return number_of_polymorphic_sites

def get_number_of_mutations(pileup):
    """""
    #========================================================================
    
    GET NUMBER OF MUTATIONS
    PURPOSE
    -------
    
    Gets the number of polymorphic sites

    INPUT
    -------
    
    [PILEUP] [pileup]
        A pileup object that represents the pileup of aligned reads

    RETURN
    -------

    [INT] [number_of_mutations]

    COMMENTS
    -------
    [N/A]

    #========================================================================
    """

    number_of_mutations = pileup.count_unique_mutations()
    return number_of_mutations

