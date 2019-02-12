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
import math
import quasitools.calculate as calculate
import quasitools.haplotype as haplotype
from quasitools.parsers.reference_parser import parse_references_from_fasta

#Remove this
import time

from quasitools.parsers.mapped_read_parser \
        import parse_pileup_from_bam, parse_haplotypes_from_fasta, parse_haplotypes_called


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
    start_time = time.time()
    click.echo("\nStarting...")

    click.echo("Using file %s as reference" % reference)
    click.echo("Reading input from file(s)  %s" % bam)
  
    complexity(ctx, reference, bam)

    click.echo("\nComplete!")
    print("--- %s seconds ---" % (time.time() - start_time))


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
    pileup = parse_pileup_from_bam(references, bam)
    consensus = pileup.build_consensus_from_range(1,50)
    haplotypes = parse_haplotypes_called(references, reference, bam, consensus,50,50)
    
    count = 1
    for i in haplotypes:
        #print((len(i)))
        count+=1

    #print(len(consensus))

     

    distance_matrix = haplotype.build_distiance_matrix(haplotypes)
    counts = haplotype.build_counts(haplotypes)
    frequencies = haplotype.build_frequencies(haplotypes)
    


    # Will create an array of dictionaries to hold each measurment
    measurements = [] 
    
    '''
    Set the Incidence - Entity Level
    '''

    number_of_haplotypes = get_number_of_haplotypes(haplotypes)
    measurements.append({'Number of haplotypes': number_of_haplotypes})
    
    number_of_polymorphic_sites = get_number_of_polymorphic_sites(pileup)
    measurements.append({'Number of polymorphic sites': number_of_polymorphic_sites})
    
    number_of_mutations = get_number_of_mutations(pileup)
    measurements.append({'Number of mutation': number_of_mutations})

    '''
    Set the Abundance - Molecular Level
    '''
    
    measurements = get_shannon_entropy(haplotypes, frequencies, measurements)
    
    simpson_index = get_simpson_index(frequencies)
    measurements.append({'Simpson Index': simpson_index})
    
    gini_simpson_index = get_gini_simpson_index(frequencies)
    measurements.append({'Gini Simpson Index': gini_simpson_index})
    

    measurements = get_hill_numbers(measurements, frequencies)
    
    '''
    Functional,  Indidence - Entity Level
    '''

    minimum_mutation = get_minimum_mutation_frequency(pileup, haplotypes)
    measurements.append({'Minimum mutation frequency': minimum_mutation}) 
    
    measurements = get_mutation_frequency(distance_matrix,measurements)
    measurements = get_FAD(distance_matrix,measurements)
    measurements = get_sample_nucleotide_diversity_entity(distance_matrix, frequencies, measurements)

    '''
    Functional, Abundance - Molecular Level
    '''
    measurements = get_maximum_mutation_frequency(counts, distance_matrix, frequencies, measurements)
    measurements = get_population_nucleotide_diversity(distance_matrix, frequencies, measurements)
    
    '''
    Other
    '''
    measurements = get_sample_nucleotide_diversity(distance_matrix, frequencies, haplotypes, measurements)

    #for x in range(len(measurements)):
        #print(measurements[x])
    for measurement in measurements:
        print(str(measurement).strip("{}").replace("'",""))

def get_sample_nucleotide_diversity(distance_matrix,frequencies, haplotypes, measurements):
    N = haplotype.calculate_total_clones(haplotypes)
    H = len(frequencies)
    P = frequencies
    D = distance_matrix

    SND = calculate.sample_nucleotide_diversity(N,H,P,D)
    measurements.append({'Sample Nucleotide Diversity (^PI)': SND})

    return measurements

def get_population_nucleotide_diversity(distance_matrix, frequencies, measurements):
    
    H = len(frequencies)
    P = frequencies
    D = distance_matrix

    PND = calculate.population_nucleotide_diversity(H,P,D)
    measurements.append({'Population Nucleotide Diversity (PI)': PND})

    return measurements
def get_maximum_mutation_frequency(counts, distance_matrix, frequencies, measurements):
    
    H = len(counts)
    F = frequencies
    D = distance_matrix

    maximum_mutation_frequency = calculate.maximum_mutation_frequency(H,F,D)
    measurements.append({'Maximum mutation frequency (Mf max)': maximum_mutation_frequency})
    
    return measurements


def get_sample_nucleotide_diversity_entity(distance_matrix, frequencies, measurements):
    H = len(frequencies)
    D = distance_matrix

    PND = calculate.sample_nucleotide_diversity_entity(H,D)
    measurements.append({'Population Nucleotide Diversity': PND})
    
    return measurements
def get_FAD(distance_matrix, measurements):
    
    H = len(distance_matrix)
    D = distance_matrix

    FAD = calculate.FAD(H,D)
    
    measurements.append({'Functional attribute diversity': FAD})
    
    return measurements

def get_mutation_frequency(distance_matrix,measurements):

    H = len(distance_matrix)
    D = distance_matrix
    mutation_frequency = calculate.mutation_frequency(H,D)

    measurements.append({'Mutation frequency': mutation_frequency})
    return measurements

def get_minimum_mutation_frequency(pileup, haplotypes):

    M = pileup.count_unique_mutations()
    N = haplotype.calculate_total_clones(haplotypes)
    a = len(haplotypes[0].sequence)

    minimum_mutation_frequency = calculate.minimum_mutation_frequency(M,N,a)

    return minimum_mutation_frequency
    

def get_number_of_haplotypes(haplotypes): 
   
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
    
    haplotype_num = len(haplotypes)


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

def get_shannon_entropy(haplotypes, frequencies, measurements):
    """""
    #========================================================================
    
    GET SHANNON ENTROPY

    PURPOSE
    -------

    Reports the shannon Entropy of the haplotpyes

    INPUT
    -------
    [HAPLOTYPES] [haplotypes]
    [FLOAT LIST] [frequencies]
    RETURN
    -------
    {DICTIONARY}[shannon_entropy_numbers]

    COMMENTS
    -------
    [N/A]
    #========================================================================
    """


    Hs = calculate.shannon_entropy(frequencies)
    H = len(haplotypes)
    N = haplotype.calculate_total_clones(haplotypes)
    
    if float(math.log(N)) != 0:
        Hsn = float (Hs) / float (math.log(N)) # entropy localized to N
        Hsh = float (Hs) / float (math.log(H)) # entropy localized to H
    else:
        Hsn = 0
        Hsh = 0
    
    measurements.append({'Shannon Entropy (Hs)': Hs})
    measurements.append({'Shannon Entropy Localized to N (Hsn)': Hsn})
    measurements.append({'Shannon Entropy Localized to H (Hsh)': Hsh})
   
    return measurements

def get_simpson_index(frequencies):
    """""
    #========================================================================
    REPORT SIMPSON INDEX

    PURPOSE
    -------
    
    Reports the Simpson Index

    INPUT
    -------
    [FLOAT LIST] frequencies
        a list of (relative) frequencies) of the haplotpyes.

    RETURN
    -------
    [STR] [simpson_index]

    COMMENTS
    -------
    [N/A]

    #========================================================================
    """

    H = len(frequencies)
    P = frequencies

    simpson_index = calculate.simpson_index(H,P)
    
    return simpson_index

def get_gini_simpson_index(frequencies):
    """""
    #========================================================================
    
    GET GINI SIMPSON INDEX
    
    PURPOSE
    -------
    Reports the Gini-Simpson Index

    INPUT
    -------

    [FLOAT LIST] [frequencies]
        a list of relative frequencies of the haplotypes

    RETURN
    -------
    [INT] gini_simpson_index

    COMMENTS
    -------

    #========================================================================
    """

    H = len(frequencies)
    P =  frequencies
    gini_simpson_index = calculate.gini_simpson_index(H,P)

    return gini_simpson_index

def get_hill_numbers(measurements, frequencies):
    """""
    #========================================================================
    
    GET HILL NUMBERS

    PURPOSE
    -------
    
    Reports the 0, 1, 2, and 3 hill numbers.

    INPUT
    -------
    [FLOAT LIST] [frequencies]
        a list of relative frequencies of the haplotypes.

    RETURN
    -------
    [Dictionary List] [measurements]

    COMMENTS
    -------

    #========================================================================
    """

    # H = len(frequencie)
    P = frequencies
    H = 50
    #print("Hill debug")
    #print(H)
    #print(P)

    hill0 = calculate.hill_number(H,P,0)
    hill1 = calculate.hill_number(H,P,1)
    hill2 = calculate.hill_number(H,P,2)
    hill3 = calculate.hill_number(H,P,3)
    
    measurements.append({'Hill Number 0': hill0})
    measurements.append({'Hill Number 1': hill1})
    measurements.append({'Hill Number 2': hill2})
    measurements.append({'Hill Numbwe 3': hill2})

    return measurements

