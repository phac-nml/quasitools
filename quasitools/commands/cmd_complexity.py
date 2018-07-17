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

__version__ = '0.1.0'

import click
import time

import os
import sys

import calculate
import pileup
import haplotype

from quasitools.parsers.mapped_read_parser import parse_pileup_from_fasta

BASES = ['A', 'C', 'T', 'G']
GAP = '-'

"""
# ============================================================================

CLICK

# ============================================================================
"""
@click.command('complexity', short_help='Calculates various quasispecies \
        complexity measures.')

@click.argument('fasta', nargs=-1,
                type=click.Path(exists=True, file_okay=True, dir_okay=False))

@click.pass_context
def cli(ctx, fasta):
    """
    # ========================================================================

    CLICK PARSER / COMPLEXITY


    PURPOSE
    -------

    Parses the command line parameters and initiates computation of
    complexity.


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

    click.echo("\n")
    click.echo("Starting...")

    click.echo("Calculating the complexity from file: %s" % (fasta))

    complexity(ctx, fasta)

    click.echo("Complete!")


def complexity(ctx, fasta):
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

    pileup = parse_pileup_from_fasta(reads_location)
    consensus = pileup.build_consensus()
    haplotypes = haplotype.build_from_reads(reads_location, consensus)

    distance_matrix = haplotype.build_distiance_matrix(haplotypes)
    counts = haplotype.build_counts(haplotypes)
    frequencies = haplotype.build_frequencies(haplotypes)

    click.echo("\n")
    click.echo("Incidence - Entity Level")
    click.echo("-------------------------")
    report_number_of_haplotypes(haplotypes)
    report_number_of_polymorphic_sites(pileup)
    report_number_of_mutations(pileup)

    click.echo("\n")
    click.echo("Abundance - Molecular Level")
    click.echo("---------------------------")
    report_shannon_entropy(haplotypes, frequencies)
    report_simpson_index(frequencies)
    report_gini_simpson_index(frequencies)
    report_hill_numbers(frequencies)

    click.echo("\n")
    click.echo("Functional, Indidence - Entity Level")
    click.echo("-------------------------------------")
    report_minimum_mutation_frequency(pileup, haplotypes)
    report_mutation_frequency(distance_matrix)
    report_FAD(distance_matrix)
    report_sample_nucleotide_diversity_entity(frequencies, distance_matrix)

    click.echo("\n")
    click.echo("Functional, Abundance - Molecular Level")
    click.echo("----------------------------------------")
    report_maximum_mutation_frequency(counts, frequencies, distance_matrix)
    report_population_nucleotide_diversity(frequencies, distance_matrix)

    click.echo("\n")
    click.echo("Other")
    click.echo("------")
    report_sample_nucleotide_diversity(haplotypes, frequencies, distance_matrix)

"""
# =============================================================================
# =============================================================================
"""
def report_shannon_entropy(haplotypes, frequencies):

    Hs = calculate.shannon_entropy(frequencies)
    H = len(haplotypes)
    N = haplotype.calculate_total_clones(haplotypes)

    Hsn = float(Hs) / float(math.log(N))   # entropy normalized to N
    Hsh = float(Hs) / float(math.log(H))   # entropy normalized to H

    print("Shannon Entropy (Hs) : " + str(Hs))
    print("Shannon Entropy, normalized to log(N) (Hsn) : " + str(Hsn))
    print("Shannon Entropy, normalized to log(H) (Hsh) : " + str(Hsh))

"""
# =============================================================================
# =============================================================================
"""
def report_minimum_mutation_frequency(pileup, haplotypes):

    M = pileup.count_unique_mutations()
    N = haplotype.calculate_total_clones(haplotypes)
    a = len(haplotypes[0].sequence)

    minimum_mutation_frequency = calculate.minimum_mutation_frequency(M, N, a)

    print("Minimum Mutation Frequency (Mf min) : " + str(minimum_mutation_frequency))

"""
# =============================================================================
# =============================================================================
"""
def report_mutation_frequency(distance_matrix):

    H = len(distance_matrix)
    D = distance_matrix

    mutation_frequency = calculate.mutation_frequency(H, D)

    print("Mutation Frequency (Mfe) : " + str(mutation_frequency))


"""
# =============================================================================
# =============================================================================
"""
def report_maximum_mutation_frequency(counts, frequencies, distance_matrix):

    H = len(counts)
    F = frequencies
    D = distance_matrix

    maximum_mutation_frequency = calculate.maximum_mutation_frequency(H, F, D)

    print("Maximum Mutation Frequency (Mf max) : " + str(maximum_mutation_frequency))

"""
# =============================================================================
# =============================================================================
"""
def report_sample_nucleotide_diversity_entity(frequencies, distance_matrix):

    H = len(frequencies)
    D = distance_matrix

    sample_nucleotide_diversity = calculate.sample_nucleotide_diversity_entity(H, D)

    print("Sample Nucleotide Diversity, Entity Level (^PIe) : " + str(sample_nucleotide_diversity))


"""
# =============================================================================
# =============================================================================
"""
def report_population_nucleotide_diversity(frequencies, distance_matrix):

    H = len(frequencies)
    P = frequencies
    D = distance_matrix

    population_nucleotide_diversity = calculate.population_nucleotide_diversity(H, P, D)

    print("Population Nucleotide Diversity (PI) : " + str(population_nucleotide_diversity))

"""
# =============================================================================
# =============================================================================
"""
def report_sample_nucleotide_diversity(haplotypes, frequencies, distance_matrix):

    N = haplotype.calculate_total_clones(haplotypes)
    H = len(frequencies)
    P = frequencies
    D = distance_matrix

    sample_nucleotide_diversity = calculate.sample_nucleotide_diversity(N, H, P, D)

    print("Sample Nucleotide Diversity (^PI) : " + str(sample_nucleotide_diversity))

"""
# =============================================================================
# =============================================================================
"""
def report_simpson_index(frequencies):

    H = len(frequencies)
    P = frequencies

    simpson_index = calculate.simpson_index(H, P)

    print("Simpson Index (Hsi) : " + str(simpson_index))

"""
# =============================================================================
# =============================================================================
"""
def report_gini_simpson_index(frequencies):

    H = len(frequencies)
    P = frequencies

    gini_simpson_index = calculate.gini_simpson_index(H, P)

    print("Gini-Simpson Index (Hgs) : " + str(gini_simpson_index))

"""
# =============================================================================
# =============================================================================
"""
def report_hill_numbers(frequencies):

    H = len(frequencies)
    P = frequencies

    q0 = calculate.hill_number(H, P, 0)
    q1 = calculate.hill_number(H, P, 1)
    q2 = calculate.hill_number(H, P, 2)
    q3 = calculate.hill_number(H, P, 3)

    print("Hill numbers")
    print("  q = 0 : " + str(q0))
    print("  q = 1 : " + str(q1))
    print("  q = 2 : " + str(q2))
    print("  q = 3 : " + str(q3))

"""
# =============================================================================
# =============================================================================
"""
def report_FAD(distance_matrix):

    H = len(distance_matrix)
    D = distance_matrix

    FAD = calculate.FAD(H, D)

    print("Functional Attribute Diversity (FAD) : " + str(FAD))

"""
# =============================================================================
# =============================================================================
"""
def report_number_of_haplotypes(haplotypes):

    H = len(haplotypes)

    print("Number of haplotypes (H) : " + str(H))

"""
# =============================================================================
# =============================================================================
"""
def report_number_of_mutations(pileup):

    M = pileup.count_unique_mutations()

    print("Number of unique mutations (M) : " + str(M))

"""
# =============================================================================
# =============================================================================
"""
def report_number_of_polymorphic_sites(pileup):

    P = pileup.count_polymorphic_sites()

    print("Number of polymorphic sites (P) : " + str(P))
