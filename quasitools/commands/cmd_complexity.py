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

import math

import quasitools.calculate as calculate
import quasitools.haplotype as haplotype

from quasitools.parsers.mapped_read_parser \
    import parse_pileup_from_fasta, parse_haplotypes_from_fasta

BASES = ['A', 'C', 'T', 'G']
GAP = '-'


@click.command(
    'complexity', short_help='Calculates various quasispecies complexity \
    measures.')
@click.argument('fasta', nargs=1,
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

    click.echo("\nStarting...")
    click.echo("\nCalculating the complexity from file: %s" % (fasta))

    complexity(ctx, fasta)

    click.echo("\nComplete!")


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

    pileup = parse_pileup_from_fasta(fasta)
    consensus = pileup.build_consensus()
    haplotypes = parse_haplotypes_from_fasta(fasta, consensus)

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
    report_sample_nucleotide_diversity(
        haplotypes, frequencies, distance_matrix)


def report_shannon_entropy(haplotypes, frequencies):
    """
    # ========================================================================

    REPORT SHANNON ENTROPY


    PURPOSE
    -------

    Reports the Shannon entropy of the haplotypes.


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        A list of Haplotype objects, for which to report the Shannon entropy.

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [NONE]


    POST
    ----

    The Shannon entropy, normalized to the number of haplotypes, and the
    Shannon entropy, normalized to the total number of clones, will be
    reported to output.

    # ========================================================================
    """

    Hs = calculate.shannon_entropy(frequencies)
    H = len(haplotypes)
    N = haplotype.calculate_total_clones(haplotypes)

    Hsn = float(Hs) / float(math.log(N))   # entropy normalized to N
    Hsh = float(Hs) / float(math.log(H))   # entropy normalized to H

    click.echo("Shannon Entropy (Hs) : " + str(Hs))
    click.echo("Shannon Entropy, normalized to log(N) (Hsn) : " + str(Hsn))
    click.echo("Shannon Entropy, normalized to log(H) (Hsh) : " + str(Hsh))


def report_minimum_mutation_frequency(pileup, haplotypes):
    """
    # ========================================================================

    REPORT MINIMUM MUTATION FREQUENCY


    PURPOSE
    -------

    Reports the minimum mutation frequency of the haplotypes.


    INPUT
    -----

    [PILEUP] [pileup]
        A Pileup object, which represents the pileup of aligned reads.

    [HAPLOTYPE LIST] [haplotypes]
        A list of Haplotype objects.


    RETURN
    ------

    [NONE]


    POST
    ----

    The minimum mutation frequency will be reported to output.

    # ========================================================================
    """

    M = pileup.count_unique_mutations()
    N = haplotype.calculate_total_clones(haplotypes)
    a = len(haplotypes[0].sequence)

    minimum_mutation_frequency = calculate.minimum_mutation_frequency(M, N, a)

    click.echo(
        "Minimum Mutation Frequency (Mf min) : "
        + str(minimum_mutation_frequency))


def report_mutation_frequency(distance_matrix):
    """
    # ========================================================================

    REPORT MUTATION FREQUENCY


    PURPOSE
    -------

    Reports the mutation frequency of the haplotypes.


    INPUT
    -----

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [NONE]


    POST
    ----

    The mutation frequency will be reported to output.

    # ========================================================================
    """

    H = len(distance_matrix)
    D = distance_matrix

    mutation_frequency = calculate.mutation_frequency(H, D)

    click.echo("Mutation Frequency (Mfe) : " + str(mutation_frequency))


def report_maximum_mutation_frequency(counts, frequencies, distance_matrix):
    """
    # ========================================================================

    REPORT MAXMIMUM MUTATION FREQUENCY


    PURPOSE
    -------

    Reports the maximum mutation frequency of the haplotypes.


    INPUT
    -----

    [INT LIST] [counts]
        A sorted list of haplotype counts, from the counts of the most
        abundant to the counts of the least abundant haplotype.

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)

    RETURN
    ------

    [NONE]


    POST
    ----

    The maximum mutation frequency will be reported to output.

    # ========================================================================
    """

    H = len(counts)
    F = frequencies
    D = distance_matrix

    maximum_mutation_frequency = calculate.maximum_mutation_frequency(H, F, D)

    click.echo(
        "Maximum Mutation Frequency (Mf max) : "
        + str(maximum_mutation_frequency))


def report_sample_nucleotide_diversity_entity(frequencies, distance_matrix):
    """
    # ========================================================================

    REPORT SAMPLE NUCLEOTIDE DIVERSITY (ENTITY LEVEL)


    PURPOSE
    -------

    Reports the entity-level sample nucleotide diversity.


    INPUT
    -----

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [NONE]


    POST
    ----

    The entity-level sample nucleotide diversity will be reported to output.

    # ========================================================================
    """

    H = len(frequencies)
    D = distance_matrix

    sample_nucleotide_diversity = \
        calculate.sample_nucleotide_diversity_entity(H, D)

    click.echo(
        "Sample Nucleotide Diversity, Entity Level (^PIe) : "
        + str(sample_nucleotide_diversity))


def report_population_nucleotide_diversity(frequencies, distance_matrix):
    """
    # ========================================================================

    REPORT POPULATION NUCLEOTIDE DIVERSITY


    PURPOSE
    -------

    Reports the population nucleotide diversity.


    INPUT
    -----

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [NONE]


    POST
    ----

    The population sample nucleotide diversity will be reported to output.

    # ========================================================================
    """

    H = len(frequencies)
    P = frequencies
    D = distance_matrix

    population_nucleotide_diversity = \
        calculate.population_nucleotide_diversity(H, P, D)

    click.echo(
        "Population Nucleotide Diversity (PI) : "
        + str(population_nucleotide_diversity))


def report_sample_nucleotide_diversity(
        haplotypes, frequencies, distance_matrix):
    """
    # ========================================================================

    REPORT SAMPLE NUCLEOTIDE DIVERSITY


    PURPOSE
    -------

    Reports the sample nucleotide diversity.


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        A list of Haplotype objects.

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [NONE]


    POST
    ----

    The sample nucleotide diversity will be reported to output.

    # ========================================================================
    """

    N = haplotype.calculate_total_clones(haplotypes)
    H = len(frequencies)
    P = frequencies
    D = distance_matrix

    sample_nucleotide_diversity = \
        calculate.sample_nucleotide_diversity(N, H, P, D)

    click.echo(
        "Sample Nucleotide Diversity (^PI) : "
        + str(sample_nucleotide_diversity))


def report_simpson_index(frequencies):
    """
    # ========================================================================

    REPORT SIMPSON INDEX


    PURPOSE
    -------

    Reports the simpson index.


    INPUT
    -----

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [NONE]


    POST
    ----

    The Simpson index will be reported to output.

    # ========================================================================
    """

    H = len(frequencies)
    P = frequencies

    simpson_index = calculate.simpson_index(H, P)

    click.echo("Simpson Index (Hsi) : " + str(simpson_index))


def report_gini_simpson_index(frequencies):
    """
    # ========================================================================

    REPORT GINI-SIMPSON INDEX


    PURPOSE
    -------

    Reports the Gini-Simpson index.


    INPUT
    -----

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [NONE]


    POST
    ----

    The Gini-Simpson index will be reported to output.

    # ========================================================================
    """

    H = len(frequencies)
    P = frequencies

    gini_simpson_index = calculate.gini_simpson_index(H, P)

    click.echo("Gini-Simpson Index (Hgs) : " + str(gini_simpson_index))


def report_hill_numbers(frequencies):
    """
    # ========================================================================

    REPORT HILL NUMBERS


    PURPOSE
    -------

    Reports the 0, 1, 2, and 3 Hill numbers.


    INPUT
    -----

    [FLOAT LIST] [frequencies]
        A list of (relative) frequencies of the Haplotypes.


    RETURN
    ------

    [NONE]


    POST
    ----

    The 0, 1, 2, and 3 Hill numbers will be reported to output.

    # ========================================================================
    """

    H = len(frequencies)
    P = frequencies

    q0 = calculate.hill_number(H, P, 0)
    q1 = calculate.hill_number(H, P, 1)
    q2 = calculate.hill_number(H, P, 2)
    q3 = calculate.hill_number(H, P, 3)

    click.echo("Hill numbers")
    click.echo("  q = 0 : " + str(q0))
    click.echo("  q = 1 : " + str(q1))
    click.echo("  q = 2 : " + str(q2))
    click.echo("  q = 3 : " + str(q3))


def report_FAD(distance_matrix):
    """
    # ========================================================================

    REPORT FUNCTIONAL ATTRIBUTE DIVERSITY


    PURPOSE
    -------

    Reports the functional attribute diversity.


    INPUT
    -----

    [2D ARRAY] [distance_matrix]
        A two dimensional array, representing the distance matrix of distances
        between the haplotypes.

        This is expected to be calculated in a similar manner as:
            haplotype.build_distiance_matrix(haplotypes)


    RETURN
    ------

    [NONE]


    POST
    ----

    The functional attribute diversity will be reported to output.

    # ========================================================================
    """

    H = len(distance_matrix)
    D = distance_matrix

    FAD = calculate.FAD(H, D)

    click.echo("Functional Attribute Diversity (FAD) : " + str(FAD))


def report_number_of_haplotypes(haplotypes):
    """
    # ========================================================================

    REPORT NUMBER OF HAPLOTYPES


    PURPOSE
    -------

    Reports the number of (unique) haplotypes.


    INPUT
    -----

    [HAPLOTYPE LIST] [haplotypes]
        A list of Haplotype objects.


    RETURN
    ------

    [NONE]


    POST
    ----

    The number of (unique) haplotypes will be reported to output.

    # ========================================================================
    """

    H = len(haplotypes)

    click.echo("Number of haplotypes (H) : " + str(H))


def report_number_of_mutations(pileup):
    """
    # ========================================================================

    REPORT NUMBER OF MUTATIONS


    PURPOSE
    -------

    Reports the number of mutations.


    INPUT
    -----

    [PILEUP] [pileup]
        A Pileup object, which represents the pileup of aligned reads.


    RETURN
    ------

    [NONE]


    POST
    ----

    The number of mutations will be reported to output.

    # ========================================================================
    """

    M = pileup.count_unique_mutations()

    click.echo("Number of unique mutations (M) : " + str(M))


def report_number_of_polymorphic_sites(pileup):
    """
    # ========================================================================

    REPORT NUMBER OF POLYMORPHIC SITES


    PURPOSE
    -------

    Reports the number of polymorphic sites.


    INPUT
    -----

    [PILEUP] [pileup]
        A Pileup object, which represents the pileup of aligned reads.


    RETURN
    ------

    [NONE]


    POST
    ----

    The number of polymorphic sites will be reported to output.

    # ========================================================================
    """

    P = pileup.count_polymorphic_sites()

    click.echo("Number of polymorphic sites (P) : " + str(P))
