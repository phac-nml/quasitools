UNDEFINED = ""

# Measurement Positions
NUMBER_OF_HAPLOTYPES = 0
HAPLOTYPE_POPULATION = 1
NUMBER_OF_POLYMORPHIC_SITES = 2
NUMBER_OF_MUTATIONS = 3
SHANNON_ENTROPY_NUMBER = 4
SHANNON_ENTROPY_NUMBER_NORMALIZED_TO_N = 5
SHANNON_ENTROPY_NUMBER_NORMALIZED_TO_H = 6
SIMPSON_INDEX = 7
GINI_SIMPSON_INDEX = 8
HILL_NUMBER_0 = 9
HILL_NUMBER_1 = 10
HILL_NUMBER_2 = 11
HILL_NUMBER_3 = 12
MINIMUM_MUTATION_FREQUENCY = 13
MUTATION_FREQUENCY = 14
FUNCTIONAL_ATTRIBUTE_DIVERSITY = 15
SAMPLE_NUCLEOTIDE_DIVERSITY_ENTITY = 16
MAXIMUM_MUTATION_FREQUENCY = 17
POPULATION_NUCLEOTIDE_DIVERSITY = 18
SAMPLE_NUCLEOTIDE_DIVERSITY = 19


# Dictionary of Names
MEASUREMENTS_NAMES = {
    NUMBER_OF_HAPLOTYPES: "Number of Haplotypes",
    HAPLOTYPE_POPULATION: "Haplotype Population",
    NUMBER_OF_POLYMORPHIC_SITES: "Number of Polymorphic Sites",
    NUMBER_OF_MUTATIONS: "Number of Mutations",
    SHANNON_ENTROPY_NUMBER: "Shannon Entropy",
    SHANNON_ENTROPY_NUMBER_NORMALIZED_TO_N: "Shannon Entropy Normalized to N",
    SHANNON_ENTROPY_NUMBER_NORMALIZED_TO_H: "Shannon Entropy Normalized to H",
    SIMPSON_INDEX: 'Simpson Index',
    GINI_SIMPSON_INDEX: "Gini-Simpson Index",
    HILL_NUMBER_0: "Hill Number #0",
    HILL_NUMBER_1: "HIll Number #1",
    HILL_NUMBER_2: "Hill Number #2",
    HILL_NUMBER_3: "Hill Number #3",
    MINIMUM_MUTATION_FREQUENCY: "Minimum Mutation Frequency",
    MUTATION_FREQUENCY: "Mutation Frequency",
    FUNCTIONAL_ATTRIBUTE_DIVERSITY: "Functional Attribute Diversity",
    SAMPLE_NUCLEOTIDE_DIVERSITY_ENTITY: "Sample Nucleotide Diversity (Entity)",
    MAXIMUM_MUTATION_FREQUENCY: "Maximum Mutation Frequency",
    POPULATION_NUCLEOTIDE_DIVERSITY: "Population Nucleotide Diversity",
    SAMPLE_NUCLEOTIDE_DIVERSITY: "Sample Nucleotide Diversity",
}

# The number of hill numbers we want to calculate.
# if we want more hill numbers the measurement positions
# and dictionary names need to be updated.
HILL_NUMBER_LENGTH = 4
