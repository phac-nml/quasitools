UNDEFINED = ""

# Measurement Positions
NUMBER_OF_HAPLOTYPES = 0
NUMBER_OF_POLYMORPHIC_SITES = 1
NUMBER_OF_MUTATIONS = 2
SHANNON_ENTROPY_NUMBER = 3
SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_N = 4
SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_H = 5
SIMPSON_INDEX = 6
GINI_SIMPSON_INDEX = 7
HILL_NUMBER_0 = 8
HILL_NUMBER_1 = 9
HILL_NUMBER_2 = 10
HILL_NUMBER_3 = 11
MINIMUM_MUTATION_FREQUENCY = 12
MUTATION_FREQUENCY = 13
FUNCTIONAL_ATTRIBUTE_DIVERSITY = 14
SAMPLE_NUCLEOTIDE_DIVERSITY_Entity = 15
MAXIMUM_MUTATION_FREQUENCY = 16
POPULATION_NUCLEOTIDE_DIVERSITY = 17
SAMPLE_NUCLEOTIDE_DIVERSITY = 18


# Dictionary of Names
MEASUREMENTS_NAMES = {
    NUMBER_OF_HAPLOTYPES: "Number of Haplotypes",
    NUMBER_OF_POLYMORPHIC_SITES: "Number of Polymorphic Sites",
    NUMBER_OF_MUTATIONS: "Number of Mutations",
    SHANNON_ENTROPY_NUMBER: "Shannon Entropy",
    SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_N: "Shannon Entropy Localized to N",
    SHANNON_ENTROPY_NUMBER_LOCALIZED_TO_H: "Shannon Entropy Localized to H",
    SIMPSON_INDEX: 'Simpson Index',
    GINI_SIMPSON_INDEX: "Gini Simpson Index",
    HILL_NUMBER_0: "Hill Number #0",
    HILL_NUMBER_1: "HIll Number #1",
    HILL_NUMBER_2: "Hill Number #2",
    HILL_NUMBER_3: "Hill Number #3",
    MINIMUM_MUTATION_FREQUENCY: "Minimum Mutation Frequency",
    MUTATION_FREQUENCY: "Mutation Frequency",
    FUNCTIONAL_ATTRIBUTE_DIVERSITY: "Functional Attribute Diversity",
    SAMPLE_NUCLEOTIDE_DIVERSITY_Entity: "Sample Nucleotide Diversity Entity",
    MAXIMUM_MUTATION_FREQUENCY: "Maximum Mutation Frequency",
    POPULATION_NUCLEOTIDE_DIVERSITY: "Population Nucleotide Diversity",
    SAMPLE_NUCLEOTIDE_DIVERSITY: "Sample Nucleotide Diversity",
}

# The number of hill numbers we want.
# if we want more hill numbers the measurement positions
# and dictionary names need to be updated.
HILL_NUMBER_LENGTH = 4
