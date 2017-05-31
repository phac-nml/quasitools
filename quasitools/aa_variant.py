"""
Copyright Government of Canada 2017

Written by: Eric Chubaty, National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import re
from datetime import date
from collections import defaultdict
from Bio.Seq import Seq
from quasitools.aa_census import CONFIDENT, UNCONFIDENT
from quasitools.variant import Variant, VariantCollection
from numpy import array as np_array


class AAVariant(Variant):

    def __init__(self, gene=".", freq=0, coverage=0, census_ind=0, info=dict(), **kwargs):
        """Add additional fields to Variant for AAVariant"""
        super(AAVariant, self).__init__(**kwargs)

        self.info = info
        self.gene = gene
        self.freq = freq
        self.coverage = coverage
        self.census_ind = census_ind

    def __info_to_str(self):
        """Convert info dict to info string for hmcf entry."""
        return "WC=%s;MC=%s;MCF=%s;CAT=%s;SRVL=%s" % (
            self.info['WC'],
            self.info['MC'],
            self.info['MCF'],
            self.info['CAT'],
            self.info['SRVL']
        )

    def to_hmcf_entry(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%0.4f\t%s\t%s\n" % (
            self.chrom, self.gene, self.id,
            self.ref, self.pos, self.alt,
            self.filter, self.freq, self.coverage, self.__info_to_str()
        )


class AAVariantCollection(VariantCollection):

    def __init__(self, references):
        """Add additional field to VariantCollection for AAVariantCollection"""
        super(AAVariantCollection, self).__init__(references)
        self.variants = defaultdict(
            lambda: defaultdict(lambda: defaultdict(dict)))

    @classmethod
    def from_aacensus(cls, aa_census, *args):
        """Build the AAVariantCollection from any number of
        AACensus objects"""

        # Handle aa_census being an array or single element nicely
        aa_census = np_array([aa_census]).flatten()

        var_collect = cls(aa_census)

        # Build up the Collection of AAVariants from many census
        for census_ind, census in enumerate(aa_census):
            for frame in census.frames:
                census_genes = census.genes
                ref_seq = census.mapped_read_collections[0].reference.seq
                ref_seq = ref_seq[:(len(ref_seq) - (len(ref_seq) % 3))]
                ref_aa = Seq(ref_seq).translate()

                # Used for WC in info field for each variant
                ref_codon_array = re.findall(".{3}", ref_seq)

                # Build up what will be the key to the dictionary
                for ref_codon_pos in range(0, len(ref_aa)):
                    coverage = census.coverage_at(frame, ref_codon_pos)

                    for confidence in (CONFIDENT, UNCONFIDENT):
                        for aa in census.aminos_at(
                            frame, ref_codon_pos, confidence
                        ):

                            if aa != ref_aa[ref_codon_pos]:
                                gene = None

                                # Start retrieving values for this AAVariant
                                for name in census.genes:
                                    if (ref_codon_pos >=
                                            census.genes[name]["start"] // 3
                                            and ref_codon_pos <=
                                            (census.genes[name]["end"] - 2) // 3):

                                        gene = census.genes[name]
                                        gene_name = name

                                if (gene is not None
                                        and gene['frame'] == frame):

                                    frequency = census.amino_frequency_at(
                                        frame, ref_codon_pos, aa, confidence
                                    ) / float(coverage)

                                    chrom = census.genes[gene_name]['chrom']

                                    # Find MC and MCF
                                    mc = ""
                                    mcf = ""

                                    for codon in census.amino_to_codons_at(
                                            frame, ref_codon_pos, aa, confidence
                                    ):

                                        mc += "%s," % codon

                                        freq_mcf = census.codon_frequency_for_amino_at(
                                            frame, ref_codon_pos,
                                            aa, confidence, codon
                                        )

                                        mcf += "%0.4f," % (float(freq_mcf) /
                                                           coverage)

                                    # Create AAVariant & slap it in the
                                    # collection
                                    mutation = AAVariant(chrom=chrom,
                                                         gene=gene_name,
                                                         id="mutation",
                                                         ref=ref_aa[
                                                             ref_codon_pos],
                                                         alt=aa,
                                                         freq=frequency,
                                                         coverage=coverage,
                                                         census_ind=census_ind,
                                                         pos=(ref_codon_pos - (
                                                            gene['start'] // 3
                                                         ) + 1),
                                                         info={
                                                             'WC': ref_codon_array[ref_codon_pos].lower(),
                                                             'MC': mc[:-1],
                                                             'MCF': mcf[:-1],
                                                             'CAT': ".",
                                                             'SRVL': "."
                                                         })

                                    var_collect.variants[chrom][ref_codon_pos][
                                        confidence][aa] = mutation

        return var_collect

    def to_hmcf_file(self, confidence, mutation_db=None):
        """Build a string representation of our AAVariant objects
        (i.e. a hmcf file)."""

        # Init
        d = date.today()

        # Header
        report = "##fileformat=HMCFv2\n"
        report += "##fileDate=%s\n" % (d.strftime("%Y%m%d"))
        report += "##source=quasitools\n"

        # Could have many reference files (aa_census)
        for refs in self.references:
            report += "##reference=%s\n" % (refs.ref_file)

        info_line = "##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n"

        report += info_line % ("WC", ".", "String", "WildType Codon")
        report += info_line % ("MC", ".", "String", "Mutation Codon")
        report += info_line % ("MCF", ".", "String", "Mutant Codon Frequency, "
                               "for each Mutant Codon, in the same order as "
                               "listed.")
        report += info_line % ("CAT", ".", "String",
                               "Drug Resistance Category")
        report += info_line % ("SRVL", ".", "String",
                               "Drug Resistance Surveillance")

        filter_line = "##FILTER=<ID=%s,Description=\"Set if %s; %s\">\n"

        for id, filter in self.filters.items():
            report += filter_line % (id, filter['result'],
                                     filter['expression'])

        report += "#CHROM\tGENE\tTYPE\tWILDTYPE\tPOS\t" \
            "MUTANT\tFILTER\tMUTANT_FREQ\tCOVERAGE\tINFO\n"

        # Body
        for chrom in self.variants:
            for ref_codon_pos in self.variants[chrom]:

                # Work with or without mutation db
                if mutation_db is not None:
                    dr_mutations = mutation_db.mutations_at(ref_codon_pos)
                else:
                    dr_mutations = None

                for aa in self.variants[chrom][ref_codon_pos][confidence]:
                    # Init drug resistance variables
                    dr_mutation = None
                    category = "."
                    surveillance = "."

                    if dr_mutations is not None and aa in dr_mutations:
                        dr_mutation = dr_mutations[aa]
                        category = dr_mutation.category
                        surveillance = dr_mutation.surveillance

                    # Ignore incomplete codons (ones with a gap)
                    if aa.lower() != "x":
                        mutation = self.variants[chrom][
                            ref_codon_pos][confidence][aa]

                        mutation.info['CAT'] = category
                        mutation.info['SRVL'] = surveillance

                        # Add this mutation to the report!
                        report += mutation.to_hmcf_entry()

        # Return string of report without nl char at the end
        return report[:-1]

    def report_dr_mutations(self, mutation_db, reporting_threshold):
        """Builds a report of all drug resistant amino acid mutations present
        in the AACensus object using a MutationDB object. A string containing
        the report is then returned.
        """

        report = ("Chromosome,Gene,Category,"
                  "Surveillance,Wildtype,Position,Mutation,"
                  "Mutation Frequency,Coverage\n")

        # Loop through the mutation database and report on present mutations
        for census in self.references:

            for dr_mutation_pos in mutation_db.positions():
                dr_mutations = mutation_db.mutations_at(dr_mutation_pos)
                #coverage = census.coverage_at(self.frame, dr_mutation_pos)

                for name in census.genes:
                    if (dr_mutation_pos >= census.genes[name]['start'] // 3
                            and dr_mutation_pos <=
                            (census.genes[name]['end'] - 2) // 3):

                        gene_name = name

                chrom = census.genes[gene_name]['chrom']

                for dr_mutation in dr_mutations:

                    if dr_mutation_pos in self.variants[chrom]:
                        if CONFIDENT in self.variants[chrom][dr_mutation_pos]:
                            if (dr_mutation in
                                    self.variants[chrom][
                                        dr_mutation_pos][CONFIDENT] and
                                    self.variants[chrom][
                                        dr_mutation_pos][CONFIDENT]
                                    [dr_mutation].filter == "PASS"):

                                mutation_freq = (
                                    self.variants[chrom][
                                        dr_mutation_pos][
                                        CONFIDENT][
                                        dr_mutation].freq
                                ) * 100

                                coverage = self.variants[chrom][
                                    dr_mutation_pos][
                                    CONFIDENT][
                                    dr_mutation].coverage

                                if mutation_freq > reporting_threshold:
                                    report += (
                                        "%s,%s,%s,%s,%s,%s,%s,%0.2f,%s\n"
                                        % (chrom,
                                           dr_mutations[dr_mutation].gene,
                                           dr_mutations[dr_mutation].category,
                                           dr_mutations[
                                               dr_mutation].surveillance,
                                           dr_mutations[dr_mutation].wildtype,
                                           dr_mutations[dr_mutation].gene_pos,
                                           dr_mutation,
                                           mutation_freq,
                                           coverage))

        return report[:-1]

    def filter(self, id, expression, result):
        """Apply filter to variants given an id, expression and result."""
        self.filters[id] = {'expression': expression, 'result': result}

        # only allow simple expressions for the time being i.e. DP>30
        (attribute, operator, value) = re.split('([><=!]+)', expression)

        for chrom in self.variants:
            for ref_codon_pos in self.variants[chrom]:

                for confidence in self.variants[chrom][ref_codon_pos]:
                    for aa in self.variants[chrom][ref_codon_pos][confidence]:
                        attribute_value = None

                        variant = self.variants[chrom][
                            ref_codon_pos][confidence][aa]

                        if hasattr(variant, attribute.lower()):
                            attribute_value = eval(
                                "variant.%s" % attribute.lower())
                        else:
                            attribute_value = variant.info[attribute.upper()]

                        if eval("%s %s %s" % (
                            attribute_value, operator, value
                        )) != result:
                            if variant.filter == '.':
                                variant.filter = 'PASS'
                        else:
                            if variant.filter == '.' or \
                                    variant.filter == 'PASS':
                                variant.filter = id
                            else:
                                variant.filter += ";%s" % id
