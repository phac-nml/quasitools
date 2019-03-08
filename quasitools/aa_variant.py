"""
Copyright Government of Canada 2017

Written by: Eric Chubaty, National Microbiology Laboratory,
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

import re
import datetime
from PyAAVF.model import AAVF, Record, Info, Filter
from collections import defaultdict
from Bio.Seq import Seq
from quasitools.aa_census import CONFIDENT, UNCONFIDENT
from quasitools.variant import Variant, VariantCollection
from numpy import array as np_array


class AAVariant(Variant):

    def __init__(self, gene=".", freq=0,
                 coverage=0, census_ind=0,
                 info=dict(), **kwargs):
        """Add additional fields to Variant for AAVariant"""
        super(AAVariant, self).__init__(**kwargs)

        self.info = info
        self.gene = gene
        self.freq = freq
        self.coverage = coverage
        self.census_ind = census_ind

    def info_to_dict(self):
        info_dict = {}

        info_dict['RC'] = self.info['RC']
        info_dict['AC'] = self.info['AC']
        info_dict['ACF'] = self.info['ACF']
        info_dict['CAT'] = self.info['CAT']
        info_dict['SRVL'] = self.info['SRVL']

        return info_dict


class AAVariantCollection(VariantCollection):

    def __init__(self, references):
        """Add additional field to VariantCollection for AAVariantCollection"""
        super(AAVariantCollection, self).__init__(references)
        self.variants = defaultdict(
            lambda: defaultdict(lambda: defaultdict(dict)))

    @classmethod
    def from_aacensus(cls, aa_census):
        """Build the AAVariantCollection from any number of
        AACensus objects"""

        # Handle aa_census being an array or single element nicely
        aa_census = np_array([aa_census]).flatten()

        var_collect = cls(aa_census)

        # Build up the Collection of AAVariants from many census
        for census_ind, census in enumerate(aa_census):
            # For each gene in this census
            for gene_key in census.genes:
                gene = census.genes[gene_key]
                frame = gene['frame']

                # Find reference sequence for this frame
                ref_seq = census.mapped_read_collections[0].reference.seq
                ref_seq = ref_seq[
                    frame:((len(ref_seq) - frame) -
                           (len(ref_seq) - frame) % 3 + frame)
                ]

                # Turn sequence to amino acids
                ref_aa = Seq(ref_seq).translate()

                # Used for WC in info field for each variant
                ref_codon_array = re.findall(".{3}", ref_seq)

                # Build up what will be the key to the dictionary
                for ref_codon_pos in range(gene['start'] // 3,
                                           gene['end'] // 3 - 2):
                    coverage = census.coverage_at(frame, ref_codon_pos)

                    for confidence in (CONFIDENT, UNCONFIDENT):
                        for aa in census.aminos_at(
                                frame, ref_codon_pos, confidence
                        ):

                            # Only add it if it is a mutation (it differs)
                            if aa != ref_aa[ref_codon_pos]:

                                # Grab values for this mutation

                                frequency = census.amino_frequency_at(
                                    frame, ref_codon_pos, aa, confidence
                                ) / float(coverage)

                                chrom = gene['chrom']

                                # Find AC (Allele Count) and
                                # ACF (Allele Count Frequency)
                                ac = ""
                                acf = ""

                                for codon in census.amino_to_codons_at(
                                        frame, ref_codon_pos, aa, confidence
                                ):

                                    ac += "%s," % codon

                                    freq_acf = census\
                                        .codon_frequency_for_amino_at(
                                            frame, ref_codon_pos,
                                            aa, confidence, codon
                                        )

                                    acf += "%0.4f," % (float(freq_acf) /
                                                       coverage)

                                # Create AAVariant & slap it in the
                                # collection
                                mutation = AAVariant(chrom=chrom,
                                                     gene=gene_key,
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
                                                         'RC': ref_codon_array[
                                                             ref_codon_pos
                                                         ].lower(),
                                                         'AC': ac[:-1],
                                                         'ACF': acf[:-1],
                                                         'CAT': ".",
                                                         'SRVL': "."
                                                     })

                                var_collect.variants[chrom][ref_codon_pos][
                                    confidence][aa] = mutation

        return var_collect

    def to_aavf_obj(self, source_command, reference_filename, confidence):
        metadata = {}
        infos = {}
        filters = {}
        column_headers = ["CHROM", "GENE", "POS", "REF", "ALT", "FILTER",
                          "ALT_FREQ", "COVERAGE", "INFO"]
        record_list = []

        # Build metadata
        metadata['fileformat'] = 'AAVFv1.0'
        metadata['fileDate'] = datetime.date.today().strftime("%Y%m%d")
        metadata['source'] = "quasitools:%s" % source_command
        metadata['reference'] = [reference_filename]

        # Build infos
        infos["RC"] = Info("RC", "1", "String", "Reference Codon",
                           source_command, "1.0")
        infos["AC"] = Info("AC", ".", "String", "Alternate Codon",
                           source_command, "1.0")
        infos["ACF"] = Info("ACF", ".", "Float", ("Alternate Codon Frequency,"
                                                  "for each Alternate Codon,"
                                                  "in the same order as"
                                                  "listed."),
                            source_command, "1.0")
        infos["CAT"] = Info("CAT", ".", "String", "Drug Resistance Category",
                            source_command, "1.0")
        infos["SRVL"] = Info("SRVL", ".", "String",
                             "Drug Resistance Surveillance",
                             source_command, "1.0")

        # Build filters
        filters["af0.01"] = Filter("af0.01", "Set if True; alt_freq<0.01")

        # Build record_list from self.variants
        for chrom in self.variants:
            for ref_codon_pos in self.variants[chrom]:
                for aa in self.variants[chrom][ref_codon_pos][confidence]:
                    variant = (self.variants[chrom][ref_codon_pos]
                               [confidence][aa])

                    if aa.lower() != "x":
                        curr_record = Record(variant.chrom, variant.gene,
                                             variant.pos, variant.ref,
                                             variant.alt, variant.filter,
                                             "%0.4f" % variant.freq,
                                             variant.coverage,
                                             variant.info_to_dict())

                        record_list.append(curr_record)

        return AAVF(metadata, infos, filters, column_headers, record_list)

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

    def apply_mutation_db(self, mutation_db):
        """Apply the mutation database to the variant collection
        to update each variants category and surveillance variable within it.

        Assumes mutation_db != None
        """

        # Iterate over the keys in variants
        for chrom in self.variants:
            for ref_codon_pos in self.variants[chrom]:

                dr_mutations = mutation_db.mutations_at(ref_codon_pos)

                for confidence in (CONFIDENT, UNCONFIDENT):
                    for aa in self.variants[chrom][
                            ref_codon_pos][confidence]:

                        # Init drug resistance variables
                        dr_mutation = None
                        category = "."
                        surveillance = "."

                        # Assign cat and srvl if it's in the mutation db
                        if aa in dr_mutations:
                            dr_mutation = dr_mutations[aa]
                            category = dr_mutation.category
                            surveillance = dr_mutation.surveillance

                        # Get the mutation to update
                        mutation = self.variants[chrom][
                            ref_codon_pos][confidence][aa]

                        # Update the mutation
                        mutation.info['CAT'] = category
                        mutation.info['SRVL'] = surveillance

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
