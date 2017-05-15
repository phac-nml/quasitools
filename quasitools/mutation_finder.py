"""
Copyright Government of Canada 2017

Written by: Cole Peters, National Microbiology Laboratory, Public Health Agency of Canada

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
from quasitools.mutations import Mutation


class MutationFinder(object):

    def __init__(self, aa_census, min_freq, frame, filters={}):
        self.aa_census = aa_census
        self.min_freq = min_freq
        self.frame = frame
        self.filters = filters

        self.mutations = defaultdict(lambda: defaultdict(
            lambda: defaultdict(dict)))

        self._build()

    def _build(self):
        """Builds a dictionary of amino acid mutations from an AACensus object
        """

        ref_seq = self.aa_census.mapped_read_collections[0].reference.seq

        ref_seq = ref_seq[:(len(ref_seq) - (len(ref_seq) % 3))]

        ref_aa = Seq(ref_seq).translate()

        if self.min_freq is not None:
            filter_id = "mf%s" % self.min_freq
            self.filters[filter_id] = {'expression': 'mf<%s' % (self.min_freq),
                                       'result': True}
        else:
            raise ValueError("min_freq must be defined.")

        for i in range(0, len(ref_aa)):
            coverage = self.aa_census.coverage_at(self.frame, i)
            for confidence in (CONFIDENT, UNCONFIDENT):
                for aa in self.aa_census.aminos_at(self.frame, i, confidence):
                    if aa != ref_aa[i]:
                        gene = None

                        for name in self.aa_census.genes:
                            if (i >= self.aa_census.genes[name]["start"] // 3
                                and i <= (self.aa_census.genes[name]["end"]
                                          - 2) // 3):

                                gene = self.aa_census.genes[name]
                                gene_name = name

                        if gene is not None:
                            filter = "."
                            frequency = self.aa_census.amino_frequency_at(
                                self.frame, i, aa, confidence
                            ) / float(coverage)

                            if frequency < self.min_freq:
                                filter = filter_id
                            else:
                                filter = "PASS"

                            mutation = Mutation(gene_name, ref_aa[i],
                                                (i - (gene["start"] // 3) + 1),
                                                frequency, filter)

                            self.mutations[i][confidence][aa] = mutation

    def to_hmcf_file(self, confidence, mutation_db=None):
        """Builds a report of all amino acid mutations found from the AACensus
        object with the given confidence. Returns a string containing the
        report.
        """

        ref_seq = self.aa_census.mapped_read_collections[0].reference.seq

        ref_codon_array = re.findall(".{3}", ref_seq)

        d = date.today()

        report = "##fileformat=HMCFv1\n"
        report += "##fileDate=%s\n" % (d.strftime("%Y%m%d"))
        report += "##source=quasitools\n"
        report += "##reference=%s\n" % (self.aa_census.ref_file)

        info_line = "##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n"

        report += info_line % ("WC", ".", "String", "WildType Codon")
        report += info_line % ("MC", ".", "String", "Mutation Codon")
        report += info_line % ("MCF", ".", "String", "Mutant Codon Frequency, "
                               "for each Mutant Codon, in the same order as "
                               "listed.")

        filter_line = "##FILTER=<ID=%s,Description=\"Set if %s; %s\">\n"

        for id, filter in self.filters.items():
            report += filter_line % (id, filter['result'],
                                     filter['expression'])

        report += "#GENE\tCATEGORY\tSURVEILLANCE\tTYPE\tWILDTYPE\tPOS\t" \
                  "MUTANT\tFILTER\tMUTANT_FREQ\tCOVERAGE\tINFO\n"

        for ref_codon_pos in self.mutations:
            if mutation_db is not None:
                dr_mutations = mutation_db.mutations_at(ref_codon_pos)
            else:
                dr_mutations = None

            for aa in self.mutations[ref_codon_pos][confidence]:
                dr_mutation = None

                if dr_mutations is not None and aa in dr_mutations:
                    dr_mutation = dr_mutations[aa]

                # Ignore incomplete codons (ones with a gap)
                if aa.lower() != "x":
                    coverage = \
                        self.aa_census.coverage_at(self.frame, ref_codon_pos)

                    mutation = self.mutations[ref_codon_pos][confidence][aa]

                    category = "."
                    surveillance = "."

                    if dr_mutation is not None:
                        category = dr_mutation.category
                        surveillance = dr_mutation.surveillance

                    report += ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%0.4f\t%s\t"
                               % (mutation.gene, category,
                                  surveillance, "mutation", mutation.wildtype,
                                  mutation.gene_pos, aa, mutation.filter,
                                  mutation.frequency, coverage))

                    report += "WC=%s;" % ref_codon_array[ref_codon_pos].lower()

                    mc = "MC="
                    mcf = "MCF="
                    for codon in self.aa_census.amino_to_codons_at(
                            self.frame, ref_codon_pos, aa, confidence):

                        mc += "%s," % codon

                        frequency = \
                            self.aa_census.codon_frequency_for_amino_at(
                                self.frame, ref_codon_pos, aa, confidence,
                                codon)

                        mcf += "%0.4f," % (float(frequency) / coverage)

                    # Remove trailing comma
                    mc = mc[:-1]
                    mcf = mcf[:-1]

                    report += "%s;%s\n" % (mc, mcf)

        return report

    def report_dr_mutations(self, mutation_db, reporting_threshold):
        """Builds a report of all drug resistant amino acid mutations present
        in the AACensus object using a MutationDB object. A string containing
        the report is then returned.
        """

        report = ("Gene,Category,Surveillance,Wildtype,Position,Mutation,"
                  "Mutation Frequency,Coverage\n")

        # Loop through the mutation database and report on present mutations
        for dr_mutation_pos in mutation_db.positions():
            dr_mutations = mutation_db.mutations_at(dr_mutation_pos)

            for dr_mutation in dr_mutations:
                coverage = self.aa_census.coverage_at(self.frame,
                                                      dr_mutation_pos)

                if dr_mutation_pos in self.mutations:
                    if CONFIDENT in self.mutations[dr_mutation_pos]:
                        if (dr_mutation in
                                self.mutations[dr_mutation_pos][CONFIDENT] and
                                self.mutations[dr_mutation_pos][CONFIDENT]
                                [dr_mutation].filter == "PASS"):

                            mutation_freq = (
                                self.mutations[dr_mutation_pos][CONFIDENT]
                                [dr_mutation].frequency
                            ) * 100

                            if mutation_freq > reporting_threshold:
                                report += (
                                    "%s,%s,%s,%s,%s,%s,%0.2f,%s\n"
                                    % (dr_mutations[dr_mutation].gene,
                                       dr_mutations[dr_mutation].category,
                                       dr_mutations[dr_mutation].surveillance,
                                       dr_mutations[dr_mutation].wildtype,
                                       dr_mutations[dr_mutation].gene_pos,
                                       dr_mutation,
                                       mutation_freq,
                                       coverage))

        return report
