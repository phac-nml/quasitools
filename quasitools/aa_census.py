"""
Copyright Government of Canada 2017

Written by: Cole Peters, National Microbiology Laboratory,
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

from Bio.Seq import Seq
from collections import defaultdict

CONFIDENT = 0
UNCONFIDENT = 1


class AACensus(object):
    def __init__(self, ref_file, mapped_read_collections, genes, frames):
        self.ref_file = ref_file
        self.mapped_read_collections = mapped_read_collections
        self.genes = genes
        self.positional_freq = defaultdict(lambda: defaultdict(
            lambda: defaultdict(int)))
        self.frames = frames

        self._build()

    def _build(self):
        """Builds the amino acid census"""

        for frame in self.frames:
            for mrc in self.mapped_read_collections:
                # For each mapped_read translate in the frame, store the
                # frequencies of the coded amino acid.
                for name, mapped_read in mrc.mapped_reads.items():
                    # Get the start of the first codon and the end of the last
                    # codon.
                    start = mapped_read.codon_start(frame)
                    end = mapped_read.codon_end(frame)

                    confidence = {}

                    read_wo_ins = mrc.reference.sub_seq(
                        start, start + (end - start)).lower()

                    # Build up a base from our reference to apply our
                    # differences to (ignore insertions)
                    for pos in mapped_read.differences:
                        difference = (mapped_read.differences[pos])[:1]

                        # If we have a difference at this position, (the
                        # difference isn't equal to ".") apply it.
                        if difference != '.':
                            # If the difference is lower case, that means it
                            # was filtered out and as such it is unconfident.
                            if difference != difference.upper():
                                confidence[((pos - start) // 3)] = UNCONFIDENT

                            read_wo_ins = read_wo_ins[:(pos - start)] + \
                                difference.upper() + \
                                read_wo_ins[(pos - start + 1):]

                    # Translate the read
                    read_wo_ins = \
                        read_wo_ins[:len(read_wo_ins) - (len(read_wo_ins) % 3)]

                    read_aas = Seq(read_wo_ins.replace("-", "N")).translate()

                    # Calculate the first position
                    start_aa = ((start) // 3)

                    # Add our codons and aas to the aa census
                    for i in range(0, len(read_aas)):
                        # Get the confidence for this codon/aa
                        aa_confidence = CONFIDENT
                        if i in confidence:
                            aa_confidence = confidence[i]

                        # Retrieve the codon which produces the aa (have to
                        # convert aa positions to codon positions)
                        codon = read_wo_ins[i*3:(i*3)+3]
                        aa = read_aas[i]

                        # If we have a complete codon (codon doesn't contain a
                        # "-"), add it to our frequencies.
                        if codon.find("-") == -1:
                            self.positional_freq[
                                (frame, (start_aa + i), aa_confidence)
                            ][aa][codon] += 1
                        else:
                            self.positional_freq[
                                (frame, (start_aa + i), UNCONFIDENT)
                            ]['X'][codon] += 1

    def aminos_at(self, frame, position, confidence):
        """Returns a list of all the amino acids at the given position"""
        return self.positional_freq[(frame, position, confidence)].keys()

    def amino_to_codons_at(self, frame, position, aa, confidence):
        """Finds the codons that produce the amino acid at the given position
        """
        return self.positional_freq[(frame, position, confidence)][aa].keys()

    def codon_frequency_for_amino_at(self, frame, position, aa, confidence,
                                     codon):
        """Finds the frequency for a codon of an amino acid at the given
        position.
        """

        frequency = 0

        if codon in self.positional_freq[(frame, position, confidence)][aa]:
            frequency = self.positional_freq[
                (frame, position, confidence)][aa][codon]

        return frequency

    def amino_frequency_at(self, frame, position, aa, confidence):
        """Finds the frequency for an amino acid at the given position"""

        frequency = 0

        if aa in self.positional_freq[(frame, position, confidence)]:
            for codon in self.positional_freq[
                    (frame, position, confidence)][aa]:

                frequency += self.positional_freq[
                    (frame, position, confidence)][aa][codon]

        return frequency

    def coverage_at(self, frame, position):
        """Finds the coverage at a given position"""

        coverage = 0

        for confidence in (CONFIDENT, UNCONFIDENT):
            for aa in self.positional_freq[(frame, position, confidence)]:
                coverage += \
                    self.amino_frequency_at(frame, position, aa, confidence)

        return coverage

    def coverage(self, frames):
        """Calculates the coverage and returns it as a string"""
        # "//" specifies integer division for python 3
        length = len(self.mapped_read_collections[0].reference.seq) // 3

        coverage_csv = ""
        for frame in frames:
            coverage_csv += "frame: %i\n" % frame

            for i in range(0, length):
                local_coverage = self.coverage_at(frame, i)
                coverage_csv += "%s,%s\n" % (i + 1, local_coverage)

            coverage_csv += "\n"

        return coverage_csv
