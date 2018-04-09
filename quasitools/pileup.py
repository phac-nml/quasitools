"""
Copyright Government of Canada 2018

Written by: Matthew Fogel, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import numpy as np

# Quasitools parsers:
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.reference_parser import parse_references_from_fasta

BASES = ['A', 'C', 'T', 'G']
GAP = '-'


class Pileup_List(object):

    def __init__(self, pileups):

        """
        Creates a array of Pileup objects.

        INPUT:
            [ARRAY OF PILEUPS] [pileups] - list of Pileup objects
        RETURN:
            [None]
        POST:
            Pileup list is constructed.
        """
        self.pileups = pileups
        self.left_pos_truncated = 0
        self.right_pos_truncated = 0

    def all_have_coverage(self, position):

        """
        Determines whether all Pileup objects in the pileup list have coverage
        at the present position.

        INPUT:
            [INT] [position] - position in each pileup to check for coverage

        RETURN:
            Returns BOOL true if all pileups have coverage at the position or
            false if at least one pileup does not have coverage at the position

        POST:
            None

        """
        if any(not pileup.has_coverage(position) for pileup in self.pileups):
            return False
        else:
            return True
        # end if

    # end def

    @staticmethod
    def construct_pileup_list(file_list, reference_loc):

        """
        Creates a Pileup_List object
        INPUT:
            [FILE LOCATION TUPLE] [file_list] - files names which represent
                                                a pileup

            [FILE LOCATION] [reference_loc] - location of the reference file
        RETURN:
            [Pileup_List] - a new object containing a list of Pileup objects.
        POST:
            [None]
        """
        pileups = []
        for bam in file_list:
            pileups.append(Pileup.construct_pileup(bam, reference_loc))

        return Pileup_List(pileups)
    # end def

    def normalize_pileups(self):

        """
        This function converts the read count for each base in each four-tuple
        of bases (A, C, T, G) into a decimal proportion of the total read
        counts for that four-tuple. The bounds are between 0 and 1.
        This prevents large read counts for a single base from inflating
        the cosine simularity calculation.

        INPUT: [None]

        RETURN: [None]

        POST: All Pileups in the Pileup_List are normalized.
        """
        for pileup in self.pileups:
            pileup.normalize_pileup()
    # end def

    def get_num_left_positions_truncated(self):

        """
        This functions returns the number of left positions truncated
        since the last time truncate_output or remove_no_coverage was called.

        INPUT: [None]

        RETURN: [INT] [self.left_pos_truncated]

        POST: [None]
        """
        return self.left_pos_truncated

    def get_num_right_positions_truncated(self):

        """
        This functions returns the number of right positions truncated
        since the last time truncate_output or remove_no_coverage was called.

        INPUT: [None]

        RETURN: [INT] [self.right_pos_truncated]

        POST: [None]
        """
        return self.right_pos_truncated

    def get_pileups_as_array(self):

        """
        This function returns the pileups Pileup_List object as a
        two-dimensional array of dictionaries.

        INPUT: [None]

        RETURN: [ARRAY OF ARRAY OF DICTIONARIES] [pileup_list]

        POST: [None]
        """
        return [pileup.get_pileup_as_array_of_dictionaries()
                for pileup in self.pileups]

    def get_pileups_as_numerical_array(self):

        """
        This function returns the pileups Pileup_List object as a
        two-dimensional numerical array.

        INPUT: [None]

        RETURN: [ARRAY] [pileup_list]

        POST: [None]
        """
        return [pileup.get_pileup_as_numerical_array()
                for pileup in self.pileups]

    def get_pileup_length(self):

        """
        This function returns the (single) length of all the pileups in the
        list of pileups. The function assumes the lengths of all the pileups
        are the same.

        INPUT: [None]

        RETURN: [INT] [len(self.pileups[0])]

        POST: [None]
        """
        return len(self.pileups[0].get_pileup_as_array_of_dictionaries())

    def select_pileup_range(self, curr_start, curr_end):

        """
        Ignores all regions of the pileup before curr_start and after curr_end

        INPUT:
            [int] [curr_start] - current start position. Must be zero-indexed
            (between zero inclusive and the length of the Pileup exclusive).
            [int] [curr_end] - current end position. Must be zero-indexed
            (between zero inclusive and the length of the Pileup exclusive).
        RETURN:
            [None]
        POST:
            Positions before curr_start and after curr_end are ignored in the
            pileup list.
        """
        for pileup in self.pileups:
            pileup.select_pileup_range(curr_start, curr_end)

    def remove_no_coverage(self):

        """
        Deletes all regions of the pileup for all Pileup objects in the pileup
        list where there is no coverage for at least one pileup (all four
        bases - A, C, T, and G are absent).

        INPUT:
            [None]

        RETURN:
            [None]

        POST:
            Sections of the pileup where there is no coverage for at least one
            pileup are deleted from all pileups in the pileup list.
        """

        # First truncate end positions to determine number of contiguous
        # positions that were truncated on the left and the right.
        self.truncate_output()

        deletion_list = []

        if len(self.pileups) > 0 and self.get_pileup_length() > 0:
            # iterate through every position in reference
            for position in range(0, self.get_pileup_length()):
                # add pos'n with empty coverage in pileup to deletion_list
                if not self.all_have_coverage(position):
                    deletion_list.insert(0, position)

        for pileup in self.pileups:
            pileup.remove_pileup_positions(deletion_list)
    # end def

    def truncate_output(self):

        """
        Deletes contiguous start and end regions of the pileup for all pileups
        in the pileup list where there is no coverage (all four bases - A, C,
        T, and G are absent).

        INPUT:
            [None]

        RETURN:
            [None]

        POST:
            The pileups are truncated (sections of the pileup where there
            is no coverage are deleted from all pileups in the pileup list.
            If curr_start > curr_end the pileup_list is empty after truncation.
        """
        self.left_pos_truncated, self.right_pos_truncated = 0, 0
        deletion_list_left, deletion_list_right, deletion_list = [], [], []
        num_pos = self.get_pileup_length()

        if len(self.pileups) > 0 and self.get_pileup_length() > 0:
            # iterate through every position in reference
            for left in range(0, num_pos):
                if not self.all_have_coverage(left):
                    deletion_list_left.insert(0, left)
                    self.left_pos_truncated += 1
                else:
                    break

            for right in reversed(range(self.left_pos_truncated, num_pos)):
                if not self.all_have_coverage(right):
                    deletion_list_right.append(right)
                    self.right_pos_truncated += 1
                else:
                    break

        # example: [7 6 5 3 2 1] = [7 6 5] + [3 2 1]
        deletion_list = deletion_list_right + deletion_list_left

        for pileup in self.pileups:
            pileup.remove_pileup_positions(deletion_list)
    # end def


class Pileup(object):

    def __init__(self, pileup):

        """
        Creates a Pileup. A Pileup is a vector constructed by aligning all
        quasispecies reads to a reference and observing the number of A, G, C,
        and T nucleotides aligned at each position of the reference.

        INPUT:
            [ARRAY OF DICTIONARIES] [pileup]

        RETURN:
            [None]

        POST:
            Pileup is constructed.
        """
        self.pileup = pileup

    def has_coverage(self, position):

        """
        Determines whether the Pileup has coverage at the present position.

        INPUT:
            [INT] [position] - position in the Pileup to check for coverage

        RETURN:
            Returns BOOL true if the Pileup has coverage at the position and
            false otherwise.

        POST:
            None

        """
        curr_pos_list = [self.pileup[position].get(base, 0) for base in BASES]

        if (np.sum(curr_pos_list) == 0 or self.pileup[position] == {}):
            return False
        else:
            return True
        # end if

    # end def

    @staticmethod
    def construct_pileup(bam, reference_loc):

        """
        Creates a Pileup.
        INPUT:
            [FILE LOCATION] [bam] - file name of BAM file to create mapped
                                    read against reference

            [FILE LOCATION] [reference_loc] - location of the reference file

        RETURN:
            [ARRAY OF DICTIONARIES] [pileup] - contains read counts for each
                                               base
        POST:
            [None]
        """

        # Build the reference object.
        references = parse_references_from_fasta(reference_loc)
        new_pileup = []
        # Iterate over each reference in the reference object.
        for reference in references:
            mrc = parse_mapped_reads_from_bam(reference, bam)

            # append reads mapped against the current reference to the pileup
            # end
            new_pileup += mrc.pileup(indels=True)
        # end for
        references = None
        return Pileup(new_pileup)
    # end def

    def normalize_pileup(self):

        """
        This function converts the read count for each base in each four-tuple
        of bases (A, C, T, G) into a decimal proportion of the total read
        counts for that four-tuple. The bounds are between 0 and 1 inclusive.
        This prevents large read counts for a base from inflating
        the cosine simularity calculation.

        INPUT: [None]

        RETURN: [None]

        POST: The Pileup's values are normalized.
        """
        new_list = []
        for i in range(0, len(self.pileup)):
            curr_pos = [self.pileup[i].get(base, 0) for base in BASES]
            total = float(np.sum(curr_pos))
            items = self.pileup[i].items()

            # normalize the data for all dictionaries in the pileup
            if total > 0:
                new_list.append(
                    {key: (float(val) / total) for (key, val) in items
                     if key is not GAP})
            else:
                new_list.append(
                    {key: 0 for (key, value) in items if key is not GAP})
            # end if
        self.pileup = new_list
    # end def

    def get_pileup_as_array_of_dictionaries(self):

        """
        This function returns the pileup in the Pileup object as an array of
        dictionaries.

        INPUT: [None]

        RETURN: [ARRAY OF DICTIONARIES] [pileup]

        POST: [None]
        """
        return self.pileup

    def get_pileup_as_numerical_array(self):

        """
        This function returns the pileup in the Pileup object as a numerical
        one-dimensional array.

        INPUT: [None]

        RETURN: [ARRAY] [pileup]

        POST: [None]
        """
        # create a list of the read counts of each base at each position in the
        # pileup, zero if the base is not in the dictionary at the position
        first = 0
        last = len(self.pileup)
        numerical_array = [self.pileup[dict].get(base, 0)
                           for dict in range(first, last) for base in BASES]
        return numerical_array

    def remove_pileup_positions(self, deletion_list):

        """
        Deletes positions in the Pileup specified in deletion_list, an array
        of integers sorted in descending order.

        INPUT:
            [ARRAY] [deletion_list] - list of positions to delete in descending
                                      order of indices.

        RETURN:
            [None]

        POST:
            The specified positions in deletion_list have been removed from
            the Pileup.
        """
        for position in deletion_list:
            del self.pileup[position]
        # end for

    def select_pileup_range(self, curr_start, curr_end):

        """
        Ignores all regions of the Pileup before curr_start and after curr_end

        INPUT:
            [int] [curr_start] - current start position. Must be zero-indexed
            (between zero inclusive and the length of the Pileup exclusive).
            [int] [curr_end] - current end position. Must be zero-indexed
            (between zero inclusive and the length of the Pileup exclusive).

        RETURN:
            [None]

        POST:
            Positions before curr_start and after curr_end are ignored in the
            Pileup.
        """
        self.pileup = self.pileup[curr_start:curr_end + 1]
