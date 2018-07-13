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

BASES = ['A', 'C', 'T', 'G']
GAP = '-'


class Pileup_List(object):
    """This class contains a list of Pileups."""

    def __init__(self, pileups):
        """
        Create a array of Pileup objects.

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
        Determine if all Pileup objects have coverage at the present position.

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

    def normalize_pileups(self):
        """
        Normalize pileups.

        Convert the read count for each base in each four-tuple
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
        Return number of left positions truncated.

        Returned value was saved the most recent time truncate_output or
        remove_no_coverage was called.

        INPUT: [None]

        RETURN: [INT] [self.left_pos_truncated]

        POST: [None]

        """
        return self.left_pos_truncated

    def get_num_right_positions_truncated(self):
        """
        Return number of right positions truncated.

        Returned value was saved the most recent time truncate_output or
        remove_no_coverage was called.

        INPUT: [None]

        RETURN: [INT] [self.right_pos_truncated]

        POST: [None]

        """
        return self.right_pos_truncated

    def get_pileups_as_array(self):
        """
        Return Pileup_List object as a two-dimensional array of dictionaries.

        INPUT: [None]

        RETURN: [ARRAY OF ARRAY OF DICTIONARIES] [pileup_list]

        POST: [None]

        """
        return [pileup.get_pileup_as_array_of_dictionaries()
                for pileup in self.pileups]

    def get_pileups_as_numerical_array(self):
        """
        Return Pileup_List object as a two-dimensional numerical array.

        INPUT: [None]

        RETURN: [ARRAY] [pileup_list]

        POST: [None]

        """
        return [pileup.get_pileup_as_numerical_array()
                for pileup in self.pileups]

    def get_pileup_length(self):
        """
        Get pileup length.

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
        Ignore all regions of the pileup before curr_start and after curr_end.

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
        Remove no coverage regions.

        Delete all regions of the pileup for all Pileup objects in the pileup
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
        Truncate output.

        Delete contiguous start and end regions of the pileup for all pileups
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
    """This class stores a pileup and utilities for modifying the pileup."""

    def __init__(self, pileup):
        """
        Create a Pileup.

        The object represents the Pileup of reads
        mapped against a reference file.

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
        Determine whether the Pileup has coverage at the present position.

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

    def normalize_pileup(self):
        """
        Normalize pileup.

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
        Return the pileup in the Pileup object as an array of dictionaries.

        INPUT: [None]

        RETURN: [ARRAY OF DICTIONARIES] [pileup]

        POST: [None]

        """
        return self.pileup

    def get_pileup_as_numerical_array(self):
        """
        Return the pileup in the Pileup object as a numerical 1D array.

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
        Remove pileup positions.

        Delete positions in the Pileup specified in deletion_list, an array
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
        Ignore all regions of the Pileup before curr_start and after curr_end.

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

    """
    # =============================================================================
    # =============================================================================
    """
    def build_consensus(self):

        consensus = []

        for position in self.pileup:

            sorted_position = sorted(position, key=position.get, reverse=True)
            base = sorted_position[0]

            consensus.append(base)

        return consensus

    """
    # =============================================================================
    # =============================================================================
    """
    def count_unique_mutations(self):

        # !!This assumes there are no gaps in the passed pileup!!

        # We need the number mutations at all mutation sites.
        # These are positions in the pileup with at least 1 disagreement.
        unique_mutations = 0

        for position in self.pileup:

            unique_mutations += len(position) - 1 # Number of different bases at position.

        return unique_mutations

    """
    # =============================================================================
    # =============================================================================
    """
    def count_polymorphic_sites(self):

        # !!This assumes there are no gaps in the passed pileup!!

        # We need the number of polymorphic sites.
        # These are positions in the pileup with at least 1 disagreement.
        polymorphic_sites = 0

        for position in self.pileup:

            if len(position) > 1:

                polymorphic_sites += 1

        return polymorphic_sites



