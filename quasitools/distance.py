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
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.spatial.distance import cosine

BASES = ['A', 'C', 'T', 'G']
GAP = '-'


class Pileup_List(object):

    def __init__(self, pileups):

        """
        Creates a array of Pileup objects (which are arrays of dictionaries)

        INPUT:
            [ARRAY OF PILEUPS] [pileups] - list of Pileup objects
        RETURN:
            [None]
        POST:
            Pileup list is constructed.
        """
        self.pileups = pileups

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
    def construct_array_of_pileups(file_list, reference_loc):

        """
        Creates an array of Pileup objects (which contain arrays of
        dictionaries)
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

    def normalize_pileup(self):

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

    def get_pileups_as_array(self):

        """
        This function returns the pileups Pileup_List object as a 2D array.

        INPUT: [None]

        RETURN: [ARRAY OF ARRAY OF DICTIONARIES] [pileup_list]

        POST: [None]
        """
        return [pileup.get_pileup_as_array() for pileup in self.pileups]

    def get_pileup_length(self):

        """
        This function returns the length of the first pileup in pileups.
        The length of each pileup should be the same.

        INPUT: [None]

        RETURN: [INT] [len(self.pileups[0])]

        POST: [None]
        """
        return self.pileups[0].get_pileup_length()

    def select_pileup_range(self, curr_start, curr_end):

        """
        Ignores all regions of the pileup before curr_start and after curr_end

        INPUT:
            [int] [curr_start] - current start position. Must be between zero
                                 inclusive and the length of the pileup,
                                 exclusive.
            [int] [curr_end] - current end position. Must be between zero
                               inclusive and the length of the pileup,
                               exclusive.
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
        list where there is no coverage (all four bases - A, C, T, and G are
        absent).

        INPUT:
            [None]

        RETURN:
            [None]

        POST:
            The pileups are truncated (sections of the pileup where there
            is no coverage are deleted from all pileups in the pileup list).
        """
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
        deletion_list_left, deletion_list_right, deletion_list = [], [], []
        num_pos = self.get_pileup_length()
        if len(self.pileups) > 0 and self.get_pileup_length() > 0:
            # iterate through every position in reference
            left_key = -1
            for left in range(0, num_pos):
                if not self.all_have_coverage(left):
                    deletion_list_left.insert(0, left)
                    left_key = left
                else:
                    break
            for right in reversed(range(left_key + 1, num_pos)):
                if not self.all_have_coverage(right):
                    deletion_list_right.append(right)
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
        Creates a Pileup (which is an array of dictionaries)

        INPUT:
            [ARRAY OF DICTIONARIES] [pileup]
        RETURN:
            [None]
        POST:
            Pileup is constructed.
        """
        self.pileup = pileup
        self.pileup_length = len(pileup) - 1

    def has_coverage(self, position):

        """
        Determines whether the Pileup has coverage at the present position.

        INPUT:
            [INT] [position] - position in the Pileup to check for coverage

        RETURN:
            Returns BOOL true if the Pileup has coverage at the position

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
        counts for that four-tuple. The bounds are between 0 and 1.
        This prevents large read counts for a base from inflating
        the cosine simularity calculation.

        INPUT: [None]

        RETURN: [None]

        POST: The Pileup object's data is normalized.
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

    def get_pileup_as_array(self):

        """
        This function returns the pileup in the Pileup object.

        INPUT: [None]

        RETURN: [ARRAY] [pileup]

        POST: [None]
        """
        return self.pileup

    def get_pileup_length(self):

        """
        This function returns the length of the Pileup list in the object.

        INPUT: [None]

        RETURN: [INT] [len(self.pileups[0])]

        POST: [None]
        """
        return len(self.pileup)

    def remove_pileup_positions(self, deletion_list):

        """
        Deletes positions in the Pileup specified in deletion_list, an array
        of integers sorted in descending order.

        INPUT:
            [ARRAY] [deletion_list] - list of positions to delete in descending
                                      order

        RETURN:
            [None]

        POST:
            The specified positions in deletion_list have been removed from
            the self.pileup in the Pileup object.
        """
        for position in deletion_list:
            del self.pileup[position]
        # end for

    def select_pileup_range(self, curr_start, curr_end):

        """
        Ignores all regions of the Pileup before curr_start and after curr_end

        INPUT:
            [int] [curr_start] - current start position. Must be between zero
                                 inclusive and the length of the Pileup
                                 exclusive.
            [int] [curr_end] - current end position. Must be between zero
                               inclusive and the length of the Pileup
                               exclusive.
        RETURN:
            [None]
        POST:
            Positions before curr_start and after curr_end are ignored in the
            Pileup.
        """
        self.pileup = self.pileup[curr_start:curr_end + 1]


class DistanceMatrix(object):

    def __init__(self, pileups):

        """
            [ARRAY] [pileups] - list of pileups - each pileup is
            represented by a list of dictionaries.
        """
        self.pileups = pileups
    # end def

    def get_distance_matrix(self):

        """
        Runs the script, calculating the angular cosine distance function
        between viral quasispecies provided in the pileup.

        Angular Cosine Distance = 2 * ACOS(similarity) / PI

        INPUT:
            [None]

        RETURN:
            Returns a pairwise matrix containing the angular cosine distance
            between all viral quasispecies is returned. The first row and first
            column of the distance matrix contain labels for which quasispecies
            are to be compared in each cell corresponding to the row and column

        POST:
            The internal pileup object is not changed by this function.

        """
        matrix = self.get_similarity_matrix()
        new_matrix = 2 * np.arccos(matrix) / np.pi
        return new_matrix.tolist()
    # end def

    def get_similarity_matrix_as_csv(self, file_list):

        """
        Converts a 2D array (angular cosine distance matrix) to a csv-formatted
        string. Prints out 8 decimal places.

        INPUT:
            [FILE LOCATION TUPLE] [file_list] - files names which represent a
                                                pileup

        RETURN:
            [STRING] [csvOut] CSV representation of a pairwise similarity
            matrix
        """
        matrix = self.get_similarity_matrix()
        return self.__get_matrix_as_csv(matrix, file_list)

    def get_distance_matrix_as_csv(self, file_list):

        """
        Converts a 2D array (angular cosine distance matrix) to a csv-formatted
        string. Prints out 8 decimal places.

        INPUT:
            [FILE LOCATION TUPLE] [file_list] - files names which represent a
                                                pileup

        RETURN:
            [STRING] [csvOut] CSV representation of a pairwise similarity
            matrix
        """
        matrix = self.get_distance_matrix()
        return self.__get_matrix_as_csv(matrix, file_list)

    def __get_matrix_as_csv(self, matrix, file_list):

        """
        Converts a 2D array (cosine similarity matrix) to a csv-formatted
        string. Prints out 8 decimal places.

        INPUT:
            [ARRAY] [matrix] - 2D array (cosine similarity matrix)

            [FILE LOCATION TUPLE] [file_list] - files names which represent a
                                                pileup

        RETURN:
            [STRING] [csvOut] CSV representation of a pairwise similarity
            matrix

        POST:
            [None]

        """
        # (distMatrix[i+1]).insert(0, file_list[i])
        # convert from 2d array to csv formatted string
        files = [file for file in list(file_list)]
        csvOut = 'Quasispecies,' + ','.join(files)
        for row in range(0, len(matrix)):
            csvOut += "\n"
            currElements = ['%.08f' % element for element in matrix[row]]
            csvOut += ','.join([file_list[row]] + currElements)
        # end for
        return csvOut
    # end def

    def get_similarity_matrix(self):

        """
        Runs the script, calculating the cosine similarity function between
        viral quasispecies provided in the pileup.

        Cosine similarity = (u * v) / ( ||u|| * ||v|| )

        INPUT:
            [None]

        RETURN:
            Returns a pairwise matrix containing the cosine similarity function
            between all viral quasispecies is returned. The 1st row and column
            of the similarity matrix contain labels for which quasispecies
            are to be compared in each cell corresponding to the row and column

        POST:
            The internal pileup object is not changed by this function.

        """
        baseList = []
        first = 0  # first position of dictionaries in each pileup_list[i]
        last = len(self.pileups[0])  # position after end position
        for num in range(0, len(self.pileups)):
            baseList.append([self.pileups[num][dict].get(base, 0)
                            for dict in range(first, last) for base in BASES])
        # end for
        baseList = np.array(baseList)
        np.set_printoptions(suppress=True)
        # create distance matrix for csv file
        simi_matrix = squareform(1 - pdist(baseList, cosine))
        di = np.diag_indices(len(simi_matrix))
        simi_matrix[di] = 1.0
        simi_matrix = simi_matrix.tolist()
        return simi_matrix
