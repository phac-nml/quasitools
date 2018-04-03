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
PADDING = '-'


class Pileup_List(object):

    def __init__(self, pileup_list):

        """
        Creates a array of pileups (which are arrays of dictionaries)

        INPUT:
            [ARRAY] [pileup_list] - list of pileups - each pileup is
                                    represented by a list of dictionaries.
        RETURN:
            [None]
        POST:
            Pileup list is constructed.
        """
        self.pileup_list = pileup_list

    def all_have_coverage(self, position):
        """
        Determines whether all pileups in the pileup list have coverage at
        the present position.

        INPUT:
            [INT] [position] - position in each pileup to check for coverage

        RETURN:
            Returns BOOL true if all pileups have coverage at the position or
            false if at least one pileup does not have coverage at the position

        POST:
            None

        """

        if any((np.sum([pileup[position].get(base, 0) for base in BASES]) == 0
                or pileup[position] == {}) for pileup in self.pileup_list):
            return False
        else:
            return True
        # end if

    # end def

    @staticmethod
    def construct_array_of_pileups(viral_files, reference_loc):

        """
        Creates a array of pileups (which are arrays of dictionaries)
        INPUT:
            [TUPLE] [viral_files] - files names which represent a pileup

            [STRING] [reference_loc] - location of the reference file
        RETURN:
            [Pileup_List] - a new object containing a list of pileups
            (list of dictionaries containing read counts for each base)
        POST:
            Pileup list object is constructed.
        """

        # Build the reference object.
        references = parse_references_from_fasta(reference_loc)
        new_pileup_list = []
        # Iterate over each reference in the reference object.
        for reference in references:
            mrcList = []
            for bam in viral_files:
                mrcList.append(parse_mapped_reads_from_bam(reference, bam))
            # end for

            # pileups are a list of dictionaries. The pileup list contains
            # a list of pileups for each mapped read.
            for num in range(0, len(mrcList)):
                if len(new_pileup_list) < (num + 1):
                    # add a new pileup to the pileup list
                    new_pileup_list.append(mrcList[num].pileup(indels=True))
                else:
                    # add the positions in the pileup for this mapped read
                    # that have not been added yet to this pileup (e.g. if you
                    # have multiple references in the reference file)
                    new_pileup_list[num] += mrcList[num].pileup(indels=True)
                # end if
            # end for
        return Pileup_List(new_pileup_list)
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

        POST: The Pileup_Util object's data is normalized.
        """
        new_list = []
        for num in range(0, len(self.pileup_list)):
            new_list.append([])
            for i in range(0, len(self.pileup_list[num])):
                curr_pos = [self.pileup_list[num][i].get(base, 0)
                            for base in BASES]
                total = float(np.sum(curr_pos))
                items = self.pileup_list[num][i].items()
                # normalize the data for all samples
                if total > 0:
                    new_list[num].append(
                        {key: (float(val) / total) for (key, val) in items
                         if key is not '-'})
                else:
                    new_list[num].append(
                        {key: 0 for (key, value) in items if key is not '-'})
                # end if
        self.pileup_list = new_list
    # end def

    def get_pileup(self):

        """
        This function returns the pileup list in the pileup_list object.

        INPUT: [None]

        RETURN: [ARRAY] [pileup_list]

        POST: [None]
        """
        return self.pileup_list

    def get_pileup_length(self):

        """
        This function returns the length of the pileup list in the object.

        INPUT: [None]

        RETURN: [INT] [len(self.pileup_list[0])]

        POST: [None]
        """
        return len(self.pileup_list[0])

    def remove_pileup_positions(self, deletion_list):
        """
        Deletes positions in the pileup specified in deletion_list, an array
        of integers sorted in descending order.

        INPUT:
            [ARRAY] [deletion_list] - list of positions to delete in descending
                                      order

        RETURN:
            [None]

        POST:
            The specified positions in deletion_list have been removed from
            the self.pileup_list in the Pileup_List object.
        """
        for position2 in deletion_list:
            for pileup in self.pileup_list:
                del pileup[position2]
            # end for
        # end for

    def select_pileup_range(self, curr_start, curr_end):
        """
        Ignores all regions of the pileup before curr_start and after curr_end

        INPUT:
            [int] [curr_start] - current start position. Must be between zero
                                 and the length of the pileup.
            [int] [curr_end] - current end position. Must be between zero and
                               the length of the pileup.
        RETURN:
            [None]
        POST:
            Positions before curr_start and after curr_end are ignored in the
            pileup list.
        """
        np_pileup = np.array(self.pileup_list)
        self.pileup_list = (np_pileup[:, curr_start:curr_end + 1]).tolist()

    def truncate_all_output(self):
        """
        Deletes all regions of the pileup for all pileups in the pileup list
        where there is no coverage (all four bases - A, C, T, and G) are
        absent.

        INPUT:
            [None]

        RETURN:
            [None]

        POST:
            The pileups are truncated (sections of the pileup where there
            is no coverage are deleted from all pileups in the pileup list.
        """
        deletion_list = []
        if len(self.pileup_list) > 0 and len(self.pileup_list[0]) > 0:
            # iterate through every position in reference
            for position in range(0, len(self.pileup_list[0])):
                # add pos'n with empty coverage in pileup to deletion_list
                if not self.all_have_coverage(position):
                    deletion_list.insert(0, position)
        self.remove_pileup_positions(deletion_list)
    # end def

    def truncate_output(self):
        """
        Deletes contiguous start and end regions of the pileup for all pileups
        in the pileup list where there is no coverage (all four bases - A, C,
        T, and G) are absent.

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
        num_pos = len(self.pileup_list[0])
        if len(self.pileup_list) > 0 and len(self.pileup_list[0]) > 0:
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
        self.remove_pileup_positions(deletion_list)
    # end def


class DistanceMatrix(object):

    def __init__(self, pileup_list):
        """
            [ARRAY] [pileup_list] - list of pileups - each pileup is
            represented by a list of dictionaries.
        """
        self.pileup_list = pileup_list
    # end def

    def get_angular_cosine_distance_matrix(self, startpos=None, endpos=None):

        """
        Runs the script, calculating the angular cosine distance function
        between viral quasispecies provided in the pileup.

        Angular Cosine Distance = 2 * ACOS(similarity) / PI

        INPUT:
            [INT] [startpos] -starting base position of reference to be
            compared when calculating cosine distance.

            [INT] [endpos] - last base position of reference to be compared
            when calculating cosine distance.

        RETURN:
            Returns a pairwise matrix containing the angular cosine distance
            between all viral quasispecies is returned. The first row and first
            column of the distance matrix contain labels for which quasispecies
            are to be compared in each cell corresponding to the row and column

        POST:
            The internal pileup object is not changed by this function.

        """
        matrix = self.get_cosine_similarity_matrix(startpos, endpos)
        new_matrix = 2 * np.arccos(matrix) / np.pi
        return new_matrix.tolist()
    # end def

    def get_matrix_as_csv(self, matrix, file_list):

        """
        Converts a 2D array (cosine similarity matrix) to a csv-formatted
        string. Prints out 8 decimal places.

        INPUT:
            [ARRAY] [matrix] - 2D array (cosine similarity matrix)

            [TUPLE] [file_list] - files names which represent a pileup

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

    def get_cosine_similarity_matrix(self, startpos=None, endpos=None):

        """
        Runs the script, calculating the cosine similarity function between
        viral quasispecies provided in the pileup.

        Cosine similarity = (u * v) / ( ||u|| * ||v|| )

        INPUT:
            [INT] [startpos] -starting base position of reference to be
            compared when calculating cosine similarity.

            [INT] [endpos] - last base position of reference to be compared
            when calculating cosine similarity.

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
        if startpos is not None:
            first = startpos
        # end if
        for num in range(0, len(self.pileup_list)):
            # last pos of dicts in each pileup_list[i]
            last = len(self.pileup_list[num])
            if endpos is not None and (endpos + 1) <= last:
                last = endpos + 1
            # end if
            baseList.append([self.pileup_list[num][dict].get(base, 0)
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
