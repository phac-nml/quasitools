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


class Distance(object):

    # def __init__(self):

    # end def

    def construct_pileup(self, viral_files, reference_loc):

        """
        Creates a array of pileups (which are arrays of dictionaries)

        INPUT:
            [TUPLE] [viral_files] - files names which represent a pileup

            [STRING] [reference_loc] - location of the reference file
        RETURN:
            [ARRAY] [pileup_list] - list of pileups (list of dictionaries
            containing read counts for each base)
        POST:
            Pileup list is constructed.
        """

        # Build the reference object.
        references = parse_references_from_fasta(reference_loc)
        pileup_list = []
        # Iterate over each reference in the reference object.
        for reference in references:
            mrcList = []
            for bam in viral_files:
                mrcList.append(parse_mapped_reads_from_bam(reference, bam))
            # end for

            # pileups are a list of dictionaries. The pileup list contains
            # a list of pileups for each mapped read.
            for num in range(0, len(mrcList)):
                if len(pileup_list) < (num + 1):
                    # add a new pileup to the pileup list
                    pileup_list.append(mrcList[num].pileup(indels=True))
                else:
                    # add the positions in the pileup for this mapped read
                    # that have not been added yet to this pileup (e.g. if you
                    # have multiple references in the reference file)
                    pileup_list[num] += mrcList[num].pileup(indels=True)
                # end if
            # end for
        return pileup_list
    # end def

    def truncate_output(self, pileup_list):
        """
        Deletes sections of the pileup for all pileups in the pileup list
        where there is no coverage (all four bases - A, C, T, and G) are
        absent.

        INPUT:
            [ARRAY] [pileup_list] - list of pileups - each pileup is
            represented by a list of dictionaries.

        RETURN:
            [None]

        POST:
            The pileups are truncated (sections of the pileup where there
            is no coverage are deleted from all pileups in the pileup list.
        """
        deletion_list = []
        if len(pileup_list) > 0:
            # iterate through every position in reference
            if len(pileup_list[0]) > 0:
                for position in range(0, len(pileup_list[0])):
                    # if any pileup at the current position is empty,
                    # add to deletion_list
                    if any((pileup[position] == {} for pileup in pileup_list)):
                        deletion_list.insert(0, position)
                    elif any(sum(pileup[position].values()) == 0
                                 for pileup in pileup_list):
                        deletion_list.insert(0, position)
                # end for
            # end if
        # end if
        for position in deletion_list:
            for pileup in pileup_list:
                del pileup[position]
            # end for
        # end for
        return pileup_list
    # end def

    def get_distance_matrix(self, pileup_list, normalized, startpos=None,
                            endpos=None):

        """
        Runs the script, calculating the cosine similarity function between
        viral quasispecies in viral_files.

        INPUT:
            [ARRAY] [pileup_list] - list of pileups - each pileup is
            represented by a list of dictionaries.

            [BOOL] [normalized] - determine whether to normalize data or not

            [INT] [startpos] -starting base position of reference to be
            compared when calculating cosine similarity.

            [INT] [endpos] - last base position of reference to be compared
            when calculating cosine similarity.

        RETURN:
            Returns a pairwise matrix containing the cosine similarity function
            between all viral quasispecies is returned. The first row and first
            column of the distance matrix contain labels for which quasispecies
            are to be compared in each cell corresponding to the row and column

        POST:
            When and if the data is normalized a new pileup_list object is used
            internally, so that the object whose reference was passed to this
            function is not affected.

        """
        if normalized:
            pileup_list = self.normalize_sum_to_one(pileup_list)
        # end def

        baseList = []
        first = 0  # first position of dictionaries in each pileup_list[i]
        if startpos is not None:
            first = startpos
        # end if
        for num in range(0, len(pileup_list)):
            # last pos of dicts in each pileup_list[i]
            last = len(pileup_list[num])
            if endpos is not None and (endpos + 1) <= last:
                last = endpos + 1
            # end if
            baseList.append([pileup_list[num][dict].get(base, 0)
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

    # end def

    def convert_distance_to_csv(self, matrix, file_list):

        """
        Converts a 2D array (cosine similarity matrix) to a csv-formatted
        string

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

    def normalize_sum_to_one(self, pileup_list):

        """
        This function calculates the mean of the read counts for each groupings
        of four bases (A, C, T, G) and subtracts this mean from the individual
        read counts. this prevents large read counts for a base from inflating
        the cosine simularity calculation.

        INPUT:

        [ARRAY] [pileup_list] (a list of dictionaries containing read counts
        for nucleotides at each position corresponding to the reference)

        RETURN: [ARRAY] [pileup_list] (the normalized list of pileups)

        POST: pileup_list data is normalized.
        """
        new_list = []
        for num in range(0, len(pileup_list)):
            new_list.append([])
            for i in range(0, len(pileup_list[num])):
                curr_pos = [pileup_list[num][i].get(base, 0) for base in BASES]
                total = float(np.sum(curr_pos))
                # normalize the data for all samples
                if total > 0:
                    new_list[num].append(
                        {key: (float(value) / total)
                            for (key, value) in pileup_list[num][i].items()})
                else:
                    new_list[num].append(
                        {key: 0
                            for (key, value) in pileup_list[num][i].items()})
                # end if
        return new_list
        '''
        # get the mean for sample one
        mean = 0
        for num in range(0, len(pileup_list)):
            # get the mean for sample one
            mean = 0
            for i in range(0, len(pileup_list[num])):
                mean = np.sum([pileup_list[num][i].get(base,0)
                           for base in BASES])/4
                # normalize the data for all samples
                # (centered cosine similarity)
                pileup_list[num][i] = {key: (value - mean)
                for key, value in pileup_list[num][i].items() }
                '''
    # end def
