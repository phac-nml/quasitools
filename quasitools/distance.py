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

from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.spatial.distance import cosine


class DistanceMatrix(object):

    def __init__(self, pileups, file_list):

        """
            [ARRAY] [pileups] - Two dimensional numerical array that represents
            a list of pileups. Every row represents a pileup and every four
            values in each row represents the base counts for a particular
            position for the pileup.

            [FILE LOCATION TUPLE] [file_list] - files names which represent a
                                                pileup
        """
        self.pileups = pileups
        self.file_list = file_list
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

    def get_similarity_matrix_as_csv(self):

        """
        Converts a 2D array (angular cosine distance matrix) to a csv-formatted
        string. Prints out 8 decimal places.

        INPUT:
            [None]

        RETURN:
            [STRING] [csvOut] CSV representation of a pairwise similarity
            matrix
        """
        matrix = self.get_similarity_matrix()
        return self.__get_matrix_as_csv(matrix)

    def get_distance_matrix_as_csv(self):

        """
        Converts a 2D array (angular cosine distance matrix) to a csv-formatted
        string. Prints out 8 decimal places.

        INPUT:
            [None]

        RETURN:
            [STRING] [csvOut] CSV representation of a pairwise similarity
            matrix
        """
        matrix = self.get_distance_matrix()
        return self.__get_matrix_as_csv(matrix)

    def __get_matrix_as_csv(self, matrix):

        """
        Converts a 2D array (cosine similarity matrix) to a csv-formatted
        string. Prints out 8 decimal places.

        INPUT:
            [ARRAY] [matrix] - 2D array (cosine similarity matrix)

        RETURN:
            [STRING] [csvOut] CSV representation of a pairwise similarity
            matrix

        POST:
            [None]

        """
        # (distMatrix[i+1]).insert(0, file_list[i])
        # convert from 2d array to csv formatted string
        files = [file for file in list(self.file_list)]
        csvOut = 'Quasispecies,' + ','.join(files)
        for row in range(0, len(matrix)):
            csvOut += "\n"
            currElements = ['%.08f' % element for element in matrix[row]]
            csvOut += ','.join([self.file_list[row]] + currElements)
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
        baseList = np.array(self.pileups)
        # create distance matrix for csv file
        simi_matrix = squareform(1 - pdist(baseList, cosine))
        di = np.diag_indices(len(simi_matrix))
        simi_matrix[di] = 1.0
        simi_matrix = simi_matrix.tolist()
        return simi_matrix
