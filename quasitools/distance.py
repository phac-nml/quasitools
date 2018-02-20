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

import sys
import numpy as np

# Quasitools parsers:
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.reference_parser import parse_references_from_fasta
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.spatial.distance import cosine

class Distance(object):

    BASES = ['A', 'C', 'T', 'G']

    #def __init__(self):

    #end def

    def construct_pileup(self, viral_files, reference_loc):

        """
        Creates a array of pileups (which are arrays of dictionaries)

        INPUT:
            viral_files - files names which represent a pileup

            reference_loc - location of the reference file
        RETURN:
            pileup_list - list of pileups
        POST:
            Pileup list is constructed.
        """

        # Build the reference object.
        references = parse_references_from_fasta(reference_loc)

        # Iterate over each reference in the reference object.
        for reference in references:
            mrcList = []
            for bam in viral_files:
                mrcList.append(parse_mapped_reads_from_bam(reference, bam))


            pileup_list = []
            for mrc in mrcList:
                pileup_list.append(mrc.pileup(indels=True))


            """
            pileup_list = [[{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1, 'C': 1},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test1
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1, 'C': 1000000},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test2
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1000000, 'C': 1},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test3
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1, 'C': 1},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test4
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1, 'C': 1000000},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test5
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1000000, 'C': 1},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test6
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1, 'T': 1000000, 'C': 1000000},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test7
                  [{'A': 1, 'T': 1, 'C': 1, 'G': 1}, {'A': 1000000, 'T': 1000000, 'C': 1000000},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}]] #test8

            pileup_list2 = [[{'A': 3, 'T': 3, 'C': 2, 'G': 4}, {'T': 6, 'C': 6},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test9
                  [{'A': 3, 'T': 3, 'C': 3, 'G': 3}, {'T': 6, 'C': 6},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test10
                  [{'A': 3, 'T': 3, 'C': 2, 'G': 4}, {'A': 1, 'T': 1, 'C': 10},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test11
                  [{'A': 3, 'T': 3, 'C': 2, 'G': 4}, {'C': 12},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test12
                  [{'A': 3, 'T': 3, 'C': 2, 'G': 4}, {'G': 6, 'C': 6},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test13
                  [{'A': 3, 'T': 3, 'C': 2, 'G': 4}, {'T': 6, 'G': 6},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}], #test14
                  [{'A': 3, 'T': 3, 'C': 2, 'G': 4}, {'G': 6, 'A': 6},
                   {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}, {'T': 12}, {'C': 12}, {'G': 12}, {'A': 12}]] #test15
            #pileups are a list of dictionaries containing read count for each base
            """

        return pileup_list
    #end def

    def get_distance(self,  startpos, endpos, pileup_list, normalize):

        """
        Runs the script, calculating the cosine similarity function between viral
        quasispecies in viral_files.

        INPUT:
            startpos - starting base position of reference to be compared
            when calculating distances

            endpos - last base position of reference to be compared when
            calculating distance.

            pileup_list - list of lists of dictionaries - each lists of
            dictionaries represents the pileups.

            normalize - determine whether to normalize data or not
        RETURN:
            A pairwise matrix containing the evolutionary distance between all
            viral quasispecies is returned.

        POST:
            [None]

        """

        #TODO - Update normalization to work with multiple input files
        if normalize:
            self.normalize_sum_to_one(pileup_list)

        baseList = []
        first = 0 #first position of dictionaries in each pileup_list[i]
        if startpos != -1:
            first = startpos
        #end if
        for num in range(0, len(pileup_list)):
            last = len(pileup_list[num]) #last pos of dicts in each pileup_list[i]
            if endpos != -1:
                last = endpos + 1
            #end if
            baseList.append([pileup_list[num][dict].get(base, 0)
            for dict in range(first, last) for base in self.BASES])
        #end for
        baseList = np.array(baseList)
        np.set_printoptions(linewidth=140, suppress="True", precision=6)

        distMatrix = squareform(pdist(baseList, self.get_cosine_similarity))
        print(distMatrix)
        return distMatrix

        print("Complete!")
    #end def

    def get_cosine_similarity(self, quasi1, quasi2):
        """
        Calculates the angular cosine distance. Calls cosine(u,v,w=None) from
        scipy.spatial.distance.cosine and returns the angular cosine distance.

        INPUT:
            [NUMPY ARRAY] [quasi1] (a list of read counts per base for the first viral
             quasispecies to be compared)
            [NUMPY ARRAY] [quasi2] (a list of read counts per base for the second viral
             quasispecies to be compared)

        RETURN:
            Returns the cosine similarity.

        POST:
            None.
        """
        #Determine total A, C, T, G pileup1 and pileup2
        return 1 - cosine(quasi1, quasi2)
        #commented out below:
        #computes cosine distance then converts to angular cosine distance
        #return 1 - ( np.arccos( 1 - cosine(quasi1, quasi2) ) / np.pi )
    #end def

    def normalize_sum_to_one(self, pileup_list):

        """
        This function calculates the mean of the read counts for each groupings of four
        bases (A, C, T, G) and subtracts this mean from the individual read counts. this
        prevents large read counts for a base from inflating the cosine simularity
        calculation.

        INPUT:

        [ARRAY] [pileup_list] (a list of dictionaries containing read counts for
        nucleotides at each position corresponding to the reference)

        RETURN: NONE

        POST: pileup_list data is normalized.
        """

        for num in range(0, len(pileup_list)):
            #get the mean for sample one
            mean = 0
            for i in range(0, len(pileup_list[num])):
                total = np.sum([pileup_list[num][i].get(base,0) for base in self.BASES])
                #normalize the data for all samples (centered cosine similarity)
                pileup_list[num][i] = {
                key: (value / total) for key, value in pileup_list[num][i].items() }
        '''
        for num in range(0, len(pileup_list)):
            #get the mean for sample one
            mean = 0
            for i in range(0, len(pileup_list[num])):
                mean = np.sum([pileup_list[num][i].get(base,0) for base in BASES])/4
                #normalize the data for all samples (centered cosine similarity)
                pileup_list[num][i] = {
                key: (value - mean) for key, value in pileup_list[num][i].items() }
                '''
    #end def
