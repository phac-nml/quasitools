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
            pileup_list - list of pileups (list of dictionaries containing
            read counts for each base)
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
            #end for

            #pileups are a
            for num in range(0, len(mrcList)):
                if len(pileup_list) < (num + 1):
                    pileup_list.append(mrcList[num].pileup(indels=True))
                else:
                    pileup_list[num] += mrcList[num].pileup(indels=True)
                #end if
            #end for

        return pileup_list
    #end def

    def get_distance_as_csv(self,  startpos, endpos, pileup_list, viral_files, normalize):

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

            viral_files - files names which represent a pileup

            normalize - determine whether to normalize data or not

        RETURN:
            Returns a csv string representation of a pairwise matrix containing
             the evolutionary distance between all viral quasispecies is
            returned. The first row and first column of the distance matrix
            contain labels for which quasispecies are to be compared in each
            cell corresponding to the row and column.

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

        #create distance matrix for csv file
        distMatrix = squareform(pdist(baseList, self.get_cosine_similarity)).tolist()
        distMatrix.insert(0, ["Quasispecies"] + list(viral_files))
        for i in range(0,len(viral_files)):
            (distMatrix[i+1]).insert(0, viral_files[i])

        #convert from 2d array to csv formatted string
        csvOut=""
        for arr in distMatrix:
            if csvOut != "":
                csvOut += "\n"
            #end if
            csvOut += ','.join(['%s' % element for element in arr])
        #end for
        return csvOut
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
