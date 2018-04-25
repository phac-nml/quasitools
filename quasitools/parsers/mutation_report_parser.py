"""
Copyright Government of Canada 2017

Written by: Camy Tran, National Microbiology Laboratory,
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

from quasitools.aa_variant import AAVariant, AAVariantCollection
from collections import defaultdict


# Parses an hmcf file into a mutation list containing tuples
# in the form [pos, alt]
def parse_mutations_from_hmcf(hmcf_file):
    mutation_list = defaultdict(list)
    # mutation_list = []
    

    with open(hmcf_file, "r") as input:
        for line in input:
            if line[0] != "#":
                (chrom, gene, id, ref, pos, alt,
                filter, freq, coverage, info) = \
                    line.rstrip().split("\t")
                
                mutation_list[pos].append(alt)

                # wc, mc, mcf, cat, srvl = info.rstrip().split(';')
                # print("%s %s\n" % (pos, alt))
               # mutation_list.append([pos, alt])
               # mutation_list.append( ("%s%s" % (pos, alt)) )

    return mutation_list

                
        

            


