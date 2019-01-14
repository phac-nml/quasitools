"""
Copyright Government of Canada 2017

Written by: Cole Peters, Eric Chubaty, National Microbiology Laboratory,
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

from collections import defaultdict

# BED FILE POSITIONS:
CHROM = 0
START = 1
END = 2
NAME = 3


def parse_BED4_file(BED4_file, ref_chrom):
    """Parses a BED4+ file (BED file with at least 4 columns) and returns
    a dictionary of genes. Each gene is a dictionary that contains the gene's
    starting and ending positions as well as its frame."""

    genes = defaultdict(lambda: defaultdict(dict))

    with open(BED4_file, "r") as f:
        for line in f:

            tokens = line.rstrip().split("\t")
            chrom = tokens[CHROM]
            start = tokens[START]
            end = tokens[END]
            name = tokens[NAME]
            # all other tokens are ignored

            if chrom != ref_chrom:
                raise ValueError("Chrom in genes bed file doesn't match "
                                 "reference")

            genes[name]["chrom"] = chrom
            genes[name]["start"] = int(start)
            genes[name]["end"] = int(end)
            genes[name]["frame"] = int(start) % 3

    return genes
