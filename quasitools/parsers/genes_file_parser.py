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


def parse_genes_file(genes_file, ref_chrom):
    """Parses a genes file (bed file) and returns a dictionary of genes.
    Each gene is a dictionary that contains the gene's starting and ending
    positions as well as its frame."""

    genes = defaultdict(lambda: defaultdict(dict))

    with open(genes_file, "r") as f:
        for line in f:
            chrom, start, end, name = line.rstrip().split("\t")

            if chrom != ref_chrom:
                raise ValueError("Chrom in genes bed file doesn't match "
                                 "reference")

            genes[name]["chrom"] = chrom
            genes[name]["start"] = int(start)
            genes[name]["end"] = int(end)
            genes[name]["frame"] = int(start) % 3

    return genes
