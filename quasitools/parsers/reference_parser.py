"""
Copyright Government of Canada 2015

Written by: Eric Enns, National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import re
import Bio.SeqIO.FastaIO
from quasitools.reference import Reference

def parse_reference_fasta(fasta):
    """Build the References object from a fasta file.

    >>> rs = parse_reference_fasta('tests/data/ref1.fasta')
    >>> print(len(rs))
    1
    >>> print(rs[0].seq)
    AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT
    """
    references = ()

    handle = open(fasta)

    for header, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(handle):
        name = re.search("(\S+)", header).group(0)
        references += (Reference(name, seq),)

    return references

if __name__ == '__main__':
    import doctest
    doctest.testmod()
