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

class MappedRead(object):
    def __init__(self, seq_id, query_start, query_end, differences, ref_start, ref_end, overlap, identity, strand):
        self.seq_id      = seq_id
        self.query_start = query_start
        self.query_end   = query_end
        self.differences = differences
        self.ref_start   = ref_start
        self.ref_end     = ref_end
        self.overlap     = overlap
        self.identity    = identity
        self.strand      = strand

    def query_length(self):
        """Calculate and return the length of read mapped to the reference.
        
        >>> mr = MappedRead('read1', 10, 300, {'110': 'g', '300': 'tagt'}, 100, 390, 97, 99.3127, '+')
        >>> print(mr.query_length())
        291
        """
        return self.query_end + 1 - self.query_start

    def codon_start(self, frame):        
        """Calculate and return the first position of the first codon for the given frame in reference base position.
        
        >>> mr = MappedRead('read1', 10, 300, {'110': 'g', '300': 'tagt'}, 100, 390, 96.6667, 99.3103, '+')
        >>> print(mr.codon_start(0))
        102
        >>> print(mr.codon_start(1))
        100
        >>> print(mr.codon_start(2))
        101
        """
        codon_start = self.ref_start

        while codon_start % 3 != frame:
            codon_start += 1

        return codon_start

    def codon_end(self, frame):
        """Calculate and return the last position of the last codon for the given frame in reference base position.
        
        >>> mr = MappedRead('read1', 10, 300, {'110': 'g', '300': 'tagt'}, 100, 390, 96.6667, 99.3103, '+')
        >>> print(mr.codon_end(0))
        389
        >>> print(mr.codon_end(1))
        390
        >>> print(mr.codon_end(2))
        388
        """
        codon_end = self.ref_end

        while (codon_end-2) % 3 != frame:
            codon_end -= 1

        return codon_end

if __name__ == '__main__':
    import doctest
    doctest.testmod()
