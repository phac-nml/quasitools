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

class Reference(object):
    def __init__(self, name, seq, path=None):
        self.name = name
        self.seq = seq
        self.path = path

    @classmethod
    def from_fasta(cls, fasta):
        """Build the Reference object from a fasta.

        >>> r = Reference.from_fasta('tests/data/ref1.fasta')
        >>> print(r.name)
        ref1
        >>> print(r.seq)
        gattaca
        """
        lines = [line.rstrip('\n') for line in open(fasta)]

        header = lines.pop(0)

        if not header.startswith('>'):
            raise Exception('Fasta reference is missing header line')
        elif re.match(""">(\S+)""", header) is None:
            raise Exception('Fasta header is malformed')

        name = re.findall(""">(\S+)""", header)[0]

        seq_lines = []

        for line in lines:
            if lines[0].startswith('>'):
                raise Exception('Fasta reference contains more than one entry')
            seq_lines.append(line)

        seq = ''.join(seq_lines)

        obj = cls(name, seq, fasta)

        return obj

    def sub_seq(self, start, end):
        """Returns a portion of the sequence

        >>> r = Reference('gattaca', 'gattaca')
        >>> print(r.sub_seq(1,5))
        attac
        >>> print(r.sub_seq(1,1))
        a
        """
        return self.seq[start:end+1]

if __name__ == '__main__':
    import doctest
    doctest.testmod()
