import re
import Bio.SeqIO.FastaIO
from quasitools.reference import Reference

class References(object):
    def __init__(self, path=None):
        self.references = {}
        self.path = path

    @classmethod
    def from_fasta(cls, fasta):
        """Build the References object from a fasta file.

        >>> rs = References.from_fasta('tests/data/ref1.fasta')
        >>> print(rs.references.keys())
        dict_keys(['ref1'])
        >>> print(rs.references['ref1'].seq)
        gattaca
        """
        obj = cls(fasta)

        handle = open(fasta)

        for header, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(handle):
            name = re.search("(\S+)", header).group(0)
            obj.references[name] = Reference(name, seq)

        return obj

    def sub_seq(self, id, start, end):
        return self.references[id].sub_seq(start, end)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
