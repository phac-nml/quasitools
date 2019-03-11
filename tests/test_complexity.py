import pytest
import os
import quasitools.haplotype as haplotype

same_sequence_list = ["AAAAAAA", "AAAAAAA", "AAAAAAA"]
empty_list = [""]
uneven_list = ["", "TAG", ""]

    

@pytest.fixture()
def haplotypes_same_sequence():
        
    haplotypes = []
    for sequence in same_sequence_list:           
        haplotypes.append(haplotype.Haplotype(sequence))

    return haplotypes

def test_consensus_same_haplotypes(haplotypes_same_sequence):

    result = haplotypes_same_sequence
    haplotype.build_consensus_from_haplotypes(haplotypes_same_sequence)
        
    assert result


        



    
