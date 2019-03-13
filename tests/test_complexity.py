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

@pytest.fixture()
def haplotypes_empty_list():

    haplotypes = []
    for sequence in empty_list:
        haplotypes.append(haplotype.Haplotype(sequence))
    
    return haplotypes

def test_consensus_same_haplotypes(haplotypes_same_sequence):

    result = haplotype.build_consensus_from_haplotypes(haplotypes_same_sequence)        
    assert result == "AAAAAAA"

def test_consensus_from_empty_list(haplotypes_empty_list):

    
    result = haplotype.build_consensus_from_haplotypes(haplotypes_empty_list)
    assert result == ""




#TODO Get this function running.
#def end_to_end(bam_location, reference_file):
    
    # call the entire program.
