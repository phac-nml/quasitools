import sys

# Quasitools parsers:
from quasitools.parsers.mapped_read_parser import parse_mapped_reads_from_bam
from quasitools.parsers.reference_parser import parse_references_from_fasta

# Locations of the reference and BAM files:
referenceLocation = "data/hxb2_pol.fas"
bamLocation = "BC07-303.bam"

# Build the reference object.
references = parse_references_from_fasta(referenceLocation)

# Iterate over each reference in the reference object.
for reference in references:
    mrc = parse_mapped_reads_from_bam(reference, bamLocation)
    pileup = mrc.pileup(indels=True)

    # Iterate over each position in the reference.
    #for i in range(0, len(pileup)):
    for i in range(0, 10):
        sys.stdout.write(str(pileup[i]))

print("\nDone!")

