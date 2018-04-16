"""
Copyright Government of Canada 2018

Written by: Eric Marinier and Matthew Fogel, National Microbiology Laboratory,
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

import Bio.SeqIO
from Bio.Seq import Seq

TRIMMING = "trimming"
MASKING = "masking"

# GLOBALS

MASK_CHARACTER = "N"
MINIMUM_QUALITY = "minimum_quality"
LENGTH_CUTOFF = "length_cutoff"
MEDIAN_CUTOFF = "median_cutoff"
MEAN_CUTOFF = "mean_cutoff"

# FILTERING SPECIFICATIONS

# TODO

"""
# =============================================================================

GET MEDIAN SCORE
-------------------


PURPOSE
-------

Returns the median score of a sequence record.


INPUT
-----

[BIOPYTHON SEQRECORD] [read]
    The read where the median score will be retrieved from.

Returns
----

[INT] [median_score]

# =============================================================================
"""


def get_median_score(read):
    length = len(read.seq)
    last_pos = length - 1

    scores = list(read.letter_annotations['phred_quality'])

    if length % 2 == 0:
        median_score = ((scores[int(last_pos // 2)] +
                         scores[int(last_pos // 2) + 1]) / 2)
    else:
        median_score = scores[int(last_pos // 2)]

    return median_score


"""
# =============================================================================

GET MEAN SCORE
-------------------


PURPOSE
-------

Returns the mean score of a sequence record.


INPUT
-----

[BIOPYTHON SEQRECORD] [read]
    The read where the mean score will be retrieved from.

Returns
----

[INT] [mean_score]

# =============================================================================
"""


def get_mean_score(read):

    mean_score = (float(sum(read.letter_annotations['phred_quality'])) /
                  float(len(read.letter_annotations['phred_quality'])))

    return mean_score


"""
# =============================================================================

TRIM READ
---------


PURPOSE
-------

Trims the read in an iterative manner until it meets all the filtering
criteria. The specifications of the filtering parameters can be found in this
file.


INPUT
-----

[BIOPYTHON SEQRECORD] [read]
    The read to iteratively trim until it meets filtering criteria.

[(FILTER -> VALUE) DICTIONARY] [filters]
    The filtering critiera, as a dictionary of (filter, value) pairs, which
    will be tested against the read after each iteration of read trimming.

POST
----

The read will first be evaluated against the filtering criteria. If the read
passes the filtering criteria, then the read will be returned without
modification. However, if the read does not meet the filtering criteria, it
will be trimmed iteratively.

In the latter case, the read will be trimmed iteratively, base by base, from
the tail end of the read, starting at position [len(read) - 1], until either
the read meets the filtering criteria, or the read is too short.

The trimmed and modified read will still be returned if it fails to be improved
with iterative trimming.

# =============================================================================
"""


def trim_read(read, filters):

    length = len(read.seq)
    length_cutoff = filters.get(LENGTH_CUTOFF)

    # while read has not passed all filters and is <= the length cutoff,
    # iteratively trim the read
    while not passes_filters(read) and length >= length_cutoff:
        read = read[:-1]
        length = len(read.seq)

    return read


"""
# =============================================================================

MASK_READ
---------


PURPOSE
-------

Masks the nucleotide of all low quality positions in the read with an "N".


INPUT
-----

[BIOPYTHON SEQRECORD] [read]
    The read to mask low quality positions within.

[(FILTER -> VALUE) DICTIONARY] [filters]
    The filtering critiera, as a dictionary of (filter, value) pairs. The
    minimum quality score will be taken from this dictionary as the value of
    the MINIMUM_QUALITY key.


POST
----

The nucleotide positions in the passed read will be masked with "N"s if their
PHRED quality score is below the minimum.

# =============================================================================
"""


def mask_read(read, filters):

    scores = list(read.letter_annotations['phred_quality'])
    minimum = int(filters.get(MINIMUM_QUALITY))

    # Check every quality score:
    for i in range(0, len(scores)):

        score = int(scores[i])

        # Is the score too low?
        if score < minimum:

            # Mask the base at this position:
            sequence = str(read.seq)
            sequence = sequence[:i] + MASK_CHARACTER + sequence[i + 1:]
            read.seq = Seq(sequence)

    return


"""
# =============================================================================

PASSES FILTERS
--------------


PURPOSE
-------

Determines whether or not the read passes all the filtering criteria.


INPUT
-----

[BIOPYTHON SEQRECORD] [read]
    The read that will be evaluated using the filtering criteria.

[(FILTER -> VALUE) DICTIONARY] [filters]
    The filtering critiera, as a dictionary of (filter, value) pairs, all of
    which will be tested against the read.

RETURN
------

[BOOL] [result]
    TRUE: the the read passes all the filtering criteria.
    FALSE: the read fails at least one of the filtering criteria.

# =============================================================================
"""


def passes_filters(read, filters):

    length_cutoff = filters.get(LENGTH_CUTOFF)
    median_cutoff = filters.get(MEDIAN_CUTOFF)
    mean_cutoff = filters.get(MEAN_CUTOFF)
    filter_ns = filters.get('ns')

    if length_cutoff and len(read.seq) < length_cutoff:
        return False

    if median_cutoff and get_median_score(read) < median_cutoff:
        return False

    if mean_cutoff and get_mean_score(read) < mean_cutoff:
        return False

    if filter_ns and 'n' in read.seq.lower():
        return False

    return True


"""
# =============================================================================

FILTER JOBS
-----------


PURPOSE
-------

Filters reads according to a variety of filtering criteria.


INPUT
-----

[FILE LOCATION] [reads_location]
    The file location of a FASTQ-encoded reads file.

[FILE LOCATION] [output_location]
    The output location of the filtered reads.

[(FILTER -> VALUE) DICTIONARY] [filters]
    The filtering critiera, as a dictionary of (filter, value) pairs, all of
    which will be tested against the read.


POST
----

The reads that pass the filtering criteria will be written to the
[output_location].

# =============================================================================
"""


def filter_reads(reads_location, output_location, filters):

    output_file = open(output_location, "w+")
    reads = Bio.SeqIO.parse(reads_location, "fastq")

    for read in reads:

        if filters.get(TRIMMING):

            trim_read(read, filters)

        if passes_filters(read, filters):

            if filters.get(MASKING):

                mask_read(read, filters)

            Bio.SeqIO.write(read, output_file, "fastq")

    output_file.close()
