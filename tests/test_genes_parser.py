"""
Copyright Government of Canada 2017

Written by: Cole Peters, National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import pytest
import os
from quasitools.parsers.genes_file_parser import parse_BED4_file


class TestGenesFileParser:
    def test_valid_BED4_file(self):
        """Tests to make sure that valid BED4+ files (bed files with at least 4
        columns) are parsed properly.
        """

        # Create a valid BED4 file
        valid_BED4_file = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), "data", "valid_BED4_file.bed")

        ref_name = "ref1"

        genes = {"gene1": {"start": 0, "end": 100},
                 # Spaces are allowed in the gene name
                 "gene 2": {"start": 101, "end": 200},
                 "gene3": {"start": 201, "end": 300}}

        with open(valid_BED4_file, "w+") as f:
            for gene in genes:
                f.write("%s\t%s\t%s\t%s\n" % (ref_name, genes[gene]["start"],
                                              genes[gene]["end"], gene))

        parsed_genes = parse_BED4_file(valid_BED4_file, ref_name)

        for gene in parsed_genes:
            assert gene in genes
            assert parsed_genes[gene]["start"] == genes[gene]["start"]
            assert parsed_genes[gene]["end"] == genes[gene]["end"]
            assert parsed_genes[gene]["frame"] == genes[gene]["start"] % 3

        os.remove(valid_BED4_file)

    def test_invalid_BED4_file(self):
        """Tests to make sure that an exception is raised when attempting to
        parse an invalid BED4 file.
        """

        # Create an invalid genes file
        invalid_BED4_file = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), "data", "invalid_BED4_file.bed")

        ref_name = "ref1"

        with open(invalid_BED4_file, "w+") as f:
            f.write("%s\t0\t100\t0\n" % ref_name)
            # Add a genes reference name that doesn't match
            # This should raise a ValueError
            f.write("different_reference\t101\t200\t2")

        with pytest.raises(ValueError):
            parse_BED4_file(invalid_BED4_file, ref_name)

        os.remove(invalid_BED4_file)
