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
from quasitools.parsers.genes_file_parser import parse_genes_file

class TestGenesFileParser:
    def test_valid_genes_file(self):
        """Tests to make sure that valid genes files (bed files) are parsed
        properly.
        """

        # Create a valid genes file
        valid_genes_file = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), "data", "valid_genes_file.bed")

        ref_name = "ref1"

        genes = {"gene1": {"start": 0, "end": 100},
                 "gene 2": {"start": 101, "end": 200},  # Spaces are allowed in the gene name
                 "gene3": {"start": 201, "end": 300}}

        with open(valid_genes_file, "w+") as f:
            for gene in genes:
                f.write("%s\t%s\t%s\t%s\n" % (ref_name, genes[gene]["start"],
                                            genes[gene]["end"], gene))

        parsed_genes = parse_genes_file(valid_genes_file, ref_name)

        for gene in parsed_genes:
            assert gene in genes
            assert parsed_genes[gene]["start"] == genes[gene]["start"]
            assert parsed_genes[gene]["end"] == genes[gene]["end"]
            assert parsed_genes[gene]["frame"] == genes[gene]["start"] % 3

        os.remove(valid_genes_file)

    def test_invalid_genes_file(self):
        """Tests to make sure that an exception is raised when attempting to
        parse an invalid genes file.
        """

        # Create an invalid genes file
        invalid_genes_file = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), "data", "invalid_genes_file.bed")

        ref_name = "ref1"

        with open(invalid_genes_file, "w+") as f:
            f.write("%s\t0\t100\t0\n" % ref_name)
            # Add a genes reference name that doesn't match
            # This should raise a ValueError
            f.write("different_reference\t101\t200\t2")

        with pytest.raises(ValueError):
            parse_genes_file(invalid_genes_file, ref_name)

        os.remove(invalid_genes_file)
