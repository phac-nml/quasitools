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

import os
import pytest
from click.testing import CliRunner
from quasitools import cli

# GLOBALS
TEST_PATH = os.path.dirname(os.path.abspath(__file__))
REF = TEST_PATH + '/data/hxb2_pol.fas'
BAM1 = TEST_PATH + '/data/quasi1.bam'
BAM2 = TEST_PATH + '/data/quasi2.bam'
FORWARD = TEST_PATH + '/data/forward.fastq'
REVERSE = TEST_PATH + '/data/reverse.fastq'

@pytest.fixture
def runner():
    return CliRunner()

def test_cli(runner):
    result = runner.invoke(cli.cli)
    assert result.exit_code == 0
    assert not result.exception
    assert result.output.split('\n', 1)[0].strip() == 'Usage: cli [OPTIONS] COMMAND [ARGS]...'

def test_cli_distance(runner):
    """
    test_cli_distance - Checks that if there are errors the correct exit
    code is raised.

    INPUT:
        [None]
    RETURN:
        [None]
    POST:
        [None]
    """

    # tests that are expected to pass

    result = runner.invoke(cli.cli, ['distance', REF])
    assert result.exit_code == 2 # test expected to fail
    assert result.exception # exception should be raised

    result = runner.invoke(cli.cli, ['distance', REF, BAM1])
    assert result.exit_code == 2 # test expected to fail
    assert result.exception # exception should be raised

    result = runner.invoke(cli.cli, ['distance', REF, BAM1, BAM2])
    assert result.exit_code == 0 # test expected to pass
    assert not result.exception # no exception should be raised

    result = runner.invoke(cli.cli, ['distance', REF, BAM1, BAM2, '--normalize', '--output_distance', '--output', '--keep_no_coverage'])
    assert result.exit_code == 0 # test expected to pass
    assert not result.exception # no exception should be raised

    result = runner.invoke(cli.cli, ['distance', REF, BAM1, BAM2, '--normalize', '--output_distance', '-s 1', '-e 1', '--output', '--keep_no_coverage'])
    assert result.exit_code == 0 # test expected to pass
    assert not result.exception # no exception should be raised

    result = runner.invoke(cli.cli, ['distance', REF, BAM1, BAM2, '--normalize', '--output_distance', '-s 1', '-e 2844', '--output', '--keep_no_coverage'])
    assert result.exit_code == 0 # test expected to succeed, UsageError raised
    assert not result.exception # no exception should be raised

    # tests that are expected to fail
    result = runner.invoke(cli.cli, ['distance', REF, BAM1, BAM2, '--normalize', '--output_distance', '-s 0', '-e 2844', '--output', '--keep_no_coverage'])
    assert result.exit_code == 2 # test expected to fail, UsageError raised
    assert result.exception # no exception should be raised

    result = runner.invoke(cli.cli, ['distance', REF, BAM1, BAM2, '--normalize', '--output_distance', '-s -1', '-e 0', '--output', '--keep_no_coverage'])
    assert result.exit_code == 2 # test expected to fail, UsageError raised
    assert result.exception # no exception should be raised

    result = runner.invoke(cli.cli, ['distance', REF, BAM1, BAM2, '--normalize', '--output_distance', '-s 1', '-e 2845', '--output', '--keep_no_coverage'])
    assert result.exit_code == 2 # test expected to fail, UsageError raised
    assert result.exception # no exception should be raised

    result = runner.invoke(cli.cli, ['distance', REF, BAM1, BAM2, '--normalize', '--output_distance', '-s 2845', '-e 2846', '--output', '--keep_no_coverage'])
    assert result.exit_code == 2 # test expected to fail, UsageError raised
    assert result.exception # no exception should be raised

    result = runner.invoke(cli.cli, ['distance', REF, BAM1, BAM2, '--normalize', '--output_distance', '-s 100', '-e 99', '--output', '----keep_no_coverage'])
    assert result.exit_code == 2 # test expected to fail, UsageError raised
    assert result.exception # no exception should be raised

    result = runner.invoke(cli.cli, ['distance', REF, '--normalize', '--output_distance', '-s 1', '-e 2844', '--output', '--keep_no_coverage'])
    assert result.exit_code == 2 # test expected to fail, UsageError raised
    assert result.exception # no exception should be raised

def test_cli_quality(runner):
    """
    test_cli_quality - Adds tests to cli

    INPUT:
        [None]
    RETURN:
        [None]
    POST:
        [None]
    """

    result = runner.invoke(cli.cli, ['quality', FORWARD, REVERSE, '-o', 'test_cli_quality'])
    assert result.exit_code == 0 # test expected to pass
    assert not result.exception # exception should not be raised

    result = runner.invoke(cli.cli, ['quality', FORWARD, REVERSE, '-o', 'test_cli_quality', '-tr', '-mr', '-me', '-n', '-sc', '30', '-lc', '100', '-mq', '30'])
    assert result.exit_code == 0 # test expected to pass
    assert not result.exception # exception should not be raised
