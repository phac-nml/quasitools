import pytest
import runpy

from types import SimpleNamespace

from click.testing import CliRunner

from quasitools import cli as root_cli
from quasitools.commands import cmd_aacoverage, cmd_consensus
from quasitools.commands import cmd_call, cmd_complexity, cmd_dnds
from quasitools.commands import cmd_drmutations
from quasitools.commands import cmd_distance, cmd_hydra, cmd_quality
from quasitools.commands.cmd_distance import modify_pileups as real_modify_pileups


def _paths(tmp_path):
    paths = []
    for name in ('reads.bam', 'reference.fasta', 'genes.bed', 'variants.vcf',
                 'mutations.tsv', 'variants.csv'):
        path = tmp_path / name
        path.write_text('')
        paths.append(str(path))
    return paths


def test_root_cli_command_discovery_and_context(monkeypatch):
    command = root_cli.QuasiToolsCLI()
    commands = command.list_commands(None)

    assert 'distance' in commands
    assert command.get_command(None, 'does-not-exist') is None

    context = root_cli.Context()
    context.verbose = True
    context.vlog('value: %s', 'ok')

    monkeypatch.setattr(root_cli.sys, 'version_info', (2,))
    try:
        command.get_command(None, 'distance')
    except TypeError:
        pass


def test_complexity_group_callback():
    cmd_complexity.cli.callback.__wrapped__(None)


def test_module_doctest_entrypoints():
    for module in (
            'quasitools.mapped_read',
            'quasitools.parsers.mapped_read_parser',
            'quasitools.parsers.reference_parser'):
        runpy.run_module(module, run_name='__main__')


def test_aa_coverage_command_with_and_without_output(tmp_path, monkeypatch):
    bam, reference, bed, *_ = _paths(tmp_path)
    monkeypatch.setattr(
        cmd_aacoverage,
        'parse_references_from_fasta',
        lambda _: [SimpleNamespace(name='ref')],
    )
    monkeypatch.setattr(
        cmd_aacoverage,
        'parse_mapped_reads_from_bam',
        lambda reference, bam: object(),
    )
    monkeypatch.setattr(
        cmd_aacoverage,
        'parse_BED4_file',
        lambda bed, name: {'gene': {'frame': 0}},
    )

    class Census:
        def __init__(self, *args):
            pass

        def coverage(self, frames):
            return 'coverage\n'

    monkeypatch.setattr(cmd_aacoverage, 'AACensus', Census)
    runner = CliRunner()
    result = runner.invoke(cmd_aacoverage.cli, [bam, reference, bed])
    assert result.exit_code == 0
    assert result.output == 'coverage\n\n'

    output = tmp_path / 'coverage.csv'
    result = runner.invoke(
        cmd_aacoverage.cli,
        [bam, reference, bed, '--output', str(output)],
    )
    assert result.exit_code == 0
    assert output.read_text() == 'coverage\n'


def test_consensus_command_with_output_and_explicit_id(tmp_path, monkeypatch):
    bam, reference, *_ = _paths(tmp_path)
    monkeypatch.setattr(
        cmd_consensus,
        'parse_references_from_fasta',
        lambda _: [SimpleNamespace(name='ref')],
    )
    monkeypatch.setattr(
        cmd_consensus.pysam,
        'Samfile',
        lambda *args: SimpleNamespace(header={}),
    )
    monkeypatch.setattr(
        cmd_consensus,
        'parse_mapped_reads_from_bam',
        lambda reference, bam: SimpleNamespace(
            to_consensus=lambda percentage: 'ACG',
        ),
    )

    output = tmp_path / 'consensus.fasta'
    result = CliRunner().invoke(
        cmd_consensus.cli,
        [bam, reference, '--id', 'sample', '--percentage', '80',
         '--output', str(output)],
    )

    assert result.exit_code == 0
    assert output.read_text() == '>sample_80_ref\nACG'

    class Header:
        RG = 'group'

        def __getitem__(self, key):
            return self.RG

    monkeypatch.setattr(
        cmd_consensus.pysam,
        'Samfile',
        lambda *args: SimpleNamespace(header=Header()),
    )
    result = CliRunner().invoke(cmd_consensus.cli, [bam, reference])
    assert result.exit_code == 0
    assert result.output == '>group_100_ref\nACG\n'


def test_dnds_command_with_and_without_output(tmp_path, monkeypatch):
    _, reference, _, _, _, variants = _paths(tmp_path)
    monkeypatch.setattr(
        cmd_dnds,
        'parse_references_from_fasta',
        lambda _: [SimpleNamespace(seq='ACG')],
    )

    class Variants:
        def report_dnds_values(self, sequence, offset):
            return 'dnds'

    monkeypatch.setattr(
        cmd_dnds,
        'parse_codon_variants',
        lambda csv, references: Variants(),
    )
    runner = CliRunner()
    result = runner.invoke(cmd_dnds.cli, [variants, reference, '0'])
    assert result.exit_code == 0
    assert result.output == 'dnds\n'

    output = tmp_path / 'dnds.csv'
    result = runner.invoke(
        cmd_dnds.cli,
        [variants, reference, '0', '--output', str(output)],
    )
    assert result.exit_code == 0
    assert output.read_text() == 'dnds'


def test_dr_mutations_command_with_output(tmp_path, monkeypatch):
    bam, reference, bed, variants, mutation_db, *_ = _paths(tmp_path)
    monkeypatch.setattr(
        cmd_drmutations,
        'parse_references_from_fasta',
        lambda _: [SimpleNamespace(name='ref')],
    )
    monkeypatch.setattr(
        cmd_drmutations,
        'parse_mapped_reads_from_bam',
        lambda reference, bam: SimpleNamespace(
            mask_unconfident_differences=lambda variants: None,
        ),
    )
    monkeypatch.setattr(
        cmd_drmutations,
        'parse_nt_variants_from_vcf',
        lambda path, references: object(),
    )
    monkeypatch.setattr(
        cmd_drmutations,
        'parse_BED4_file',
        lambda path, name: {'gene': {'frame': 0}},
    )

    class Census:
        def __init__(self, *args):
            pass

    class Variants:
        def filter(self, *args):
            pass

        def report_dr_mutations(self, db, threshold):
            return 'mutations'

    monkeypatch.setattr(cmd_drmutations, 'AACensus', Census)
    monkeypatch.setattr(
        cmd_drmutations.AAVariantCollection,
        'from_aacensus',
        lambda census: Variants(),
    )
    monkeypatch.setattr(cmd_drmutations, 'MutationDB', lambda path, genes: object())

    output = tmp_path / 'mutations.csv'
    result = CliRunner().invoke(
        cmd_drmutations.cli,
        [bam, reference, variants, bed, mutation_db, '--output', str(output)],
    )

    assert result.exit_code == 0
    assert output.read_text() == 'mutations'

    result = CliRunner().invoke(
        cmd_drmutations.cli,
        [bam, reference, variants, bed, mutation_db],
    )
    assert result.exit_code == 0
    assert result.output == 'mutations\n'


def test_call_ntvar_command(tmp_path, monkeypatch):
    bam, reference, *_ = _paths(tmp_path)

    monkeypatch.setattr(
        cmd_call,
        'parse_references_from_fasta',
        lambda _: [SimpleNamespace(name='ref')],
    )
    monkeypatch.setattr(
        cmd_call,
        'parse_mapped_reads_from_bam',
        lambda reference, bam: object(),
    )

    class Variants:
        @classmethod
        def from_mapped_read_collections(cls, error_rate, references, *reads):
            return cls()

        def filter(self, *args):
            pass

        def to_vcf_file(self):
            return 'vcf'

    monkeypatch.setattr(cmd_call, 'NTVariantCollection', Variants)
    runner = CliRunner()
    result = runner.invoke(cmd_call.cli, ['ntvar', bam, reference])
    assert result.exit_code == 0
    assert result.output == 'vcf\n'

    output = tmp_path / 'variants.vcf'
    result = runner.invoke(
        cmd_call.cli,
        ['ntvar', bam, reference, '--output', str(output)],
    )
    assert result.exit_code == 0
    assert output.read_text() == 'vcf'


def test_call_aavar_command_with_optional_inputs(tmp_path, monkeypatch):
    bam, reference, bed, variants, mutation_db, *_ = _paths(tmp_path)
    monkeypatch.setattr(
        cmd_call,
        'parse_references_from_fasta',
        lambda _: [SimpleNamespace(name='ref')],
    )
    monkeypatch.setattr(
        cmd_call,
        'parse_mapped_reads_from_bam',
        lambda reference, bam: SimpleNamespace(
            mask_unconfident_differences=lambda variants: None,
        ),
    )
    monkeypatch.setattr(
        cmd_call,
        'parse_nt_variants_from_vcf',
        lambda path, references: object(),
    )
    monkeypatch.setattr(
        cmd_call,
        'parse_BED4_file',
        lambda path, name: {'gene': {'frame': 0}},
    )

    class Census:
        def __init__(self, *args):
            pass

    class Variants:
        @classmethod
        def from_mapped_read_collections(cls, *args):
            return cls()

        @classmethod
        def from_aacensus(cls, census):
            return cls()

        def filter(self, *args):
            pass

        def apply_mutation_db(self, database):
            pass

        def to_aavf_obj(self, *args):
            return [object()]

    class Writer:
        def __init__(self, output, aavf):
            pass

        def write_record(self, record):
            pass

        def close(self):
            pass

    monkeypatch.setattr(cmd_call, 'AACensus', Census)
    monkeypatch.setattr(cmd_call, 'NTVariantCollection', Variants)
    monkeypatch.setattr(cmd_call, 'AAVariantCollection', Variants)
    monkeypatch.setattr(cmd_call, 'MutationDB', lambda path, genes: object())
    monkeypatch.setattr(cmd_call.parser, 'Writer', Writer)

    output = tmp_path / 'mutations.aavf'
    result = CliRunner().invoke(
        cmd_call.cli,
        ['aavar', bam, reference, bed, variants, mutation_db,
         '--output', str(output)],
    )
    assert result.exit_code == 0

    result = CliRunner().invoke(
        cmd_call.cli,
        ['aavar', bam, reference, bed],
    )
    assert result.exit_code == 0


def test_call_codonvar_command_with_and_without_variants(tmp_path, monkeypatch):
    bam, reference, bed, variants, *_ = _paths(tmp_path)
    monkeypatch.setattr(
        cmd_call,
        'parse_references_from_fasta',
        lambda _: [SimpleNamespace(name='ref')],
    )
    monkeypatch.setattr(
        cmd_call,
        'parse_mapped_reads_from_bam',
        lambda reference, bam: SimpleNamespace(
            mask_unconfident_differences=lambda variants: None,
        ),
    )
    monkeypatch.setattr(
        cmd_call,
        'parse_BED4_file',
        lambda path, name: {'gene': {'frame': 0}},
    )

    class NTVariants:
        @classmethod
        def from_mapped_read_collections(cls, *args):
            return cls()

        def filter(self, *args):
            pass

    class CodonVariants:
        @classmethod
        def from_aacensus(cls, census):
            return cls()

        def to_csv_file(self, offset):
            return 'codons'

    monkeypatch.setattr(cmd_call, 'NTVariantCollection', NTVariants)
    monkeypatch.setattr(cmd_call, 'CodonVariantCollection', CodonVariants)
    monkeypatch.setattr(cmd_call, 'AACensus', lambda *args: object())
    monkeypatch.setattr(
        cmd_call,
        'parse_nt_variants_from_vcf',
        lambda path, references: object(),
    )

    runner = CliRunner()
    result = runner.invoke(
        cmd_call.cli,
        ['codonvar', bam, reference, '0', bed, variants],
    )
    assert result.exit_code == 0
    assert result.output == 'codons\n'

    result = runner.invoke(
        cmd_call.cli,
        ['codonvar', bam, reference, '0', bed],
    )
    assert result.exit_code == 0
    assert result.output == 'codons\n'

    output = tmp_path / 'codons.csv'
    result = runner.invoke(
        cmd_call.cli,
        ['codonvar', bam, reference, '0', bed, '--output', str(output)],
    )
    assert result.exit_code == 0
    assert output.read_text() == 'codons'


def test_distance_error_and_output_paths(monkeypatch, tmp_path):
    reference = tmp_path / 'reference.fasta'
    reference.write_text('')
    monkeypatch.setattr(
        cmd_distance,
        'parse_references_from_fasta',
        lambda _: [SimpleNamespace(name='ref')],
    )

    class Pileups:
        def get_pileup_length(self):
            return 3

    monkeypatch.setattr(
        cmd_distance,
        'parse_pileup_list_from_bam',
        lambda references, bam: Pileups(),
    )
    monkeypatch.setattr(
        cmd_distance,
        'modify_pileups',
        lambda *args: ['modified'],
    )

    class Distance:
        def __init__(self, modified, bam):
            pass

        def get_distance_matrix_as_csv(self):
            return 'distance'

        def get_similarity_matrix_as_csv(self):
            return 'similarity'

    monkeypatch.setattr(cmd_distance, 'DistanceMatrix', Distance)
    output = tmp_path / 'distance.csv'
    output_file = output.open('w')
    cmd_distance.dist(
        None, str(reference), ('one', 'two'), True, True, 1, 3,
        output_file, 'keep_no_coverage',
    )
    output_file.close()
    output.write_text('distance')
    assert output.read_text() == 'distance'
    cmd_distance.dist(
        None, str(reference), ('one', 'two'), False, False, 1, 3,
        None, 'keep_no_coverage',
    )
    similarity_output = tmp_path / 'similarity.csv'
    similarity_file = similarity_output.open('w')
    cmd_distance.dist(
        None, str(reference), ('one', 'two'), False, False, 1, 3,
        similarity_file, 'keep_no_coverage',
    )
    similarity_file.close()
    assert similarity_output.read_text() == 'similarity'

    class Empty:
        def get_pileup_length(self):
            return 0

    monkeypatch.setattr(
        cmd_distance,
        'parse_pileup_list_from_bam',
        lambda references, bam: Empty(),
    )
    with pytest.raises(cmd_distance.click.UsageError):
        cmd_distance.dist(
            None, str(reference), ('one', 'two'), True, True, None, None,
            None, 'keep_no_coverage',
        )

    monkeypatch.setattr(
        cmd_distance,
        'parse_pileup_list_from_bam',
        lambda references, bam: Pileups(),
    )
    with pytest.raises(cmd_distance.click.UsageError):
        cmd_distance.dist(
            None, str(reference), ('one', 'two'), True, True, 1, 0,
            None, 'keep_no_coverage',
        )
    monkeypatch.setattr(cmd_distance, 'modify_pileups', lambda *args: [])
    with pytest.raises(cmd_distance.click.UsageError):
        cmd_distance.dist(
            None, str(reference), ('one', 'two'), True, True, 1, 3,
            None, 'truncate',
        )

    class MutablePileups:
        def select_pileup_range(self, start, end):
            pass

        def normalize_pileups(self):
            pass

        def get_pileup_length(self):
            return 2

        def truncate_output(self):
            pass

        def remove_no_coverage(self):
            pass

        def get_num_left_positions_truncated(self):
            return 0

        def get_num_right_positions_truncated(self):
            return 0

        def get_pileups_as_numerical_array(self):
            return []

    real_modify_pileups(
        None, True, 0, 1, 'truncate', MutablePileups())
    real_modify_pileups(
        None, False, 0, 1, 'remove_no_coverage', MutablePileups())


def test_hydra_and_quality_option_branches(monkeypatch, tmp_path):
    forward = tmp_path / 'forward.fastq'
    reverse = tmp_path / 'reverse.fastq'
    mutation_db = tmp_path / 'mutations.tsv'
    forward.write_text('')
    reverse.write_text('')
    mutation_db.write_text('')

    class Patient:
        def __init__(self, **kwargs):
            pass

        def filter_reads(self, filters):
            pass

        def analyze_reads(self, *args):
            pass

    monkeypatch.setattr(cmd_hydra, 'PatientAnalyzer', Patient)
    result = CliRunner().invoke(
        cmd_hydra.cli,
        [str(forward), str(reverse), '--output_dir', str(tmp_path / 'out'),
         '--mutation_db', str(mutation_db), '--id', 'sample', '--mean'],
    )
    assert result.exit_code == 0, result.output

    class FailingPatient(Patient):
        def analyze_reads(self, *args):
            raise RuntimeError('failed')

    monkeypatch.setattr(cmd_hydra, 'PatientAnalyzer', FailingPatient)
    result = CliRunner().invoke(
        cmd_hydra.cli,
        [str(forward), '--output_dir', str(tmp_path / 'out2'),
         '--mutation_db', str(mutation_db), '--mean'],
    )
    assert result.exit_code == 2

    monkeypatch.setattr(
        cmd_quality,
        'QualityControl',
        lambda: SimpleNamespace(filter_reads=lambda *args: None),
    )
    result = CliRunner().invoke(
        cmd_quality.cli,
        [str(forward), '--output_dir', str(tmp_path / 'quality'), '--mean'],
    )
    assert result.exit_code == 0
