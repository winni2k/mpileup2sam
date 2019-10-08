import re
import subprocess
from pathlib import Path

import pytest


class MpileupExpectation:

    def __init__(self, expected, actual):
        self.expected = expected
        self.actual = actual

    def are_equal(self):
        with open(self.expected, 'r') as fh:
            expected = strip_start_end_chars(fh.read())
        with open(self.actual, 'r') as fh:
            actual = strip_start_end_chars(fh.read())

        assert actual == expected


def strip_start_end_chars(text):
    return re.sub(r'(\^.|\$)', '', text)


def run_conversion(input_mpileup, tmpdir):
    reference = 'test/cases/simulated.fa'
    out_prefix = tmpdir / "split_%#.%."
    bam_list = tmpdir / 'bam.list'
    out_sam = tmpdir / 'out.sam'
    out_mpileup = tmpdir / 'out.mpileup'

    subprocess.run(
        ['mpileup2sam', '--reference', reference, str(input_mpileup), str(out_sam)],
        check=True)
    subprocess.run(['samtools', 'split', '-f', str(out_prefix), str(out_sam)])
    n_samples = len(list(Path(tmpdir).glob('split_*.bam')))
    with open(bam_list, 'w') as fh:
        for samp_idx in range(n_samples):
            fh.write(f'{tmpdir}/split_{samp_idx}.bam\n')

    subprocess.run(
        ['samtools', 'mpileup', '-A', '-B',
         '-d', '100000',
         '-f', str(reference),
         '-b', str(bam_list),
         '-o', str(out_mpileup)],
        check=True)
    return MpileupExpectation(input_mpileup, out_mpileup)


@pytest.mark.parametrize('head_len', [1, 100])
def test_headX(tmpdir, head_len):
    # given
    orig_mpileup = 'test/cases/simulated.head100.mpileup'
    input_mpileup = tmpdir / 'simulated.mpileup'

    with open(orig_mpileup) as ifh:
        with open(input_mpileup, 'w') as ofh:
            for idx, line in enumerate(ifh):
                if idx == head_len:
                    break
                ofh.write(line)

    # when
    expect = run_conversion(input_mpileup, tmpdir)

    # then
    expect.are_equal()
