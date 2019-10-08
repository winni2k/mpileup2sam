import difflib
import io
import re

from mpileup2sam import Converter


def make_spaces_tabs(input):
    return re.sub(r' +', '\t', input)


def make_displayable(input):
    return re.sub(r'\t', ' ', input)


def test_convert_two_samples_one_read():
    # given
    input = make_spaces_tabs(
        '1       10000   N       0       *       *       0       *       *\n'
        '1       10005   C       0       *       *       1       ^I.     F'
    )
    in_fh = io.StringIO(input)
    converter = Converter(fh=in_fh, sample_ids=['samp_0', 'samp_1'], ref_seqs=['1'],
                          ref_seq_sizes=[10005])

    # when
    lines = list(converter.lines())

    # then
    assert ''.join(lines) == make_spaces_tabs('@HD VN:1.6 SO:unknown\n'
                                              '@SQ SN:1 LN:10005\n'
                                              '@RG ID:0 SM:samp_0\n'
                                              '@RG ID:1 SM:samp_1\n'
                                              'r0 0 1 10005 0 1M * 0 0 C F RG:Z:1\n'), difflib.Differ()


def test_convert_one_sample_one_read():
    # given
    input = make_spaces_tabs(
        '1       10005   C       1       ^I.     F\n'
    )
    in_fh = io.StringIO(input)
    converter = Converter(fh=in_fh, sample_ids=['samp_0'], ref_seqs=['1'],
                          ref_seq_sizes=[10005])

    # when
    lines = list(converter.lines())

    # then
    assert ''.join(lines) == make_spaces_tabs('@HD VN:1.6 SO:unknown\n'
                                              '@SQ SN:1 LN:10005\n'
                                              '@RG ID:0 SM:samp_0\n'
                                              'r0 0 1 10005 0 1M * 0 0 C F RG:Z:0\n'), difflib.Differ()


def test_convert_one_sample_two_reads():
    # given
    input = make_spaces_tabs(
        '1       10005   C       2       ^I..$     FG\n'
    )
    in_fh = io.StringIO(input)
    converter = Converter(fh=in_fh, sample_ids=['samp_0'], ref_seqs=['1'],
                          ref_seq_sizes=[10005])

    # when
    lines = list(converter.lines())

    # then
    assert ''.join(lines) == make_spaces_tabs(
        '@HD VN:1.6 SO:unknown\n'
        '@SQ SN:1 LN:10005\n'
        '@RG ID:0 SM:samp_0\n'
        'r0 0 1 10005 0 1M * 0 0 C F RG:Z:0\n'
        'r1 0 1 10005 0 1M * 0 0 C G RG:Z:0\n'
    ), difflib.Differ()


def test_two_samples_five_reads_two_sites():
    # given
    input = make_spaces_tabs(
        '1       10003   A       3       ^I,t.$     FG5      0       *         *\n'
        '1       10005   C       0       *          *        2       ^I.A$     LM\n'
    )
    in_fh = io.StringIO(input)
    converter = Converter(fh=in_fh, sample_ids=['samp_0', 'samp_1'], ref_seqs=['1'],
                          ref_seq_sizes=[10005])

    # when
    lines = list(converter.lines())

    # then
    assert ''.join(lines) == make_spaces_tabs(
        '@HD VN:1.6 SO:unknown\n'
        '@SQ SN:1 LN:10005\n'
        '@RG ID:0 SM:samp_0\n'
        '@RG ID:1 SM:samp_1\n'
        'r0 16 1 10003 0 1M * 0 0 A F RG:Z:0\n'
        'r1 16 1 10003 0 1M * 0 0 T G RG:Z:0\n'
        'r2 0 1 10003 0 1M * 0 0 A 5 RG:Z:0\n'
        'r3 0 1 10005 0 1M * 0 0 C L RG:Z:1\n'
        'r4 0 1 10005 0 1M * 0 0 A M RG:Z:1\n'
    ), difflib.Differ()
