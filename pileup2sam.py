import itertools
import re

import click


class BufferedLineReader:

    def __init__(self, fh):
        self.fh = fh
        self.buffer = []

    def peek_line(self):
        self.buffer.append(self.fh.readline())
        return self.buffer[-1]

    def __iter__(self):
        return itertools.chain(self.buffer, self.fh.readlines())


def generate_sample_names_from_pileup_line(line):
    return [f'sample_{idx}' for idx in range((len(line.split()) // 3) - 1)]


def get_ref_seqs_and_sizes_from_ref_filehandle(ref_file):
    import pysam

    ref = pysam.FastaFile(filename=ref_file)
    return ref.references, ref.lengths


class Converter:

    def __init__(self, fh, sample_ids, ref_seqs, ref_seq_sizes):
        self.sample_ids = sample_ids
        self.fh = fh
        self.ref_seqs = ref_seqs
        self.ref_seq_sizes = ref_seq_sizes

    @classmethod
    def from_pileup_file_handle_and_reference(cls, fh, ref_file):
        fh = BufferedLineReader(fh)
        sample_ids = generate_sample_names_from_pileup_line(fh.peek_line())
        ref_seqs, ref_seq_sizes = get_ref_seqs_and_sizes_from_ref_filehandle(ref_file)
        return cls(fh=fh, sample_ids=sample_ids, ref_seqs=ref_seqs, ref_seq_sizes=ref_seq_sizes)

    def lines(self):
        yield from self.header_lines()
        yield from self.body_lines()

    def header_lines(self):
        yield '@HD\tVN:1.6\tSO:unknown\n'
        for idx, ref_seq in enumerate(self.ref_seqs):
            yield f'@SQ\tSN:{ref_seq}\tLN:{self.ref_seq_sizes[idx]}\n'
        for idx, sample in enumerate(self.sample_ids):
            yield f'@RG\tID:{idx}\tSM:{sample}\n'

    def body_lines(self):
        read_num = 0
        read_start_ends = re.compile(r'(\^.|\$)')
        for line in self.fh:
            fields = line.split()
            n_samples = (len(fields) // 3) - 1
            assert n_samples == len(self.sample_ids), (n_samples, len(self.sample_ids))
            for samp_idx, sample_id in enumerate(self.sample_ids):
                sample_offset = 3 + samp_idx * 3
                if fields[sample_offset] == '0':
                    continue
                read_chars = re.sub(read_start_ends, '', fields[sample_offset + 1])
                assert len(read_chars) == int(fields[sample_offset])
                assert len(read_chars) == len(fields[sample_offset + 2])
                for read_char_idx, read_char in enumerate(read_chars):
                    bit_mask, seq = read_char_to_seq(read_char, fields[2])
                    yield (f'r{read_num}\t{bit_mask}\t{fields[0]}\t{fields[1]}\t0\t1M\t*\t0\t0\t'
                           f'{seq}\t{fields[sample_offset + 2][read_char_idx]}\tRG:Z:{samp_idx}\n')
                    read_num += 1


def read_char_to_seq(char, ref):
    if char == '.':
        return 0x0, ref
    if char == ',':
        return 0x10, ref
    if char.islower():
        return 0x10, char.upper()
    return 0x0, char


@click.command()
@click.argument('input', type=click.File('rt'))
@click.argument('output', type=click.File('wt'))
@click.option('-r', '--reference', required=True, type=click.Path(exists=True),
              help='Indexed FASTA reference file that was used to create the pileup input file')
def cli(input, output, reference):
    converter = Converter.from_pileup_file_handle_and_reference(fh=input, ref_file=reference)
    for line in converter.lines():
        output.write(line)
    return 0


if __name__ == '__main__':
    cli()
