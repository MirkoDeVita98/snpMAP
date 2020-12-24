def file_length(fname: str) -> int:
    return sum(1 for line in open(fname))


def no_info_lines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            if l.startswith('#CHROM'):
                return i + 1
    return -1