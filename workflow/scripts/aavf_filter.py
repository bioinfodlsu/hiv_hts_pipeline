import argparse

def parse_info_field(info_str):
    """
    Parse the INFO column of an AAVF line into a dict of key: value
    """
    info = {}
    for entry in info_str.split(';'):
        if '=' in entry:
            key, value = entry.split('=', 1)
            info[key] = value
        else:
            info[entry] = True
    return info


def filter_aavf(input_file, output_file,
                min_acc, max_acc, min_acf, max_acf):
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            # Write meta and header lines unchanged
            if line.startswith('#'):
                fout.write(line)
                continue

            cols = line.rstrip('\n').split('\t')
            info_dict = parse_info_field(cols[8])  # INFO column at index 8

            # Extract ACC and ACF
            acc = info_dict.get('ACC')
            acf = info_dict.get('ACF')

            try:
                acc_val = int(acc) if acc is not None else None
            except ValueError:
                acc_val = None
            try:
                acf_val = float(acf) if acf is not None else None
            except ValueError:
                acf_val = None

            # Apply filters
            if min_acc is not None and (acc_val is None or acc_val < min_acc):
                continue
            if max_acc is not None and (acc_val is None or acc_val > max_acc):
                continue
            if min_acf is not None and (acf_val is None or acf_val < min_acf):
                continue
            if max_acf is not None and (acf_val is None or acf_val > max_acf):
                continue

            fout.write(line)


def main():
    parser = argparse.ArgumentParser(
        description='Filter AAVF by ACC and ACF values in INFO column.')
    parser.add_argument('-i', '--input', required=True,
                        help='Input AAVF file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output filtered AAVF file')
    parser.add_argument('--min-acc', type=int, default=None,
                        help='Minimum alternate codon coverage (ACC) to keep')
    parser.add_argument('--max-acc', type=int, default=None,
                        help='Maximum alternate codon coverage (ACC) to keep')
    parser.add_argument('--min-acf', type=float, default=None,
                        help='Minimum alternate codon frequency (ACF) to keep')
    parser.add_argument('--max-acf', type=float, default=None,
                        help='Maximum alternate codon frequency (ACF) to keep')

    args = parser.parse_args()

    filter_aavf(
        input_file=args.input,
        output_file=args.output,
        min_acc=args.min_acc,
        max_acc=args.max_acc,
        min_acf=args.min_acf,
        max_acf=args.max_acf
    )

if __name__ == '__main__':
    main()
