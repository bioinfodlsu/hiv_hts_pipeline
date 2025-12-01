import argparse
import sys

def parse_info_field(info_str):
    """
    Parse the INFO column of a VCF line into a dict of key: value
    """
    info = {}
    for entry in info_str.split(';'):
        if '=' in entry:
            key, value = entry.split('=', 1)
            info[key] = value
        else:
            info[entry] = True
    return info


def filter_vcf(input_vcf, output_vcf, min_af, max_af, min_dp, max_dp):
    with (open(input_vcf, 'r')) as fin, \
         (open(output_vcf, 'w')) as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue

            cols = line.rstrip('\n').split('\t')
            info_dict = parse_info_field(cols[7])  # INFO column

            # Fetch allele frequency (AF) and depth (DP) from INFO
            af = info_dict.get('AF')
            dp = info_dict.get('DP')

            try:
                af_val = float(af) if af is not None else None
            except ValueError:
                af_val = None

            try:
                dp_val = int(dp) if dp is not None else None
            except ValueError:
                dp_val = None

            # Apply filters
            if min_af is not None and (af_val is None or af_val < min_af):
                continue
            if max_af is not None and (af_val is None or af_val > max_af):
                continue
            if min_dp is not None and (dp_val is None or dp_val < min_dp):
                continue
            if max_dp is not None and (dp_val is None or dp_val > max_dp):
                continue

            fout.write(line)


def main():
    parser = argparse.ArgumentParser(
        description='Filter VCF by allele frequency (AF) and depth (DP) from INFO column.')
    parser.add_argument('-i', '--input', required=True,
                        help='Input VCF file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output filtered VCF file')
    parser.add_argument('--min-af', type=float, default=None,
                        help='Minimum allele frequency to keep')
    parser.add_argument('--max-af', type=float, default=None,
                        help='Maximum allele frequency to keep')
    parser.add_argument('--min-dp', type=int, default=None,
                        help='Minimum depth to keep')
    parser.add_argument('--max-dp', type=int, default=None,
                        help='Maximum depth to keep')

    args = parser.parse_args()

    filter_vcf(
        input_vcf=args.input,
        output_vcf=args.output,
        min_af=args.min_af,
        max_af=args.max_af,
        min_dp=args.min_dp,
        max_dp=args.max_dp,
    )

if __name__ == '__main__':
    main()
