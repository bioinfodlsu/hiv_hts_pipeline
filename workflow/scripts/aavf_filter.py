import argparse
import datetime
import pandas as pd

def copy_metadata(aavf_file, out_file):
    out_file.write("##fileformat=AAVFv1.0\n")
    out_file.write(f"##fileDate={datetime.datetime.now().strftime('%Y%m%d')}\n")
    out_file.write("##source=aavf_filter.py\n")

    for line in aavf_file:
        if line.startswith("##reference") or line.startswith("##INFO") or line.startswith("##FILTER"):
            out_file.write(line)
        else:
            break

def add_filter_metadata(aavf_file, filter_params, out_file):
    out_file.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
    filter_ids = [
        "AF<" + filter_params["max_AF"] if filter_params["max_AF"] else None,
        "AF>" + filter_params["min_AF"] if filter_params["min_AF"] else None,
        "DP<" + filter_params["max_cov"] if filter_params["max_cov"] else None,
        "DP>" + filter_params["min_cov"] if filter_params["min_cov"] else None
    ]

    if filter_params["max_AF"]:
        out_file.write(f"##FILTER=<ID={filter_ids[0]},Description=\"Allele frequency less than {filter_params['max_AF']}\">\n")

    if filter_params["min_AF"]:
        out_file.write(f"##FILTER=<ID={filter_ids[1]},Description=\"Allele frequency greater than {filter_params['min_AF']}\">\n")

    if filter_params["max_cov"]:
        out_file.write(f"##FILTER=<ID={filter_ids[2]},Description=\"Coverage less than {filter_params['max_cov']}\">\n")

    if filter_params["min_cov"]:
        out_file.write(f"##FILTER=<ID={filter_ids[3]},Description=\"Coverage greater than {filter_params['min_cov']}\">\n")

    return filter_ids

def parse_aavf(aavf_file):
    # store info into a dataframe
    # for line in aavf_file:
    #     if line.startswith("#"):
    #         continue
    aavf_headers = ["CHROM", "GENE", "POS", "REF", "ALT", "FILTER", "ALT_FREQ", "COVERAGE", "INFO"]
    df = pd.read_csv(aavf_file, sep="\t", comment="#", header=None, names=aavf_headers)

    info_df = df['INFO'].str.split(';', expand=True)
    print(info_df)
    split_data = info_df.apply(lambda x: x.str.split('=', expand=True))
    print(split_data)
    aavf_df = pd.concat([df, pd.DataFrame(info_df.tolist(), index=df.index)], axis=1)
    aavf_df = aavf_df.drop(columns=['INFO'])

    return aavf_df        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--aavf", required=True, help="AAVF input to filter")
    parser.add_argument("--out", required=True, help="Output file")
    parser.add_argument("-A", required=False, help="Maximum allele frequency")
    parser.add_argument("-a", required=False, help="Minimum allele frequency")
    parser.add_argument("-C", required=False, help="Maximum coverage")
    parser.add_argument("-c", required=False, help="Minimum coverage")
    args = parser.parse_args()

    filter_params = {
        "max_AF": None,
        "min_AF": None,
        "max_cov": None,
        "min_cov": None
    }

    if args.A:
        filter_params["max_AF"] = args.A

    if args.a:
        filter_params["min_AF"] = args.a
    
    if args.C:
        filter_params["max_cov"] = args.C
    
    if args.c:
        filter_params["min_cov"] = args.c

    

    with open(args.aavf) as aavf_file, open(args.out, 'w') as out_file:
        copy_metadata(aavf_file, out_file)
        filter_ids = add_filter_metadata(aavf_file, filter_params, out_file)
        
        out_file.write("#CHROM\tGENE\tPOS\tREF\tALT\tFILTER\tALT_FREQ\tCOVERAGE\tINFO\n")

        aavf_data = parse_aavf(aavf_file)
        print(aavf_data)

if __name__ == "__main__":
    main()