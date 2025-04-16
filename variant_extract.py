import glob
import os
from cyvcf2 import VCF
import pandas as pd
import argparse
import shutil
import requests
#read each vcf file and save fields to csv file using cyvcf2
def vcf_extract_csv(file,temp_dir, chrom, start,end):
    basename = os.path.basename(file)
    filename = basename.split('.vcf')[0] #extract basename of file

    vcf = VCF(file)
    region = f"{chrom}:{start}-{end}"
    variants = list(vcf(region))  # Evaluate the iterator once
    
    try:
        df = pd.DataFrame({
            'Chrom': [variant.CHROM for variant in variants],
            'Reference': [variant.REF for variant in variants],
            'Alternate': [variant.ALT[0] if variant.ALT else None for variant in variants],
            'Position': [variant.POS for variant in variants],
            'Depth': [variant.format('DP')[0] if variant.format('DP') else None for variant in variants],
            'Quality': [variant.QUAL for variant in variants],
            'Allelic_depth': [variant.format('AD')[0].tolist() if variant.format('AD') is not None and len(variant.format('AD')[0]) > 0 else None for variant in variants]})
    except ValueError as e:
        print(f"ValueError during DataFrame creation: {e}")
        return pd.DataFrame()
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return pd.DataFrame()

    vcf.close() #close file

    #check if file exists and if so add 1
    if os.path.exists(f'{temp_dir}/{filename}.csv'):
        filename += '-1'

    df.to_csv(f"{temp_dir}/{filename}.csv")


#main function to extract vcf data into temp csv files
#todo: implement multiprocessing or convert to bash for parallel jobs
def extract_variants(temp_dir, path, chrom,start="",end=""):
    path = '{}/*.vcf.gz'.format(path)
    all_files = glob.glob(path)
    os.makedirs(temp_dir, exist_ok=True)
    for file in all_files:
       vcf_extract_csv(file, temp_dir,chrom, start,end)

def extract_annotations(chrom, start, end):
    #use myvariant.info API to get annotations which has an aggregate of different databases
    base_url = "https://myvariant.info/v1/query"
    query_string = f"{chrom}:{start}-{end}"
    params = {
        'q': query_string,
        'assembly': 'hg38',  # Explicitly specify GRCh38
    }

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Raise an exception for bad status codes
        annotations = response.json()

    except requests.exceptions.RequestException as e:
        print(f"Error fetching annotations: {e}")
        return
    except Exception as e:
        print(f"Error parsing JSON response: {e}")
        return
    print(annotations)


#function call and params
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract variant data from VCF files to CSV files.")
    parser.add_argument("--path", help="Path to the directory containing VCF files.")
    parser.add_argument("--chrom", required=True, help="Chromosome to extract (e.g., chrX).")
    parser.add_argument("--start", default="", help="Start position of the region to extract.")
    parser.add_argument("--end", default="", help="End position of the region to extract.")

    args = parser.parse_args()
    temp_dir = 'temp_csvs'

    extract_variants(temp_dir,args.path, args.chrom, args.start, args.end)

    extract_annotations(args.chrom, args.start, '101397903')
