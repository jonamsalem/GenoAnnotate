import glob
import os
from cyvcf2 import VCF, Writer
import argparse
import logging
from multiprocessing import Pool, cpu_count
import subprocess
from report import generate_report

logging.basicConfig(
    level=logging.INFO)

DEFAULT_VAF = 0.20 
DP = 20 

#read each vcf file and save fields to csv file using cyvcf2
def vcf_extract(file,temp_dir, chrom, start,end,provided_vaf):
    
    basename = os.path.basename(file)
    filename = basename.split('.vcf')[0] #extract basename of file 
    
    #check if file is gzipped and indexed using tabix
    if not file.endswith('.vcf.gz'):
        logging.error(f"File {file} is not a gzipped VCF file.")
        return None, False  
    if not os.path.exists(file + '.tbi'):
        logging.error(f"File {file} is not indexed. Please index the file before proceeding.")
        return None, False  

    if provided_vaf == "":
        provided_vaf = DEFAULT_VAF
    else:
        try:
            provided_vaf = float(provided_vaf)
        except ValueError:
            logging.error("Invalid VAF value. Please provide a numeric value.")
            return


    #check if file exists and if so add -1 to filename
    output_file_name = f'{temp_dir}/{filename}'

    if os.path.exists(output_file_name + '.vcf'):
        logging.warning(f"File {filename}.vcf already exists. Adding '-1' to filename.")
        output_file_name += '-1'
    

    vcf = VCF(file)
    region = f"{chrom}:{start}-{end}"
    variants = list(vcf(region))  # Evaluate the iterator once
    w = Writer(output_file_name, vcf)

    wrote_to_file = False
    variants_passed = 0
    # Filter variants
    for variant in variants:
        pass_filters = True 

        try:
            dp_values = variant.format('DP')
            if dp_values[0][0] < DP:
                pass_filters = False

            #calculate VAF from AD
            ref_freq = variant.format('AD')[0][0]
            alt_freq = variant.format('AD')[0][1]
           
            if ref_freq + alt_freq > 0:
                vaf = alt_freq / (ref_freq + alt_freq)
            else:
                vaf = 0

            if vaf < provided_vaf: 
                pass_filters = False
           
            if variant.QUAL < 30:  
                pass_filters = False   
             
        except Exception as e:
            logging.error(f"Error processing variant {variant}: {e}")
            pass_filters = False

        if pass_filters:
            w.write_record(variant)
            wrote_to_file = True
            variants_passed += 1
        
    vcf.close()
    w.close()

    # Check if any variants were written to the file
    if not wrote_to_file:
        logging.warning(f"No variants passed the filters for {file}")
        os.remove(output_file_name)

    logging.info(f"Extracted {variants_passed} variants from {file}")

    return output_file_name, wrote_to_file



#annotate vcf files using annovar
def annovar_annotate(annovar_path, outputs_dir, file, ref):
    basename = file.split('/')[-1].split('.vcf')[0]  # Extract the base name of the file
    command_annovar = f"{annovar_path}/table_annovar.pl {file} {annovar_path}/humandb/ -buildver {ref} -protocol refGene,clinvar_20240917,revel -operation g,f,f -nastring . -vcfinput -out {outputs_dir}/annotated_{basename}"
    
    try:
        subprocess.run(
            command_annovar,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"Annovar failed for {file} (exit {e.returncode}).")



def annotate(file, temp_dir, outputs_dir, chrom, start, end, vaf, annovar_path, ref):
    if ref != 'hg38' and ref != 'hg19':
        logging.warning(f"Invalid reference genome version {ref}. Defaulting to hg38.")
        ref = 'hg38'

    output_file, wrote_to_file = vcf_extract(file, temp_dir,chrom, start,end, vaf)
    if wrote_to_file:
        annovar_annotate(annovar_path, outputs_dir, output_file,ref)

        #extract the relevant columns from the output text file
        basename = os.path.basename(output_file)
        filename = basename.split('.vcf')[0] 

        output_file_name = f"{outputs_dir}/annotated_{filename}.{ref}_multianno"

        #pull out the relevant columns into csv file and delete the temp files
        #todo: find better way to do this
        command_txt = (
            f"awk -F'\\t' 'BEGIN {{OFS=\"\\t\"}} {{print $1, $2, $3, $4, $5, $6, $7, $9, $15, $24}}' \"{output_file_name}.txt\" > \"{output_file_name}.filtered.{chrom}.txt\""
        )
        
        subprocess.run(command_txt, shell=True)
        logging.info(f"Generated Annotated File: {output_file_name}.txt")

        command_cleanup = f"find . -name 'annotated_{basename}*' -not -name 'annotated_{basename}*.filtered*' -delete"
        subprocess.run(command_cleanup, shell=True)



def annotate_wrapper(args):
    return annotate(*args)

#main function to extract vcf data into temp csv files
def annotate_variants(temp_dir, outputs_dir, path, chrom,start="",end="", vaf="", annovar_path="", ref="hg38"):
    if annovar_path == "":
        logging.error("Invalid Annovar path. Please provide a valid path.")
        return
    
    command_find = 'find ' + path + ' -name "*.vcf.gz"' 
    all_files = subprocess.check_output(command_find, shell=True).decode('utf-8').splitlines()

    
    args_list = [(file, temp_dir,outputs_dir, chrom, start, end, vaf, annovar_path, ref) for file in all_files]

    #multiprocess the files
    with Pool(processes=cpu_count()) as pool:
        pool.map(annotate_wrapper, args_list)
    
    logging.info(f"{len(all_files)} VCF files processed.")

                



#function call and params
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract variant data from VCF files to CSV files.")
    parser.add_argument("--path", help="Path to the directory containing VCF files.")
    parser.add_argument("--chrom", required=True, help="Chromosome to extract (e.g., chrX).")
    parser.add_argument("--start", default="", help="Start position of the region to extract.")
    parser.add_argument("--end", default="", help="End position of the region to extract.")
    parser.add_argument("--vaf", default="", help="VAF threshold for filtering variants.")
    parser.add_argument("--annovar", default="", help="Annovar path to the directory containing annovar files.")
    parser.add_argument("--ref", default="hg38", help="Human reference genome version (default: hg38).")
    parser.add_argument("--report", default="False", help="Generate aggregate report (default: False).")

    args = parser.parse_args()
    temp_dir = 'temp_vcfs' #temporary directory to store intermediate files
    os.makedirs(temp_dir, exist_ok=True)

    outputs_dir = 'outputs'
    os.makedirs(outputs_dir, exist_ok=True)

    annotate_variants(temp_dir=temp_dir, outputs_dir=outputs_dir, path=args.path, chrom=args.chrom, start=args.start, end=args.end, vaf=args.vaf, annovar_path=args.annovar, ref=args.ref)

    #delete temp directory
    command_vcf_cleanp = f"rm -rf {temp_dir}"
    subprocess.run(command_vcf_cleanp, shell=True)

    #generate report
    if args.report.lower() == "true":
        generate_report(outputs_dir)
        logging.info("Generated HTML report.")

