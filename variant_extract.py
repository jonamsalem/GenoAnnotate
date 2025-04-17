import glob
import os
from cyvcf2 import VCF, Writer
import argparse
import logging



#read each vcf file and save fields to csv file using cyvcf2
def vcf_extract(file,temp_dir, chrom, start,end,provided_vaf):
    
    basename = os.path.basename(file)
    filename = basename.split('.vcf')[0] #extract basename of file

    if provided_vaf == "":
        provided_vaf = 0.05
    else:
        try:
            provided_vaf = float(provided_vaf)
        except ValueError:
            logging.error("Invalid VAF value. Please provide a numeric value.")
            return


    #check if file exists and if so add -1 to filename
    output_file = f'{temp_dir}/{filename}.vcf'

    if os.path.exists(output_file):
        logging.warning(f"File {filename}.vcf already exists. Adding '-1' to filename.")
        output_file += '-1'
    

    vcf = VCF(file)
    region = f"{chrom}:{start}-{end}"
    variants = list(vcf(region))  # Evaluate the iterator once
    w = Writer(output_file, vcf)

    pass_filters = True
    wrote_to_file = False

    # Filter variants based on quality and depth and write to VCF
    for variant in variants:
        pass_filters = True  # Reset for the next variant

        try:
            if 'DP' in variant.FORMAT:
                dp_values = variant.format('DP')
                if dp_values[0][0] < 20:  # Minimum depth of 20
                    pass_filters = False

            #calculate VAF from AD
            ref_freq = variant.format('AD')[0][0]
            alt_freq = variant.format('AD')[0][1]
           
            if ref_freq + alt_freq > 0:
                vaf = alt_freq / (ref_freq + alt_freq)
            else:
                vaf = 0
            if vaf < provided_vaf:  # Minimum VAF of 5%
                pass_filters = False
           
            # Filter based on quality
            if variant.QUAL < 30:  # Minimum quality of 20
                pass_filters = False   
             
        except Exception as e:
            logging.error(f"Error processing variant {variant}: {e}")
            pass_filters = False

        if pass_filters:
            w.write_record(variant)
            wrote_to_file = True
        
    vcf.close()
    w.close()

    # Check if any variants were written to the file
    if not wrote_to_file:
        logging.warning(f"No variants passed the filters for {file}")
        os.remove(output_file)


    return output_file, wrote_to_file, filename


#annotate vcf files using annovar
def annovar_annotate(annovar_path, file):
    basename = file.split('/')[-1].split('.vcf')[0]  # Extract the base name of the file
    command = f"{annovar_path}/table_annovar.pl  {file}  {annovar_path}/humandb/ -buildver hg38 -protocol refGene,clinvar_20240917,revel -operation g,f,f -nastring . -vcfinput -out annotated_{basename}"
    os.