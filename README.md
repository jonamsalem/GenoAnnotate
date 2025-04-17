# GenoAnnotate

A Python tool for identifying pathogenic and novel genetic variants from per-sample VCF files, with automated report generation.
## Overview

GenoAnnotate helps researchers and clinicians analyze genetic variants by:
- Extracting variants from specific chromosomes and regions
- Filtering variants based on customizable criteria
- Annotating variants with pathogenicity scores and functional predictions
- Identifying novel variants not present in population databases
- Generates text files with gene annotations, clinical significance from ClinVar (CLNSIG), and REVEL scores (pathogenicity predictions for missense variants).
- Generate **aggregate** CSV and HTML reports of all annotated variant text files.
- 🚀 Coming Soon: No-Code Cloud Version for browser-based analysis

## Requirements

- Python 3.9+
- cyvcf2 (`pip install cyvcf2`)
- Properly formatted VCF files (bgzipped and indexed with tabix)!
- [`Annovar`](#annovar-instructions) (required for local annotation). See instructions below.


## Installation

```bash
# Clone the repository
git clone https://github.com/jonamsalem/GenoAnnotate.git
cd GenoAnnotate

# Install dependencies and activate the environment (conda suggested)
conda env create -f env.yaml
```

## Preparing VCF Files

Before using GenoAnnotate, ensure your VCF files are properly compressed and indexed:

```bash
# Compress with bgzip (if not already compressed)
bgzip -c your_file.vcf > your_file.vcf.gz

# Index with tabix
tabix -p vcf your_file.vcf.gz
```


### Example Usage

```bash
$ python3 variant_extract.py --path './vcfs' --chrom 'chrX' --annovar '~/annovar' --vaf '0.30' --start '101397803' --end '101407925'
```

Command-line Flags

| Flag         | Required | Description                                                                 |
|--------------|----------|-----------------------------------------------------------------------------|
| `--path`     | Yes      | Path to the directory containing `.vcf.gz` files. VCF files can be nested in sub-directories.                          |
| `--chrom`    | Yes      | Chromosome to extract from (e.g., `chrX`, `chr1`).                 |
| `--annovar`  | Yes      | Full path to the ANNOVAR directory containing `table_annovar.pl` and `humandb/`. |
| `--vaf`      | No       | Variant Allele Frequency threshold (e.g., `0.05`). Defaults to `0.05`.     |
| `--start`    | No       | Start genomic coordinate (e.g., `101397803`).|
| `--end`      | No       | End genomic coordinate (e.g., `101407925`).|
| `--ref`      | No       | Human reference genome version (hg38/hg19). Default is hg38. |
| `--report`   | No       |Generate a CSV and HTML summary report. Default is false. |




## Annovar Instructions

To run the program locally, [Annovar](https://annovar.openbioinformatics.org/en/latest/) must be installed.

Note that HPC environments may already have Annovar installed. 

**Option 1 (Preferred)**

Register on the Annovar website to download the most up-to-date version of Annovar.

**Option 2**

Download Annovar directly:

```bash
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar -xzvf annovar.latest.tar.gz
```

After installation, navigate to the `annovar` directory and execute the following commands:

For hg38:

```bash
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240917 humandb/ && 
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar revel humandb/ && 
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
```

For hg19:
```bash
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20240917 humandb/ && 
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar revel humandb/ && 
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
```

```bash
echo "export ANNOVAR_PATH_GENOANOTATE=\"$(pwd)\"" >> ~/.bashrc
source ~/.bashrc
echo $ANNOVAR_PATH_GENOANOTATE #copy the output as it will serve as the path for Annovar for GenoAnnotate
```

This script will save the path to Annovar as well as download necessary databases for annotations 

Clinvar provides clinically validated information about known variants

Revel helps predict the impact of novel or rare variants 

Refgene has gene definitions and information about gene structures in the human genome


## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
