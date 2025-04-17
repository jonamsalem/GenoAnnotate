# GenoAnnotate

A Python tool for identifying pathogenic and novel genetic variants from VCF files.

## Overview

GenoAnnotate helps researchers and clinicians analyze genetic variants by:
- Extracting variants from specific chromosomes and regions
- Filtering variants based on customizable criteria
- Annotating variants with pathogenicity scores and functional predictions
- Identifying novel variants not present in population databases

## Requirements

- Python 3.6+
- cyvcf2 (`pip install cyvcf2`)
- Properly formatted VCF files (bgzipped and indexed with tabix)
- [`Annovar`](#annovar-instructions) (for local execution). See instructions below.

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/GenoAnnotate.git
cd GenoAnnotate

# Install dependencies
pip install -r requirements.txt 
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
$ python3 variant_extract.py --path './vcfs' --chrom 'chrX' --annovar '~/annovar' --vaf '0.05' --start '101397803' --end '101407925'
```

Command-line Flags

| Flag         | Required | Description                                                                 |
|--------------|----------|-----------------------------------------------------------------------------|
| `--path`     | Yes      | Path to the directory containing `.vcf.gz` files.                           |
| `--chrom`    | Yes      | Chromosome to extract from (e.g., `chrX`, `chr1`).                 |
| `--annovar`  | Yes      | Full path to the ANNOVAR directory containing `table_annovar.pl` and `humandb/`. |
| `--vaf`      | No       | Variant Allele Frequency threshold (e.g., `0.05`). Defaults to `0.05`.     |
| `--start`    | No       | Start genomic coordinate (e.g., `101397803`). Optional if `--end` is not used. |
| `--end`      | No       | End genomic coordinate (e.g., `101407925`). Optional if `--start` is not used. |




## Annovar Instructions


To run the program locally, [Annovar](https://annovar.openbioinformatics.org/en/latest/) must be installed.

**Option 1 (Preferred)**

Register on the Annovar website to download the most up-to-date version of Annovar.

**Option 2**

Download Annovar directly:

```bash
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar -xzvf annovar.latest.tar.gz
```

After installation, navigate to the `annovar` directory and execute the following commands:

```bash
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240917 humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar revel humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_exome humandb/
```

These databases are used to annotate variant files. 


## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

