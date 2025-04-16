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


### Basic Usage

```bash
$ python3 variant_extract.py --path 'vcfs' --chrom 'chrX' 

```

### Advanced Filtering

```python
from genoAnnotate import filter_variants

# Filter variants based on quality and depth
filter_variants('input.vcf.gz', 'output.vcf', min_quality=30, min_depth=20)
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

