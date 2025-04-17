import pandas as pd
from pathlib import Path
from pathlib import Path

def generate_report(output_dir):
    files = Path(output_dir).glob("*.txt")
    df = pd.DataFrame()
    
    #aggregate all files into a single DataFrame
    for file in files:
        temp_df = pd.read_csv(file, low_memory=False, sep='\t')
        df = pd.concat([df, temp_df], ignore_index=True)
    
    data = df.copy()
    
    if 'REVEL' in data.columns:
        data['REVEL'] = pd.to_numeric(data['REVEL'], errors='coerce')
    
    transitions = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
    transversions = [('A', 'C'), ('C', 'A'), ('G', 'T'), ('T', 'G'),
                     ('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')]
    
    data['mutation'] = list(zip(data['Ref'], data['Alt']))

    # Only classify SNVs (single nucleotide variants)
    snv_mask = (data['Ref'].str.len() == 1) & (data['Alt'].str.len() == 1)
    data['is_transition'] = False
    data['is_transversion'] = False
    data.loc[snv_mask, 'is_transition'] = data.loc[snv_mask, 'mutation'].isin(transitions)
    data.loc[snv_mask, 'is_transversion'] = data.loc[snv_mask, 'mutation'].isin(transversions)
        
    report = {}
    
    # 1. Total mutations
    report['Total Variants'] = len(data)
    
    # 2. Mutation types
    report['Exonic Variants'] = sum(data['Func.refGene'] == 'exonic')
    report['Intronic Variants'] = sum(data['Func.refGene'] == 'intronic')
    report['UTR Variants'] = sum(data['Func.refGene'].str.contains('UTR', na=False))
    report['Splicing Variants'] = sum(data['Func.refGene'] == 'splicing')
    
    # 3. Clinical significance
    if 'CLNSIG' in data.columns and not data['CLNSIG'].isna().all():
        report['Benign Variants'] = sum(data['CLNSIG'].str.contains('Benign', na=False) &
                                       ~data['CLNSIG'].str.contains('Likely_benign', na=False))
        report['Likely Benign Variants'] = sum(data['CLNSIG'].str.contains('Likely_benign', na=False))
        report['Pathogenic Variants'] = sum(data['CLNSIG'].str.contains('Pathogenic', na=False) &
                                           ~data['CLNSIG'].str.contains('Likely_pathogenic', na=False))
        report['Likely Pathogenic Variants'] = sum(data['CLNSIG'].str.contains('Likely_pathogenic', na=False))
        report['VUS Variants'] = sum(data['CLNSIG'].str.contains('uncertain', na=False, case=False))
        report['Conflicting Variants'] = sum(data['CLNSIG'].str.contains('Conflicting', na=False))
    
    # 4. Mutation characteristics
    report['Transitions'] = sum(data['is_transition'])
    report['Transversions'] = sum(data['is_transversion'])
    if report['Transversions'] > 0:
        report['Ts/Tv ratio'] = report['Transitions'] / report['Transversions']
    else:
        report['Ts/Tv ratio'] = float('nan')
    
    # 5. Genes with mutations
    report['Unique Genes'] = data['Gene.refGene'].nunique()

    if 'ExonicFunc.refGene' in data.columns:
        exonic_data = data[~data['ExonicFunc.refGene'].isna()]
        
        # Get counts of each exonic function type
        exon_func_counts = exonic_data['ExonicFunc.refGene'].value_counts()
        
        # Add each type to the report
        for func_type, count in exon_func_counts.items():
            report[f'Exonic {func_type}'] = count
            
        # Add specific counts for common types
        report['Nonsynonymous SNVs'] = sum(data['ExonicFunc.refGene'] == 'nonsynonymous SNV')
        report['Synonymous SNVs'] = sum(data['ExonicFunc.refGene'] == 'synonymous SNV')
        report['Frameshift Insertions'] = sum(data['ExonicFunc.refGene'].str.contains('frameshift insertion', na=False))
        report['Frameshift Deletions'] = sum(data['ExonicFunc.refGene'].str.contains('frameshift deletion', na=False))
        report['Stopgain Variants'] = sum(data['ExonicFunc.refGene'].str.contains('stopgain', na=False))
        report['Stoploss Variants'] = sum(data['ExonicFunc.refGene'].str.contains('stoploss', na=False))
    
    # 6. Variants with REVEL scores
    if 'REVEL' in data.columns:
        report['Variants with REVEL Scores'] = sum(~data['REVEL'].isna())
        high_revel = sum(data['REVEL'] > 0.5)
        report['Variants with High REVEL (>0.5)'] = high_revel
    
    # 7. Most mutated genes (top 5)
    top_genes = data['Gene.refGene'].value_counts().head(5)
    for i, (gene, count) in enumerate(top_genes.items()):
        report[f'Top Annotated Gene {i+1}'] = f"{gene} ({count})"
    

    report_df = pd.DataFrame(list(report.items()), columns=['Metric', 'Value'])
    report_df.to_csv(f'{output_dir}/summary_report.csv', index=False)
    report_df.to_html(f'{output_dir}/summary_report.html')