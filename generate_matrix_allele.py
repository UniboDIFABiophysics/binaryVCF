import io
import os
import argparse
import numpy as np
import pandas as pd

zygosity_value = {
    '0/1': 1,
    '1/1': 2
}

def extract_patient_id_from_vcf_path(path: str):
    return os.path.basename(path).split('_')[0]

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('#')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        names=['CHROM',	'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'EXTRA'],
        dtype={'CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str, 'EXTRA': str },
        sep='\t',
    )

def load_dataframe(path):
    df = read_vcf(path)
    considered_vcf_cols = ['CHROM', 'POS', 'REF', 'ALT', 'EXTRA']
    df = df[considered_vcf_cols]
    df['patientId'] = extract_patient_id_from_vcf_path(path)
    df['variant'] = df['CHROM'] + '-' + df['POS'].astype(str) + '-' + df['REF'] + '-' + df['ALT']
    df['zygosity'] = df['EXTRA'].apply(lambda extra: extra.split(':')[0])
    return df[(df['zygosity'] == '0/1') | (df['zygosity'] == '1/1')].drop(columns = considered_vcf_cols)

def main(vcf_directory, output_file):
    df = pd.concat(load_dataframe(vcf_directory + path) for path in os.listdir(vcf_directory) if path.endswith(".vcf"))
    df['zygosity_value'] = df['zygosity'].apply(lambda zyg: zygosity_value[zyg])
    matrix = pd.crosstab(df['variant'], df['patientId'], values=df['zygosity_value'], aggfunc='first')
    matrix.fillna(0).astype('int8').to_csv(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_dir')
    parser.add_argument('--out')
    args = parser.parse_args()
    main(args.vcf_dir, args.out)
