import os
import pandas as pd
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import math
from collections import Counter


#extract Cambray data of interest and output to extracted_data.csv
columns = ["gs.sequence","gs.cdsCAI","gs.utrCdsStructureMFE","gs.fivepCdsStructureMFE","gs.threepCdsStructureMFE","gs.cdsNucleotideContentAT","clean.lin.prot.mean","clean.lin.prot.acf.mean","ss.rna.dna.mean","protect.rna.dna.mean","halflife.rna.dna.mean","GC","GC3"]
df = pd.read_csv("Cambray_DataS15.csv", usecols=columns)
df.to_csv('extracted_data.csv')

dfc = pd.read_csv('extracted_data.csv')
dfg = pd.read_csv('1241934tables1.csv')

#function takes dictionary with gene sequences as keys, calculates RNAfold estimated MFE. Takes sequences from Goodman.
def calculate_G5mfe(sequence):
    fiveprime = sequence[3:]
    process = subprocess.Popen(['RNAfold', '--noPS'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=fiveprime)
    mfe = float(stdout.splitlines()[1].split()[-1].strip('()'))
    return mfe


#function takes dictionary with gene sequences as keys, calculates RNAfold estimated MFE. Takes sequences from Cambray
def calculate_C5mfe(sequence):
    fiveprime = sequence[:33]
    process = subprocess.Popen(['RNAfold', '--noPS'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=fiveprime)
    mfe = float(stdout.splitlines()[1].split()[-1].strip('()'))
    return mfe

# function uses above function to return dataframe row and calculated stability
def G5_process_row(row):
    mfe = calculate_G5mfe(row['CDS.seq'])
    return row['CDS.seq'], mfe

def C5_process_row(row):
    mfe = calculate_C5mfe(row['gs.sequence'])
    return row['gs.sequence'], mfe

#function takes dataframe using parallel processing on rows to use above functions to get stability measure and create a dictionary for pandas dataframe with columns containing sequence and stability.
def process_G5dataframe(df):
    results = []
    total_rows = len(df)
    with ThreadPoolExecutor() as executor:
        future_to_row = {executor.submit(G5_process_row, row): row for _, row in df.iterrows()}
        for index, future in enumerate(as_completed(future_to_row)):
            try:
                sequence, mfe = future.result()
                results.append({'CDS.seq': sequence, 'STR_5_G': mfe})
                progress = (index + 1) / total_rows * 100
                print(f"Progress: {progress:.2f}% ({index + 1}/{total_rows} rows processed)", end='\r')
            except Exception as exc:
                print(f'An error occurred: {exc}')
    print(f"Processing complete. Total rows processed: {total_rows}")
    result_df = pd.DataFrame(results)
    return result_df

def process_C5dataframe(df):
    results = []
    total_rows = len(df)
    with ThreadPoolExecutor() as executor:
        future_to_row = {executor.submit(C5_process_row, row): row for _, row in df.iterrows()}
        for index, future in enumerate(as_completed(future_to_row)):
            try:
                sequence, mfe = future.result()
                results.append({'gs.sequence': sequence, 'STR_5_G': mfe})
                progress = (index + 1) / total_rows * 100
                print(f"Progress: {progress:.2f}% ({index + 1}/{total_rows} rows processed)", end='\r')
            except Exception as exc:
                print(f'An error occurred: {exc}')
    print(f"Processing complete. Total rows processed: {total_rows}")
    result_df = pd.DataFrame(results)
    return result_df



G5_df = process_G5dataframe(dfg)
C5_df = process_C5dataframe(dfc)

#output to csv
G5_df.to_csv('GSTR.csv', index=False)
C5_df.to_csv('CSTR.csv', index=False)





