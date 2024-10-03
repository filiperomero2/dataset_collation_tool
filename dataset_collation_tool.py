#!/usr/env python

# Written by Filipe Moreira (PhD)
# Federal University of Rio de Janeiro, Brazil
# 2024-10-03

import os
import sys
import subprocess
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO

def get_args():

    parser = argparse.ArgumentParser(
    description='A simple script to collate datasets for virus genomic epidemiology.',
    usage='''dataset_collator.py [args]''')

    parser.add_argument('--taxid',type=int,
    help='NCBI tax id for target virus',
    required = True)

    parser.add_argument('--target',type = int,
    help='Number of sequences to include per location and period (Default = 1)',
    nargs='?',const=1, default=1)

    parser.add_argument('--length-filter',type = int,
    help='Minimum sequence length to be included in the dataset (Default = 10000)',
    nargs='?',const=1, default=10000)
    
    parser.add_argument('--ambiguity-filter',type = float,
    help='Maximum percentage [0 - 1] of ambiguity characters (Default = 0.25)',
    nargs='?',const=1, default=0.25)

    parser.add_argument('--geo-unit', type = str,
    help='Geographic level to sample (default: country).',
    choices = ['country','continent'],
    nargs='?', default = 'country')

    # if year-month option selected, it assumes year-only seqs are from Jan
    parser.add_argument('--time-unit', type = str,
    help='Temporal level to sample (default: year).',
    choices = ['year','year-month'],
    nargs='?', default = 'year')

    parser.add_argument('--output', 
    help='Complete path for output directory.',
    required = True)

    parser.add_argument('--tag', type = str,
    help='Tag (name) for subsetting scheme.',
    nargs='?', default = 'subset')

    parser.add_argument('--skip-download',
    help='Option to skip the download step (it assumes the output diretory already exists and contains the ncbi_dataset/).',
    default=False,
    action=argparse.BooleanOptionalAction)

    args = vars(parser.parse_args())
    return args

def validate_args(args):

    if args['skip_download']:
        print("--skip-download option identified. Sarching for the ncbi_dataset/ directory.")
        my_path = f"{args['output']}/ncbi_dataset"
        my_path_2 = f"{args['output']}/ncbi_dataset/{args['tag']}"
        if(os.path.isdir(my_path)):
            print('ncbi_dataset/ identified. Continuing...')
            if(os.path.isdir(my_path_2)):
                print(f'Tag ->{args['tag']}<- for naming subsetting output directory was already used. Please provide a new one.')
                exit()
            else:
                print('Valid tag will be employed.')
        else:
            print('ncbi_dataset/ not identified. Please check provided paths.')
            exit()
    else:
        if(os.path.isdir(args['output'])):
            print('Output file already exists (and --skip-download option was off). Please check paths.')
            exit()
        else:
            print('Output path does not exist. Continuing...')

    print("All arguments were succesfully verified.")

    return(args)

def download_ncbi_datasets(args):
    
    command_create_dir = f"mkdir -p {args['output']}".split()
    subprocess.run(command_create_dir)

    command_download_dataset = f"datasets download virus genome taxon {args['taxid']} --filename {args['output']}/ncbi_package.zip".split()
    subprocess.run(command_download_dataset)

    command_unzip = f"unzip {args['output']}/ncbi_package.zip -d {args['output']}".split()
    subprocess.run(command_unzip)

    outfile = open(f"{args['output']}/ncbi_dataset/data/data_report.tsv", "w")
    command_convert_metadata_to_csv = f"dataformat tsv virus-genome --inputfile {args['output']}/ncbi_dataset/data/data_report.jsonl".split()
    subprocess.run(command_convert_metadata_to_csv,stdout=outfile)

    outfile_2 = open(f"{args['output']}/ncbi_dataset/data/genomic.renamed.fasta", "w")
    command_rename_fasta = ['sed','s/ .*//g',f'{args['output']}/ncbi_dataset/data/genomic.fna']
    subprocess.run(command_rename_fasta,stdout=outfile_2)

    print("The NCBI data package was downloaded and formatted.")

    return

def filter_and_sample_sequences(args):

    geo_unit = args['geo_unit']
    time_unit = args['time_unit']

    command_create_dir = f"mkdir -p {args['output']}/ncbi_dataset/{args['tag']}".split()
    subprocess.run(command_create_dir)

    metadata_path = f"{args['output']}/ncbi_dataset/data/data_report.tsv"
    metadata = pd.read_csv(metadata_path,sep='\t', header=0)

    seq_path = f"{args['output']}/ncbi_dataset/data/genomic.renamed.fasta"
    sequences = SeqIO.to_dict(SeqIO.parse(seq_path, "fasta"))
    metadata = count_ambiguities(sequences,metadata)

    pass_length_filter = metadata['Length'] >= args['length_filter']
    metadata_filt = metadata[pass_length_filter].copy()

    pass_ambiguity_filter = metadata_filt['ambiguity'] <= args['ambiguity_filter']
    metadata_filt = metadata_filt[pass_ambiguity_filter].copy()

    country = metadata_filt['Geographic Location'].str.split(':').str[0]
    metadata_filt.loc[:,'country'] = country

    location = metadata_filt['Geographic Location'].str.split(':').str[1].str.replace(',','').str.replace(' ','_').str[1:]
    metadata_filt.loc[:,'location'] = location
    metadata_filt['location'] = metadata_filt['location'].fillna('NA')

    metadata_filt.rename(columns={'Geographic Region': 'continent'}, inplace=True)

    metadata_filt = metadata_filt.dropna(subset=[geo_unit,'Isolate Collection date']).copy()
    metadata_filt['Isolate Collection date filt'] = pd.to_datetime(metadata_filt['Isolate Collection date'],format='mixed')
    metadata_filt['year'] = metadata_filt['Isolate Collection date filt'].dt.year
    metadata_filt['month'] = metadata_filt['Isolate Collection date filt'].dt.month
    metadata_filt['year-month'] = metadata_filt['year'].astype(str) + "-" + metadata_filt['month'].astype(str)

    unique_time = np.sort(metadata_filt[time_unit].unique())

    target = args['target']
    sampled_indexes = pd.Index([])

    for time in unique_time:
        time_id = metadata_filt[time_unit] == time
        geo = metadata_filt[geo_unit][time_id]
        unique_geo = geo.unique()
        for location in unique_geo:
            df = metadata_filt[(metadata_filt[time_unit] == time) & (metadata_filt[geo_unit] == location)]
            if df.shape[0] <= target:
                my_sample = df.index
            else:
                my_sample = df.sample(n=target).index
            sampled_indexes = sampled_indexes.append(my_sample)

    sub_df = metadata_filt.loc[sampled_indexes]

    print("The dataset was filtered and sampled.")

    write_outputs(sequences,sub_df,args)

    return

def count_ambiguities(sequences,metadata):

    ambiguity = {}
    for sequence_name in sequences.keys():
        sequence = sequences[sequence_name]
        sequence_length = len(sequence.seq)
        count_A = sequence.count("A")
        count_T = sequence.count("T")
        count_C = sequence.count("C")
        count_G = sequence.count("G")
        non_ambiguity_chars = count_A + count_T + count_C + count_G
        ambiguity_percentage = ((sequence_length - non_ambiguity_chars) / sequence_length)
        ambiguity[sequence_name] = ambiguity_percentage
    
    metadata_cp = metadata.copy()
    metadata_cp['ambiguity'] = metadata_cp['Accession'].map(ambiguity)
    return metadata_cp

def write_outputs(sequences, sub_df, args):

    accessions = sub_df['Accession']

    seq_objects_list = []
    for accession in accessions:

        seq = sequences[accession]
        
        my_country = sub_df['country'][sub_df['Accession'] == accession].astype('string')
        my_country = ''.join(my_country).replace(' ','_').replace('\'','_')

        my_location = sub_df['location'][sub_df['Accession'] == accession].astype('string')
        my_location = ''.join(my_location)

        my_date = sub_df['Isolate Collection date'][sub_df['Accession'] == accession].astype('string')
        my_date = ''.join(my_date)

        my_date = sub_df['Isolate Collection date'][sub_df['Accession'] == accession].astype('string')
        my_date = ''.join(my_date)

        new_seq_name = accession + "|" + my_country + "|" + my_location + "|" + my_date
        seq.id = new_seq_name
        seq.name = ""
        seq.description = ""

        seq_objects_list.append(seq)


    output_fasta = f"{args['output']}/ncbi_dataset/{args['tag']}/sampled_sequences.fasta"
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(seq_objects_list, output_handle, "fasta")

    output_metadata = f"{args['output']}/ncbi_dataset/{args['tag']}/sampled_metadata.tsv"
    sub_df.to_csv(output_metadata, sep = "\t", index=False)

    print('Sequence and filtered output files were written to files.')

    return

def main():
    args = get_args()
    validate_args(args)
    if not args['skip_download']:
        download_ncbi_datasets(args)
    filter_and_sample_sequences(args)
    print('The end.')
    return

if __name__  == '__main__':
    sys.exit(main())
