#!/usr/bin/python3.6

import os, sys
import numpy as np
import pandas as pd
import argparse as ap
import configparser
import ast
import re
from os import path


def main(path_in_list, path_out, cols_in, saving_mode):

    if os.path.isdir(path_out) != True:
        os.makedirs(path_out)


    list_folder_tsv = []

    with open(path_in_list, "r") as path_list:

        for line in path_list.readlines():
            list_folder_tsv.append(line[:-1])
        path_list.close()

    
    for i, name_file in enumerate(list_folder_tsv):


        END_TSV_NAME = "_consensusSV_annotsv.tsv" if re.search("_consensusSV_annotsv.tsv", name_file) else "_annotsv.tsv"

        cols_in.append(name_file.split("/")[-1][:-len(END_TSV_NAME)])


        df = pd.read_csv(name_file, sep = '\t')
        df = df.reindex(columns=cols_in)

        # abs della length
        df['SV_length'] = df['SV_length'].abs()

        # apply filters
        ind_f = df['FILTER'] == 'PASS'
        ind_a = df['Annotation_mode'] == 'full'

        ind_tot = np.all([ind_f, ind_a], axis=0).tolist()

        df_filt = df.loc[ind_tot, cols_in]


        df_filt = df.reindex(index=ind_tot)
        df = df[df['FILTER'] == 'PASS']
        df = df[df['Annotation_mode'] == 'full']
        df_filt = df
        cols_in = cols_in[:-1]


        # print(cols_in)

        # drop FILTER and Annotation_mode columns & modify col names
        df_filt.drop("FILTER", axis='columns', inplace=True)
        df_filt.drop("Annotation_mode", axis='columns', inplace=True)
        col_names = list(df_filt.columns)[:7] + ["INFO"]
        df_filt.columns = col_names


        #concateno i tsv
        if i==0:
            df_tot = df_filt
        else:
            df_tot = pd.concat((df_tot, df_filt), axis = 0)

        
    # aggiorno l'id con sample ID
    df_tot['ID_unique'] = df_tot['AnnotSV_ID'] + np.array(df_tot['Samples_ID'], dtype=str)
    s_id = list(df_tot['Samples_ID'].values)

    
    # divido in base al tipo:
    df_del = df_tot.loc[df_tot['SV_type'] == 'DEL'].iloc[:, 2:]
    df_del.index = np.arange(0, len(df_del))
    df_dup = df_tot.loc[df_tot['SV_type'] == 'DUP'].iloc[:, 2:]
    df_dup.index = np.arange(0, len(df_dup))
    df_ins = df_tot.loc[df_tot['SV_type'] == 'INS'].iloc[:, 2:]
    df_ins.index = np.arange(0, len(df_ins))
    if saving_mode:
        #df_tot.to_csv(path_out + 'sv_tot_filt_aggr.csv')
        
        df_del.to_csv(path_out + 'sv_del.csv')
        df_dup.to_csv(path_out + 'sv_dup.csv')
        df_ins.to_csv(path_out + 'sv_ins.csv')
        
        
if __name__ == '__main__':
    
    os.chdir("/hpcshare/genomics/plamagna/sv_backup/")    
    parser = ap.ArgumentParser()
    parser.add_argument("-c", "--configuration", type=str, default='config.ini', help='Parameters for configuration')      
    args = vars(parser.parse_args())  
    config = configparser.ConfigParser()
    config.read(args['configuration'])

    path_in_list = str(ast.literal_eval(config['load_tsv']['path_in_list']))
    path_out = str(ast.literal_eval(config['load_tsv']['path_out']))
    cols_in = list(ast.literal_eval(config['load_tsv']['cols_in']))
    saving_mode = str(ast.literal_eval(config['load_tsv']['saving_mode']))
    
    main(path_in_list, path_out, cols_in, saving_mode)
