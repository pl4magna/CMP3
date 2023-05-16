# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 10:32:09 2022
@author: massale, plamagna
"""

""" LAUNCH
echo "singularity exec -B /hpcshare/genomics/plamagna/,/hpcshare/genomics/afant/PPMI/SV_CONSENSUS/ /hpcshare/genomics/plamagna/sv_backup/container/clustering3.sif python /hpcshare/genomics/plamagna/sv_backup/code/main_supergrouping.py -c /hpcshare/genomics/plamagna/sv_backup/code/config.ini" | qsub -q fatnodes -e /hpcshare/genomics/plamagna/sv_backup/main_supergrouping.err
"""

import os
import ast
import time
# import numpy as np
import configparser
import pandas as pd
import argparse as ap
from sv_utils import *

# def main(path_in, var_type, eps, n_samples, path_out, njobs):
def main(path_in, var_type, path_out):    

    # file_2 = open('log_times_'+var_type+'_eps_' +str(eps).replace('.', '')+ '.txt', "a+")
    # start_ = time.time()
    
    # import data csv
    data = pd.read_csv(path_in + 'sv_' +var_type+ '.csv')

    # create dictionary of intervals with the infos (chr, start, end, length, id_unique,...)
    all_intervals = create_intervals(data, var_type)
    # print(all_intervals.values()[0])
    # all_intervals_df = pd.DataFrame(all_intervals)
    # all_intervals_df.to_csv(path_out + 'intervals')
    
    # create supergroups of intervals based on overlapping of intervals
    intervals_super_group = group(all_intervals)
    intervals_super_group_df = pd.DataFrame(intervals_super_group)
    intervals_super_group_df.to_csv(path_out + 'intervals_super_group_' + var_type)
    # dict_to_file(intervals_super_group)
    groups = exctract_groups(intervals_super_group)
    with open('groups_' + var_type,'w') as data: 
      data.write(path_out + str(groups))

    groups_dataframe = groups_df(groups)    
    groups_dataframe.to_csv(path_out + 'groups_df_' + var_type, sep = '\t', index = False)


#     groups_dataframe = pd.read_csv(path_out + 'subset_5k', sep = '\t', index_col = False)

#     all_intervals = recreate_intervals(groups_dataframe)
# #    print(all_intervals.values())
#     intervals_cluster = cluster(all_intervals)

#     intervals_cluster_df = pd.DataFrame(intervals_cluster)
#     intervals_cluster_df.head(500).to_csv(path_out + 'intervals_cluster')
#     # dict_to_file(intervals_super_group)
#     clusters = exctract_groups(intervals_cluster)

#     groups_dataframe = groups_df(clusters, switch=2)    
#     groups_dataframe.to_csv(path_out + 'clusters_df', sep = '\t', index = False)

  
if __name__ == '__main__':
    
    # print(os.getcwd())
    parser = ap.ArgumentParser()
    parser.add_argument("-c", "--configuration", type=str, default='config.ini', help='Parameters for configuration')    
    # parser.add_argument("-e", "--eps", type=float, required=True, default = 0.1, help = 'Eps_value - DBSCAN parameter')
    args = parser.parse_args()    #command line parser
    # eps_value = args.eps
    args = vars(parser.parse_args())   #config_file parser
    config = configparser.ConfigParser()
    config.read(args['configuration'])      
    # path_in = str(ast.literal_eval(config['clust']['path_in']))
    # path_out = str(ast.literal_eval(config['clust']['path_out']))
    # var_type = str(ast.literal_eval(config['clust']['var_type']))
    # n_samples = int(ast.literal_eval(config['clust']['nsamples']))
    # njobs = int(ast.literal_eval(config['clust']['njobs']))
    path_in = str(ast.literal_eval(config['clust']['path_in']))
    path_out = str(ast.literal_eval(config['clust']['path_out']))
    var_type = str(ast.literal_eval(config['clust']['var_type']))
    cols_in = str(ast.literal_eval(config['load_tsv']['cols_in']))
    # saving_mode = int(ast.literal_eval(config['load_tsv']['saving_mode']))
    # print('eps = {}'.format(eps_value))
    if os.path.isdir(path_out) != True:
        os.mkdir(path_out)

    # main(path_in, var_type, eps_value, n_samples, path_out, njobs)
    main(path_in, var_type, path_out)
