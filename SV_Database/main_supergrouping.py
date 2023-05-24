""" LAUNCH
echo "singularity clustering3.sif python main_supergrouping.py -c config.ini" | qsub 
"""

import os
import ast
import time
import configparser
import pandas as pd
import argparse as ap
from sv_utils import *


def main(path_in, var_type, path_out):    

    # import data csv
    data = pd.read_csv(path_in + 'sv_' +var_type+ '.csv')
    all_intervals = create_intervals(data, var_type)
    
    # create supergroups of intervals based on overlapping of intervals
    intervals_super_group = group(all_intervals)
    intervals_super_group_df = pd.DataFrame(intervals_super_group)
    intervals_super_group_df.to_csv(path_out + 'intervals_super_group_' + var_type)
    groups = exctract_groups(intervals_super_group)
    with open('groups_' + var_type,'w') as data: 
      data.write(path_out + str(groups))

    groups_dataframe = groups_df(groups)    
    groups_dataframe.to_csv(path_out + 'groups_df_' + var_type, sep = '\t', index = False)

  
if __name__ == '__main__':
    
    parser = ap.ArgumentParser()
    parser.add_argument("-c", "--configuration", type=str, default='config.ini', help='Parameters for configuration')    
    args = parser.parse_args()
    args = vars(parser.parse_args())
    config = configparser.ConfigParser()
    config.read(args['configuration'])      
    path_in = str(ast.literal_eval(config['clust']['path_in']))
    path_out = str(ast.literal_eval(config['clust']['path_out']))
    var_type = str(ast.literal_eval(config['clust']['var_type']))
    cols_in = str(ast.literal_eval(config['load_tsv']['cols_in']))

    if os.path.isdir(path_out) != True:
        os.mkdir(path_out)

    main(path_in, var_type, path_out)
