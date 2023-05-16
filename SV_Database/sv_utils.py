#!/usr/bin/python3.6

import numpy as np
import pandas as pd


THRESHOLD_SUPERGROUP = 4


def create_intervals(data, var_type):
    data.loc[data['SV_chrom'] == 'X', 'SV_chrom'] =  101 #assegno cromosoma X a un numero
    data.loc[data['SV_chrom'] == 'Y', 'SV_chrom'] =  102 #assegno cromosoma Y a un numero
    data.loc[data['SV_chrom'] == 'M', 'SV_chrom'] =  103 #assegno cromosoma M a un numero
    data['SV_chrom'] = data['SV_chrom'].values.astype(int) # chrom intero
    sorted_data = data.sort_values(by = ['SV_chrom', 'SV_start'], ascending=True) #ordino gli intervalli per cromosoma e poi per start
    # rimodifico chrom 101 --> X
    # sorted_data['SV_chrom'] = data['SV_chrom'].values.astype(str) # chrom intero
    sorted_data.loc[sorted_data['SV_chrom'] == 101, 'SV_chrom'] = 'X' #assegno cromosoma X a un numero
    sorted_data.loc[sorted_data['SV_chrom'] == 102, 'SV_chrom'] = 'Y' #assegno cromosoma X a un numero
    dict_intervals = {}
    
    start = sorted_data.loc[:, 'SV_start'].values
    end = sorted_data.loc[:, 'SV_end'].values
    sv_id = sorted_data.loc[:, 'ID_unique'].values
    ch = sorted_data.loc[:, 'SV_chrom'].values
    l = sorted_data.loc[:, 'SV_length'].values
    ind_i = np.array(sorted_data.index)
    sorted_data.reset_index(inplace=True)
    ind_sort = np.array(sorted_data.index)
    
    list_keys = list(sv_id)
    list_values = [{'chrom':ch[i],'start':start[i],'end':end[i],'length':l[i],'ind_i':ind_i[i],'sv_id':sv_id[i],'ind_sort_i':ind_sort[i]} for i in range(len(sorted_data))]    
    for key, value in zip(list_keys, list_values):
        dict_intervals[key] = value
    return dict_intervals


def recreate_intervals(data):
    #add group length col
    data["group_length"] = data["end"] - data["start"]

    data.loc[data['chrom'] == 'X', 'chrom'] =  101 #assegno cromosoma X a un numero
    data.loc[data['chrom'] == 'Y', 'chrom'] =  102 #assegno cromosoma Y a un numero
    data['chrom'] = data['chrom'].values.astype(int) # chrom intero
    sorted_data = data.sort_values(by = ['chrom', 'start'], ascending=True) #ordino gli intervalli per cromosoma e poi per start
    # rimodifico chrom 101 --> X
    # sorted_data['SV_chrom'] = data['SV_chrom'].values.astype(str) # chrom intero
    sorted_data.loc[sorted_data['chrom'] == 101, 'chrom'] = 'X' #assegno cromosoma X a un numero
    sorted_data.loc[sorted_data['chrom'] == 102, 'chrom'] = 'Y' #assegno cromosoma X a un numero
    dict_intervals = {}
    
    start = sorted_data.loc[:, 'start'].values
    end = sorted_data.loc[:, 'end'].values
    g_id = sorted_data.loc[:, 'group_id'].values
    ch = sorted_data.loc[:, 'chrom'].values
    l = sorted_data.loc[:, 'group_length'].values
    ind_i = np.array(sorted_data.index)
    sorted_data.reset_index(inplace=True)
    ind_sort = np.array(sorted_data.index)
    
    list_keys = list(g_id)
    list_values = [{'chrom':ch[i],'start':start[i],'end':end[i], 'length':l[i], 'ind_i':ind_i[i],'g_id':g_id[i],'ind_sort_i':ind_sort[i]} for i in range(len(sorted_data))]    
    for key, value in zip(list_keys, list_values):
        dict_intervals[key] = value
    return dict_intervals



def overlap(interval_i, group_j, threshold=THRESHOLD_SUPERGROUP):

    threshold=THRESHOLD_SUPERGROUP
    # check the chromosome:
    if interval_i['chrom'] == group_j['chrom']:
        # calculate the overlap size
        # overlap_size = min(interval_i['end'], group_j['inner_end']) - max(interval_i['start'], group_j['inner_start'])
        # # overlap_size = group_j['inner_end'] - group_j['inner_start']
        # interval_length = interval_i['end'] - interval_i['start']


        if (abs(interval_i['start'] - group_j['inner_start']) <= threshold) and (abs(interval_i['end'] - group_j['inner_end']) <= threshold): #se c'è sovrapposizione
            return True
        # if abs(interval_length - overlap_size) <= interval_length / 3: #se c'è sovrapposizione
        #     return True        
        return False
    return False

def clustering(interval_i, group_j):
    # check the chromosome:
    if interval_i['chrom'] == group_j['chrom']:
        # calculate the overlap size
        # overlap_size = min(interval_i['end'], group_j['inner_end']) - max(interval_i['start'], group_j['inner_start'])
        overlap_size = group_j['inner_end'] - group_j['inner_start']
        interval_length = interval_i['end'] - interval_i['start']
        # if ((interval_i['start'] >= group_j['inner_start'] - overlap_size) and interval_i['end'] in range(group_j['inner_start'], group_j['inner_end'])) or \
        # ((interval_i['end'] >= group_j['inner_end'] + overlap_size) and interval_i['start'] in range(group_j['inner_start'], group_j['inner_end'])):
        #     return 0
        if abs(interval_length - overlap_size) <= interval_length / 2: #se c'è sovrapposizione > 50% della lunghezza 
            return True
        return False
    return False



def does_interval_belong_to_group(interval_i, group_i):
        # None: intervals empty, True: yes, False: no
        if not group_i['intervals']: #se il gruppo è vuoto, ovvero intervals == []
            return None
        return overlap(interval_i, group_i) #calcola se c'è overlap tra l'intervallo e il gruppo preesistente: True o False

def does_interval_belong_to_cluster(interval_i, group_i):
        # None: intervals empty, True: yes, False: no
        if not group_i['intervals']: #se il gruppo è vuoto, ovvero intervals == []
            return None
        return clustering(interval_i, group_i) #calcola se c'è overlap tra l'intervallo e il gruppo preesistente: True o False

# crea una lista di dizionari: ogni dict è un supergruppo 
def add_interval(interval, group):
    if group['intervals']: # if the group already exists, so it containts intervals --> update the group infos
        group['inner_start'] = max(group['inner_start'], interval['start'])
        group['inner_end'] = min(group['inner_end'], interval['end'])
        group['outer_start'] = min(group['outer_start'], interval['start'])
        group['outer_end'] = max(group['outer_end'], interval['end'])
    else: # if the group does not exist: creation of a new group with related infos
        group['chrom'] = interval['chrom']
        group['inner_end'] = group['outer_end'] = interval['end']
        group['inner_start'] = group['outer_start'] = interval['start']
        
    group['inner_id'] = str(group['chrom'])+ '_' +str(group['inner_start'])+ '_' +str(group['inner_end'])
    group['intervals'].append(interval)
    # group['intervals'].sort()
    return group #info aggiornate

def add_group(interval, group, update):
    if update == 1:
        if group['intervals']: # if the group already exists, so it containts intervals --> update the group infos
            group['inner_start'] = max(group['inner_start'], interval['start'])
            group['inner_end'] = min(group['inner_end'], interval['end'])
            group['outer_start'] = min(group['outer_start'], interval['start'])
            group['outer_end'] = max(group['outer_end'], interval['end'])
        else: # if the group does not exist: creation of a new group with related infos
            group['chrom'] = interval['chrom']
            group['inner_end'] = group['outer_end'] = interval['end']
            group['inner_start'] = group['outer_start'] = interval['start']
    elif update == 0:       
        if not group['intervals']: # if the group does not exist: creation of a new group with related infos
            group['chrom'] = interval['chrom']
            group['inner_end'] = group['outer_end'] = interval['end']
            group['inner_start'] = group['outer_start'] = interval['start']
    
    group['inner_id'] = str(group['chrom'])+ '-' +str(group['inner_start'])+ '-' +str(group['inner_end'])
    group['intervals'].append(interval)
    # group['intervals'].sort()
    return group #info aggiornate


def group(all_intervals):
    list_intervals = list(all_intervals.values())
    groups = [] # initialize the list that contains all the supergroups
    new_group = {'chrom':0, 'inner_start':0, 'inner_end':0, 'intervals':[]} # initialize new group with its characteristics
    for interval in list_intervals:
        check_interval = does_interval_belong_to_group(interval, new_group)
        if check_interval != False: # if interval belong to group
            new_group = add_interval(interval, new_group) # add interval to new_group and update group infos
        else: # if the interval does not belong to group: closing the group and re-initialize the infos of a new group
            groups.append(new_group)
            new_group = {'chrom':0, 'inner_start':0, 'inner_end': 0, 'intervals':[], 'inner_id': None}
            new_group = add_interval(interval, new_group) # add interval to new_group and update group infos
    groups.append(new_group)
    return groups

def cluster(all_intervals):
    list_intervals = list(all_intervals.values())
    clusters = [] # initialize the list that contains all the supergroups
    new_cluster = {'chrom':0, 'inner_start':0, 'inner_end':0, 'intervals':[]} # initialize new cluster with its characteristics
    for interval in list_intervals:
        check_interval = does_interval_belong_to_cluster(interval, new_cluster)
        if check_interval == True: # if interval belong to group
            new_cluster = add_group(interval, new_cluster, update=1) # add interval to new_group and update group infos
        elif check_interval == 0:
            new_cluster = add_group(interval, new_cluster, update=0) # add interval to new_group
        else: # if the interval does not belong to group: closing the group and re-initialize the infos of a new group
            clusters.append(new_cluster)
            new_cluster = {'chrom':0, 'inner_start':0, 'inner_end': 0, 'intervals':[], 'inner_id': None}
            new_cluster = add_group(interval, new_cluster, update=1) # add interval to new_group and update group infos
    clusters.append(new_cluster)
    return clusters

def dict_to_file(intervals_super_group, name='intervals_super_group'):
    with open(name,'w') as data: 
      data.write(str(intervals_super_group[0:4]))

def exctract_groups(intervals_super_group):
    id_group = 1
    groups = []
    for group in intervals_super_group:
        subset = {key: value for key, value in group.items() if key == 'chrom' or key == 'inner_start' or key == 'inner_end'}
        subset.update({"group_id": id_group})
        groups.append(subset)
        id_group += 1
    return groups

def groups_df(groups, switch=1):
    group_id = []
    chrom = []
    start = []
    end = []
    for group in groups:
        group_id.append(group['group_id'])
        chrom.append(group['chrom'])
        if group['inner_start'] < group['inner_end']:
            start.append(group['inner_start'])
            end.append(group['inner_end'])
        else:
            start.append(group['inner_end'])
            end.append(group['inner_start'])            
    if switch == 1:
        df = pd.DataFrame({'group_id': group_id, 'chrom': chrom, 'start': start, 'end': end})
    elif switch == 2:
        df = pd.DataFrame({'cluster_id': group_id, 'chrom': chrom, 'start': start, 'end': end})
    return df