# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 18:19:33 2022
@author: massale
"""

import numpy as np
import pandas as pd
from scipy.spatial import distance
from sklearn.cluster import DBSCAN
#import matplotlib.pyplot as plt
import time

THRESHOLD_SUPERGROUP = 4

# =============================================================================
# #create dictionary of all intervals
# def create_intervals(data, var_type):
#     data.loc[data['SV_chrom'] == 'X', 'SV_chrom'] =  101 #assegno cromosoma X a un numero
#     data.loc[data['SV_chrom'] == 'Y', 'SV_chrom'] =  102 #assegno cromosoma Y a un numero
#     data['SV_chrom'] = data['SV_chrom'].values.astype(int) # chrom intero
#     sorted_data = data.sort_values(by = ['SV_chrom', 'SV_start'], ascending=True) #ordino gli intervalli per cromosoma e poi per start
#     # rimodifico chrom 101 --> X
#     # sorted_data['SV_chrom'] = data['SV_chrom'].values.astype(str) # chrom intero
#     sorted_data.loc[sorted_data['SV_chrom'] == 101, 'SV_chrom'] = 'X' #assegno cromosoma X a un numero
#     sorted_data.loc[sorted_data['SV_chrom'] == 102, 'SV_chrom'] = 'Y' #assegno cromosoma X a un numero
#     all_intervals = {}
#     
#     for i, ind in enumerate(sorted_data.index):
#         start = sorted_data.loc[ind, 'SV_start']
#         sv_id = sorted_data.loc[ind, 'ID_unique']
#         ch = sorted_data.loc[ind, 'SV_chrom']
#         l = sorted_data.loc[ind, 'SV_length']
#         if var_type == 'ins':
#             end = np.int(start + l)
#         else:
#             end = sorted_data.loc[ind, 'SV_end']
#         all_intervals[sv_id] = {'chrom':ch,'start':start,'end':end,'length':l,'ind_i':ind,'sv_id':sv_id,'ind_sort_i': i}
#     return all_intervals
# =============================================================================

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

# =============================================================================
# 
# #create dictionary of all intervals
# def create_intervals2(data, var_type):
#     
#     data.loc[data['SV_chrom'] == 'X', 'SV_chrom'] =  101 #assegno cromosoma X a un numero
#     data.loc[data['SV_chrom'] == 'Y', 'SV_chrom'] =  102 #assegno cromosoma Y a un numero
#     data['SV_chrom'] = data['SV_chrom'].values.astype(int) # chrom intero
#     sorted_data = data.sort_values(by = ['SV_chrom', 'SV_start'], ascending=True) #ordino gli intervalli per cromosoma e poi per start
#     # sorted_data['SV_chrom'] = data['SV_chrom'].values.astype(str) # chrom intero
#     sorted_data.loc[sorted_data['SV_chrom'] == 101, 'SV_chrom'] = 'X' #assegno cromosoma X a un numero
#     sorted_data.loc[sorted_data['SV_chrom'] == 102, 'SV_chrom'] = 'Y' #assegno cromosoma X a un numero
#     
#     sorted_data['ind_i'] = np.array(sorted_data.index)
#     sorted_data.reset_index(inplace=True)
#     sorted_data['ind_sort_i'] = np.array(sorted_data.index)
#     sorted_data.rename(columns={"SV_chrom": "chrom", "SV_start": "start", "SV_end": "end", "SV_length": "length"}, inplace=True)
#     sel_cols = ["chrom","start","end","length","ind_i","ind_sort_i"]
#     sorted_data.set_index('ID_unique', inplace=True)
#     sorted_data = sorted_data.loc[:, sel_cols]
#     
#     sorted_data_T = sorted_data.T
#     all_intervals = sorted_data_T.to_dict('dict')
#     return all_intervals
# =============================================================================


# test the overlapping between the group_i and a new_interval
# def overlap(interval_i, group_j):
    # check the chromosome:
#    if interval_i['chrom'] == group_j['chrom']:
#        # calculate the overlap size
#        overlap_size = max(interval_i['end'], group_j['inner_end']) - min(interval_i['start'], group_j['inner_start'])
#        # length of the group: it represents the size of the group!
#        group_size = group_j['inner_end'] - group_j['inner_start']
#        # total size is the sum of the two intervals' length that I am considering: interval_i and group_i 
#        total_size = interval_i['length'] + group_size
#        if total_size > overlap_size: #se c'è sovrapposizione
#            return True
#    return False

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

# def overlap(interval_i, group_j):
#     # check the chromosome:
#     if interval_i['chrom'] == group_j['chrom']:
#         if interval_i['chrom'] 
#             return True
#     return False

'''
def overlap_50(interval_i, interval_j):
    overlap_size = min(interval_i['end'], interval_j['end']) - max(interval_i['start'], interval_j['start'])
    size_i = interval_i['length']
    size_j = interval_j['length']
    if overlap_size >= 0.5*size_i and overlap_size >= 0.5*size_j:
        return True
    else:
        return False

def overlap_50_incluster(count_ov50, labels_, super_group):
    for el in np.unique(labels_):
        if el >= 0 and el < 100:
            bool_ind = labels_ == el
            intervals_ = np.array(super_group['intervals'])
            intervals_el = intervals_[bool_ind]
            interval_i = intervals_el[0]
            interval_j = intervals_el[-1]
            if overlap_50(interval_i, interval_j):
                count_ov50 += 1
    return count_ov50

   
def overlapping_clusters(super_group, labels_, counts_tot):
    sg_int= list(super_group['intervals'])
    if len(np.unique(labels_)) >= 1:
        intervals_clusters = []
        list_clusters = [el for el in np.unique(labels_) if el != -1]
        for l in list_clusters:
            ind_l = np.where(labels_ == l)[0]
            sg_int = np.array(sg_int)
            intervals_clusters.append(sg_int[ind_l])
    
    l = len(intervals_clusters)
    #print('- EPS = {}: supergruppo {} con {} intervals'.format(eps, ind, l))
    if l > 1:
        lista_ind_comb=[]
        for v in np.arange(l):
            for w in np.arange(l):
                if v!=w and [w, v] not in lista_ind_comb:
                        lista_ind_comb.append([v, w])
                
        for k, j in lista_ind_comb:
            in_start1 = max([elem['start'] for elem in intervals_clusters[k]])
            in_end1 = min([elem['end'] for elem in intervals_clusters[k]])
            in_start2 = max([elem['start'] for elem in intervals_clusters[j]])
            in_end2 = min([elem['end'] for elem in intervals_clusters[j]])
            cluster1 = {'start': in_start1, 'end':in_end1, 'length': in_end1-in_start1}
            cluster2 = {'start': in_start2, 'end':in_end2, 'length': in_end2-in_start2}
            overlap = overlap_clusters(cluster1, cluster2)
            l1 = cluster1['length']
            l2 = cluster2['length']
            
            if overlap >= l1*0.5 and overlap >= l2*0.5:
                counts_tot += 1 #se overlap > 50% reciproco tra i clusters
    
    return counts_tot

# test the overlapping between the group_i and a new_interval
def overlap_clusters(cluster_i, cluster_j):
    # calculate the overlap size
    overlap_size = min(cluster_i['end'], cluster_j['end']) - max(cluster_i['start'], cluster_j['start'])
    return overlap_size
'''

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



'''
def count_clusters_tot(df2):
    #df2['id_supG'] = pd.to_numeric(df2['id_supG'])
    #df2['cluster_id'] = pd.to_numeric(df2['cluster_id'])
    bool_111 = np.any([[df2['cluster_id'].values == -111], [df2['cluster_id'].values == '-111']], axis=0).ravel()
    bool_1 = np.any([[df2['cluster_id'].values == -1], [df2['cluster_id'].values == '-1']], axis=0).ravel()
    bool_all = ~(bool_111 + bool_1)
    df2.loc[bool_111, 'cluster_id_tot'] = -111
    df2.loc[bool_1, 'cluster_id_tot'] = -1
    if type(df2.loc[bool_all, 'id_supG']) == str:
        str_idsupG = df2.loc[bool_all, 'id_supG'].values
    else:
        str_idsupG = df2.loc[bool_all, 'id_supG'].values.astype('str')
    if type(df2.loc[bool_all,'cluster_id']) == str:
        str_clId = df2.loc[bool_all,'cluster_id'].values
    else:
        str_clId = df2.loc[bool_all,'cluster_id'].values.astype('str')
    df2.loc[bool_all, 'cluster_id_tot'] = np.core.defchararray.add(str_idsupG,str_clId)
    a = np.array(df2.loc[bool_all, 'cluster_id_tot'], dtype=float)
    #b = np.array(a, dtype=float)
    df2.loc[bool_all, 'cluster_id_tot'] = np.digitize(a, bins = np.unique(a), right = True)
    return df2

def compute_frequency(df2, n_samples):
    freqs = df2['cluster_id_tot'].value_counts()
    bool_111 = df2['cluster_id_tot'] == -111
    bool_1 = df2['cluster_id_tot'] == -1
    bool_all = ~(bool_111 + bool_1)
    df2.loc[bool_111, 'freq'] = 1/n_samples
    df2.loc[bool_1, 'freq'] = 1/n_samples
    freqs.index = np.array(freqs.index, dtype=int)
    b = df2.loc[bool_all, 'cluster_id_tot'].values
    c = freqs[b].values.reshape((len(b), 1))
    df2.loc[bool_all, 'freq'] = c/n_samples
    return df2

def create_table_supergroup(ind, df, super_group):
    # create table of supergroup infos:
    df.loc[ind, 'inner_id'] = super_group['inner_id']
    df.loc[ind, 'inner_start'] = super_group['inner_start']
    df.loc[ind, 'inner_end'] = super_group['inner_end']
    df.loc[ind, 'chrom'] = super_group['chrom']
    df.loc[ind, 'n_int'] = len(super_group['intervals'])
    df.loc[ind, 'outer_start'] = super_group['outer_start']
    df.loc[ind, 'outer_end'] = super_group['outer_end']
    return df
    
def create_table_results(data,n_clusters_tot,len_sg,n_noise,n_noise_,rate, rate_tot, count_ov50,counts_tot):
    df_tot = pd.DataFrame(np.zeros((8, 1)))
    df_tot.index = ['N intervals', 'N supergroups', 'N clusters', 'N unique_intervals', 'CountOv50_in_clusters', 'CountsOv50_bw_clusters', 'Rate1', 'Rate2']
    df_tot.loc['N intervals', 0] = data.shape[0]
    df_tot.loc['N clusters', 0] = n_clusters_tot
    df_tot.loc['N supergroups', 0]  = len_sg
    df_tot.loc['N unique_intervals', 0] = n_noise
    df_tot.loc['N supergroups of 1 interval', 0] = n_noise_
    df_tot.loc['Rate1=Ncl/Nsup', 0] = rate
    df_tot.loc['Rate2=Ncl+singleCl/Nsup', 0] = rate_tot
    df_tot.loc['CountOv50_in_clusters', 0] = count_ov50
    df_tot.loc['CountsOv50_bw_clusters', 0] = counts_tot
    return df_tot

def create_table_intervals(ind, res_, df, super_group, labels_, ind_intervals, len_int):
    n_int = int(df.loc[ind, 'n_int'])
    inn_id = np.array(np.full(shape=(n_int), fill_value = super_group['inner_id']), ndmin=2).T
    cl_id =  np.array(labels_, ndmin=2).T
    start_ =  np.array([el['start'] for el in super_group['intervals']], ndmin=2).T
    sv_id =  np.array([el['sv_id'] for el in super_group['intervals']], ndmin=2).T
    end_ =  np.array([el['end'] for el in super_group['intervals']], ndmin=2).T
    chrom = np.array(np.full(shape=(n_int), fill_value = super_group['chrom']), ndmin=2).T
    ind_sg = np.array(np.full(shape=(n_int), fill_value = ind), ndmin=2).T       
    df2_arr = np.concatenate([sv_id, chrom, start_, end_, inn_id, ind_sg, cl_id], axis=1)
    res_[ind_intervals : ind_intervals+len_int, :] = df2_arr
    
    return res_
#--------------------------------------------------
#             CLUSTERING WITH DBSCAN
#--------------------------------------------------

def cluster(super_group, eps=0.5, min_samples=2, njobs=1):
    if len(super_group['intervals']) < 2:
        clusters = [-111]
        n_cl_nonoise = -111
        C = []
    else:    
        A = list([(i['start'], i['length']) for i in super_group['intervals']])
        C = DBSCAN(eps=eps, metric=interval_distance_for_cdist, min_samples=min_samples, n_jobs=njobs).fit(A) 
        clusters = C.labels_
        n_cl = np.unique(clusters)
        n_cl_nonoise = len([elem for elem in n_cl if elem != -1])
    return clusters, n_cl_nonoise, C

def dict_clustering_(model):
    clust_model = {}
    if model != []:
        clust_model['components'] = model.components_
        clust_model['core_sample_indices'] = model.core_sample_indices_
        clust_model['eps'] = model.eps
        clust_model['labels'] = model.labels_
    return clust_model

def cluster2(list_intervals, eps=0.5, min_samples=2):
    A = list([(i['start'], i['length']) for i in list_intervals])
    C = DBSCAN(eps=eps, metric=interval_distance_for_cdist, min_samples=min_samples, n_jobs=-1).fit(A)       
    return C.labels_

def interval_distance_for_cdist(u, v):
    # customised distance measure
    # 1 - (shared size / merged size)
    overlap_size = min(u[0]+u[1], v[0]+v[1]) - max(u[0], v[0])
    merged_size = max(u[0]+u[1], v[0]+v[1]) - min(u[0], v[0])
    if overlap_size < 0:
        overlap_size = 0
    distance = 1 - (overlap_size/merged_size)
    return distance

def dbscan_predict(dbscan_model, interval_i, metric=interval_distance_for_cdist):
    # Result is noise by default
    y_new = np.ones(shape=1, dtype=int)*-1
    # Iterate all input samples for a label
    #for j, x_new in enumerate(interval_i):
        # Find a core sample closer than EPS
    x = [interval_i['start'], interval_i['length']]
    if dbscan_model != {}:
        for i, x_core in enumerate(dbscan_model['components']): 
            if metric(x, x_core) < dbscan_model['eps']:
                # Assign label of x_core to x_new
                y_new = dbscan_model['labels'][dbscan_model['core_sample_indices'][i]]
                break
    return y_new

def loop_supergroup_clust_save(ind_intervals, intervals_super_group, df, eps, njobs, count_ov50, num_clust, clustering_models, file_2, res_, counts_tot):
    for ind, super_group in enumerate(intervals_super_group):
        start_iter = time.time()
        df = create_table_supergroup(ind, df, super_group)        
        len_int = len(super_group['intervals'])
       
        ## timing #1
        tic = time.time()
        ### CLUSTER SUPERGROUPS WITH DBSCAN:
        labels_, n_cl_nonoise, model = cluster(super_group, eps=eps, njobs=njobs)
        toc = time.time()
        file_2.write('\nsupergroup:{}/{}'.format(ind, len(intervals_super_group)))
        file_2.write('\ntime clustering:{}'.format(toc-tic)) 
        
        ## timing #2
        tic = time.time()
        #conteggio clusters che hanno effettivamente intervalli che si overlappano + del 50%
        count_ov50 = overlap_50_incluster(count_ov50, labels_, super_group)
        toc = time.time()
        file_2.write('\ntime countOv50:{}'.format(toc-tic))

        ### SAVING CLUSTERING MODELS FOR FUTURE PREDICTION
        clust_model_i = dict_clustering_(model)
        clustering_models.append(clust_model_i)
        
        ### PLOT SUPERGROUP - INTERVALS COLORED BY CLUSTERS
        #if os.path.isdir('figs/'+var_type+'/clusters/') != True:
            #os.makedirs('figs/'+var_type+'/clusters/')
        #plot_clusters_supergroups(labels_, ind, super_group, path='figs/'+var_type+'/clusters/')
        
        ### COUNT CLUSTERS FOR EACH SUPERGROUP
        num_clust.append(n_cl_nonoise)
        bool_nonoise = np.array(num_clust) != -111
        
        tic = time.time()
        #fill the df2 dataframe:
        res_ = create_table_intervals(ind, res_, df, super_group, labels_, ind_intervals, len_int)
        ind_intervals = ind_intervals + len_int 
        toc = time.time()
        file_2.write('\ntime count clusters_2:{}'.format(toc-tic))
   
        ## timing #4
        tic = time.time()
        # count overlapping clusters
        counts_tot = overlapping_clusters(super_group, labels_, counts_tot)
        toc = time.time()
        file_2.write('\ntime count ov. clusters:{}'.format(toc-tic))

        end_iter = time.time()
        file_2.write('\ntime tot iteration supergroup:{}'.format(end_iter-start_iter))
        
    return df, res_, counts_tot, num_clust
#--------------------------------------------------
#                FUNCTIONS OF PLOT
#--------------------------------------------------

# =============================================================================
# def plot_supergoup(ind, supergroup, path):
#     fig = plt.figure(figsize=[7,5])
#     for i, interval in enumerate(supergroup['intervals']):
#         xx = [interval['start'], interval['end']]
#         yy = np.linspace(i+1,i+1,num=len(xx))
#         plt.plot(xx,yy, linewidth=3)
#         plt.scatter(xx,yy, s=50)
#     plt.xlabel('positions', fontsize=13)
#     plt.ylabel('intervals', fontsize=13)
#     plt.yticks([])
#     plt.grid(linewidth=0.4)
#     plt.title('supergroup '+ str(ind) + ' - chr '+ str(supergroup['chrom']), fontsize=16)
#     plt.ylim([0, yy[0]+1])
#     plt.close()
#     fig.savefig(path + 'plot_supergroup_' +str(ind)+ '_chr' +str(supergroup['chrom'])+ '.png')
#     
# =============================================================================
    
# =============================================================================
# def plot_clusters(labels_, list_all_intervals, path):
#     for cl in np.unique(labels_):
#         fig = plt.figure(figsize=[7,5])
#         indices = labels_== cl
#         intervals_cl = np.array(list_all_intervals)[indices]
#         
#         for i, interval in enumerate(intervals_cl):
#             xx = [interval['start'], interval['end']]
#             yy = np.linspace(i+1, i+1, num=len(xx))
#             plt.plot(xx, yy, linewidth=3)
#             plt.scatter(xx, yy, s=50)
#         
#         plt.xlabel('positions', fontsize=13)
#         plt.ylabel('intervals', fontsize=13)
#         plt.yticks([])
#         plt.grid(linewidth=0.4)
#         plt.title('cluster '+ str(cl) + ' - chr '+ str(interval['chrom']), fontsize=16)
#         plt.ylim([0, yy[0]+1])
#         plt.close()
#         fig.savefig(path + 'plot_cluster_' +str(cl)+ '_chr' +str(interval['chrom'])+ '.png')
#         
#              
# =============================================================================
# =============================================================================
# def plot_clusters_supergroups(labels, ind, supergroup, path):
#     fig = plt.figure(figsize=[7, 5])
#     colors = ['tomato', 'limegreen', 'deepskyblue', 'purple', 'grey']
#     cc = [colors[lab] if lab != -111 else 'grey' for lab in labels]
#     for i, interval in enumerate(supergroup['intervals']):
#         xx = [interval['start'], interval['end']]
#         yy = np.linspace(i+1,i+1,num=len(xx))
#         plt.plot(xx, yy, linewidth=3, c=cc[i])
#         plt.scatter(xx, yy, s=50, c=cc[i])
#     plt.xlabel('positions', fontsize=13)
#     plt.ylabel('intervals', fontsize=13)
#     plt.yticks([])
#     plt.grid(linewidth=0.4)
#     plt.title('Clusters of supergroup '+ str(ind) + ' - chr '+ str(supergroup['chrom']), fontsize=16)
#     plt.ylim([0, yy[0]+1])
#     plt.close()
#     if len(supergroup['intervals']) == 1:
#         print(-111)
#         #fig.savefig(path + 'cl_111/' + 'plot_clusters_supergroup_' +str(ind)+ '_chr' +str(supergroup['chrom'])+ '.png')
#     else:
#         if np.all(np.unique(labels) == [-1]):
#             print(-1)
#             #fig.savefig(path +'cl_-1/'+ 'plot_clusters_supergroup_' +str(ind)+ '_chr' +str(supergroup['chrom'])+ '.png')
#         elif np.all(np.unique(labels) == [0]) or np.all(np.unique(labels)==[-1, 0]):
#             #fig.savefig(path +'1cluster/'+ 'plot_clusters_supergroup_' +str(ind)+ '_chr' +str(supergroup['chrom'])+ '.png')
#             np.save(path + '1cluster/labels_cl_'+str(ind), labels)
#         else:
#             fig.savefig(path +'2cluster/'+ 'plot_clusters_supergroup_' +str(ind)+ '_chr' +str(supergroup['chrom'])+ '.png')
#         
# =============================================================================
'''