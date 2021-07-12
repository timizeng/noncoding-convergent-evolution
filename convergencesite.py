import pandas as pd 
import numpy as np 
from collections import defaultdict

df = pd.read_csv('samples.tsv',sep='\t')
df['up'] = df['start'].apply(lambda x:int(x)-200)
df['down'] = df['start'].apply(lambda x:int(x) +200)
df.sort_values('motif_alt_id',inplace=True)
all_species = ['or.sativa', 'eu.grandis', 'po.trichocarpa','se.indicum']
any_species = ['av.marina', 'rh.apiculata','so.alba']
candidates = set(df['motif_alt_id'])
out_dfs = []

def mark_same_motifs(intervals):
    ov = defaultdict(list)
    for index,interval in enumerate(intervals):
        if index ==0:
            for sub in intervals[index+1:]:
                if sub[0] <= interval[1]:
                    ov[str(index)].append(sub)
        elif index > 0:
            for sub in intervals[index+1:]:
                if sub[0] <= interval[1]:
                    keys = list(ov.keys())
                    total_items = [item for i in keys for item in ov[i]]
                    if sub not in total_items and interval not in total_items:
                        ov[str(index)].append(sub)
    return ov 
for i in candidates:
    sub_df = df[df['motif_alt_id'] == i]
    if sub_df.shape[0] >=4:
        sub_df = sub_df.sort_values('start')
        sub_df.reset_index(inplace=True)
        intervals = [[up,down] for up,down in zip(sub_df['up'],sub_df['down'])]
        overlap_los = mark_same_motifs(intervals)
        for i in overlap_los.keys():
            same_motif_df = sub_df.iloc[int(i):int(i)+len(overlap_los[i])+1]
            if same_motif_df.shape[0] > 4:
                spe = set(same_motif_df['species'])
                print(same_motif_df.shape[0])
                print(set(same_motif_df['species']))
                if len(spe & set(all_species) ) >=4 and len(spe & set(any_species)) <=1:
                    out_dfs.append(same_motif_df)