#%%
from readpaf import parse_paf
import numpy as np
import pandas as pd

def process_paf(paf_path, reads_n, pen, avg_rl):

    with open(paf_path,'r') as handle:
        df = parse_paf(handle, dataframe=True)

    df_n = df[['query_name', 'query_length', 'query_start', 'query_end', 
        'target_name', 'target_length', 'target_start', 'target_end',
        'residue_matches', 'alignment_block_length', 'NM','AS']]
    # tmp = df_n.copy()
    df_n1 = df_n[df_n['query_start']<=20]
    df_n2 = df_n[df_n['target_start']<=20]
    df_n = pd.concat([df_n1,df_n2],ignore_index=True)
    df_n = df_n[~((df_n['target_start']<=20) & (df_n['query_start']<=20))]

    # 至少要少于read长度
    df_n = df_n[df_n['alignment_block_length'] <= df_n['query_length']]
    df_n = df_n[df_n['alignment_block_length'] <= df_n['target_length']]
    # 至少满足比其中一个少0.9
    con1 = df_n['alignment_block_length'] < df_n['query_length']*0.9
    con2 = df_n['alignment_block_length'] < df_n['target_length']*0.9
    df_n = df_n[con1 + con2]

    # 至少满足比其中一个的0.3倍长
    con1 = df_n['alignment_block_length'] > df_n['query_length']*0.3
    con2 = df_n['alignment_block_length'] > df_n['target_length']*0.3
    df_n = df_n[con1 + con2]

    df_n.sort_values(by=['AS','NM'],ascending=[False,True],inplace=True)
    df_n.drop_duplicates(subset=['query_name','target_name'],keep='first',inplace=True,ignore_index=True)

    ovlp_mat = np.zeros((reads_n+1, reads_n+1), dtype=np.float64)

    for i in range(1, reads_n+1):

        ri = f'read_{i}'

        idf = df_n[df_n['query_name']==ri]
        if not idf.empty:
            for _, rec in idf.iterrows():
                j = int(rec['target_name'].split('_')[-1])
                if ovlp_mat[j,i]!=0 or ovlp_mat[i,j]!=0:
                    continue
                NM = rec['NM']
                # abl = rec['alignment_block_length']
                # 给完全匹配的片段合适的bonus
                # 相同identity选较长的
                # scorep = avg_rl*(abl-pen*NM)/(abl+pen*NM) + rec['alignment_block_length'] - pen*NM
                # scorep = avg_rl*(1 - pen*NM/rec['alignment_block_length'])
                scorep = rec['alignment_block_length'] - pen*NM

                if rec['query_start'] <= 20:
                    ovlp_mat[j,i] = -scorep
                else:
                    ovlp_mat[i,j] = -scorep
        
        idf = df_n[df_n['target_name']==ri]
        if not idf.empty:
            for _,rec in idf.iterrows():
                j = int(rec['query_name'].split('_')[-1])
                if ovlp_mat[j,i]!=0 or ovlp_mat[i,j]!=0:
                    continue
                NM = rec['NM']
                # abl = rec['alignment_block_length']
                # 给完全匹配的片段合适的bonus
                # 相同identity选较长的
                # scorep = avg_rl*(abl-pen*NM)/(abl+pen*NM) + rec['alignment_block_length'] - pen*NM
                # scorep = avg_rl*(1 - pen*NM/rec['alignment_block_length'])
                scorep = rec['alignment_block_length'] - pen*NM
                
                if rec['query_start'] <= 20:
                    ovlp_mat[i,j] = -scorep
                else:
                    ovlp_mat[j,i] = -scorep
    
    return ovlp_mat
        



# %%
