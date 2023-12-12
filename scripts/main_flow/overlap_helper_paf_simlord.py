#%%
from readpaf import parse_paf
import numpy as np
import pandas as pd

def process_paf(paf_path, pen1, pen2):

    with open(paf_path,'r') as handle:
        df = parse_paf(handle, dataframe=True)

    df_n = df[['query_name', 'query_length', 'query_start', 'query_end', 
        'target_name', 'target_length', 'target_start', 'target_end',
        'residue_matches', 'alignment_block_length', 'NM','AS']]
    # writer = pd.ExcelWriter('df.xls', engine='openpyxl')
    # df.to_excel(writer)
    # writer.close()
    # writer = pd.ExcelWriter('df_n.xls', engine='openpyxl')
    # df_n.to_excel(writer)
    # writer.close()
    # df_n.loc[:,'query_name'] = df_n['query_name'].astype(np.int64).values.copy()
    # df_n.loc[:,'target_name'] = df_n['target_name'].astype(np.int64).values.copy()
    # df_n[['query_name','target_name']] = df_n[['query_name','target_name']].apply(pd.to_numeric)
    df_n = df_n.astype({'query_name':np.int64, 'target_name':np.int64})

    df_n1 = df_n[df_n['query_start']<=20]
    df_n2 = df_n[df_n['target_start']<=20]
    df_n = pd.concat([df_n1,df_n2],ignore_index=True)
    df_n = df_n[~((df_n['target_start']<=20) & (df_n['query_start']<=20))]

    # 至少要少于read长度
    df_n = df_n[df_n['alignment_block_length'] <= df_n['query_length']]
    df_n = df_n[df_n['alignment_block_length'] <= df_n['target_length']]
    # 至少满足比其中一个少
    con1 = df_n['alignment_block_length'] <= df_n['query_length']*0.95
    con2 = df_n['alignment_block_length'] <= df_n['target_length']*0.95
    df_n = df_n[con1 + con2]

    # df_n = df_n[df_n['NM']<=250]

    # 至少满足比其中一个的0.3倍长
    con1 = df_n['alignment_block_length'] > df_n['query_length']*0.85
    con2 = df_n['alignment_block_length'] > df_n['target_length']*0.85
    df_n = df_n[con1 + con2]

    df_n.sort_values(by=['AS','NM'],ascending=[False,True],inplace=True)
    df_n.drop_duplicates(subset=['query_name','target_name'],keep='first',inplace=True,ignore_index=True)

    reads_n = max(df_n['query_name'].max(), df_n['target_name'].max())
    ovlp_mat = np.zeros((reads_n+1, reads_n+1), dtype=np.float64)

    for i in range(1, reads_n+1):

        ri = i

        idf = df_n[df_n['query_name']==ri]
        if not idf.empty:
            for _, rec in idf.iterrows():
                j = int(rec['target_name'])
                # if i == 154 and j == 523:
                #     stop = 0
                if ovlp_mat[j,i]!=0 or ovlp_mat[i,j]!=0:
                    continue
                NM = rec['NM']
                # abl = rec['alignment_block_length']
                # scorep = avg_rl*(abl-pen*NM)/(abl+pen*NM) + rec['alignment_block_length'] - pen*NM
                # scorep = avg_rl*(1 - pen*NM/rec['alignment_block_length'])
                # scorep = rec['alignment_block_length'] - pen * NM
                scorep = round(rec['alignment_block_length'] - pen1 * NM - pen2 * NM/rec['alignment_block_length'])

                if rec['query_start'] <= 20:
                    ovlp_mat[j,i] = -scorep
                else:
                    ovlp_mat[i,j] = -scorep
        
        idf = df_n[df_n['target_name']==ri]
        if not idf.empty:
            for _,rec in idf.iterrows():
                j = int(rec['query_name'])
                # if i == 154 and j == 523:
                #     stop = 0
                if ovlp_mat[j,i]!=0 or ovlp_mat[i,j]!=0:
                    continue
                NM = rec['NM']
                # abl = rec['alignment_block_length']
                # scorep = avg_rl*(abl-pen*NM)/(abl+pen*NM) + rec['alignment_block_length'] - pen*NM
                # scorep = avg_rl*(1 - pen*NM/rec['alignment_block_length'])
                # scorep = rec['alignment_block_length'] - pen*NM
                scorep = round(rec['alignment_block_length'] - pen1 * NM - pen2 * NM/rec['alignment_block_length'])
                
                if rec['query_start'] <= 20:
                    ovlp_mat[i,j] = -scorep
                else:
                    ovlp_mat[j,i] = -scorep
    
    return ovlp_mat, df_n, df
#%%
if __name__=='__main__':
    import os
    paf_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/data/MHC/output/v30K/' + 'aln.paf'
    overlap_, paf_df = process_paf(paf_path, 2, 100000)


# %%
