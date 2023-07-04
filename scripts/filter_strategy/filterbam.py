#%%
# 将包含式的reads去掉，降低深度
import pysam
import os
from copy import deepcopy

def filter(pen:int, vn:int, out_dir):
    num_rds = []
    print("filtering...")
    for i in range(1, 3):
        inbam = pysam.AlignmentFile(out_dir + f'/r{i}.sort.bam', 'rb')
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        newbam_name = out_dir + f'/readsh{i}.bam'
        newbam = pysam.AlignmentFile(newbam_name,'wb',template=inbam)
        ALL = set()
        start, end = 0, 30000
        target = 0
        call_score = lambda x: - pen*abs(x.pos - target) - 100*x.get_tag('NM')
        ref_name = inbam.references[0]
        bam_list = list(inbam.fetch(ref_name, start, end))
        ii = 0
        while True:
            if len(bam_list) == 0:
                
                if picked.reference_end - end <= 100: 
                    break
                else:
                    end += int(qalen/20)
                    all_bam_list = list(inbam.fetch(ref_name, start, end))
                    bam_list = []
                    for r in all_bam_list:
                        if start <= r.pos <= end and r.reference_end > fe:
                            bam_list.append(r)
                    continue

            score_list = list(map(call_score, bam_list))
            max_id = score_list.index(max(score_list))
            picked = deepcopy(bam_list[max_id])
            if picked.flag == 16:
                picked.flag = 0
                picked.is_reverse = False
            ii += 1
            picked.qname = f'{ii}'
            
            newbam.write(picked)
            
            
            ALL.add(picked.qname)

            qalen = picked.query_alignment_length
            start = int(picked.pos + qalen*1/12)
            end = int(picked.pos + qalen*5/12)
            targets = int(picked.pos + qalen/4)
            targete = int(picked.reference_end + qalen/4)
            fe = int(picked.reference_end+qalen*0.05)
            call_score = lambda x: - pen*(abs(x.pos - targets) + abs(x.reference_end - targete)) - 100*x.get_tag('NM')
            all_bam_list = list(inbam.fetch(ref_name, start, end))
            bam_list = []
            for r in all_bam_list:
                if start <= r.pos <= end and r.reference_end > fe:
                    bam_list.append(r)
        inbam.close()
        newbam.close()
        num_rds.append(ii)
    print('\n')
    print('Completed,', num_rds[0], 'and', num_rds[1], 'reads have been filtered out.')
    print('')
    num_reads_file = out_dir + "/nr.txt"
    with open(num_reads_file,'w+') as f:
        f.write(str(num_rds[0])+' '+str(num_rds[1]))
