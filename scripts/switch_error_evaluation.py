#%%
import pysam
import os
import numpy as np
#%%
cur_file = os.path.abspath(__file__)
cur_dir = os.path.dirname(cur_file)

rl = 30

for i in range(1,11):
    for j in range(1,3):
        vcf_all_path = cur_dir + '/H_all.vcf.gz'
        vcf_true_all = pysam.VariantFile(vcf_all_path)
        vcf_true1_path = cur_dir + '/H1.vcf.gz'
        vcf_true1 = pysam.VariantFile(vcf_true1_path)
        vcf_true2_path = cur_dir + '/H2.vcf.gz'
        vcf_true2 = pysam.VariantFile(vcf_true2_path)
        # vcf_phased_path = cur_dir + f'/data/assemblies/hifiasm/{rl}k/{i}/hap{j}.phase.vcf.gz'
        vcf_phased_path = cur_dir + f'/data/assemblies/VRP_assembler/{rl}k/{i}/hap{j}.racon.phase.vcf.gz'
        vcf_phased = pysam.VariantFile(vcf_phased_path)
        cnt = 0
        error = 0
        hap_true_phase = []
        for record in vcf_true_all:
            
            chromosome = record.chrom
            position = record.pos
            ref_true = record.ref
            alt_true = record.alts

            recs1 = vcf_true1.fetch(chromosome, position-1, position)
            ref_true_1 = '-'
            alt_true_1 = ('-',)
            for irec in recs1:
                ref_true_1 = irec.ref
                alt_true_1 = irec.alts


            recs2 = vcf_true2.fetch(chromosome, position-1, position)
            ref_true_2 = '-'
            alt_true_2 = ('-',)
            for irec in recs2:
                ref_true_2 = irec.ref
                alt_true_2 = irec.alts

            if alt_true_1[0]==alt_true_2[0]:
                    continue

            phase_recs = vcf_phased.fetch(chromosome, position-1, position)
            ref_phase = '-'
            alt_phase = ('-',)
            for irec in phase_recs:
                ref_phase = irec.ref
                alt_phase = irec.alts
            cnt += 1

            if alt_phase[0] == alt_true_1[0]:
                hap_true_phase.append(0)
            elif alt_phase[0] == alt_true_2[0]:
                hap_true_phase.append(1)
            else:
                error += 1
                cnt -= 1

        hap_true_phase_arr = np.array(hap_true_phase)
        tmp = hap_true_phase_arr[1:]-hap_true_phase_arr[:-1]
        res = tmp==1
        se = np.count_nonzero(res)
        print(f'v{rl}{i} hap{j}: {round(se*100/cnt,5)}%')

# %%
