#%%
import subprocess
from filterbam import filter
from reads_helper import get_reads
import os
import time
import argparse
parser = argparse.ArgumentParser(description='pass commands')
parser.add_argument("-vn", type=int, help="dir name saving filtered reads")
args = parser.parse_args()
vn = args.vn

st = time.time()
dir = "/mnt/d/projects/vrpassembler"
# dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


bam_dir = dir + f'/data/yeast/BFH/sim_reads'
print('saving to', bam_dir)
penf = 500

filter(penf, bam_dir)

rn = []
for i in [1,2]:
    rf = bam_dir + f'/r{i}/readsh{i}.bam'
    # rf = bam_dir + f'final{i}.bam'
    outf = bam_dir + f'/r{i}/readsh{i}.sort.bam'
    process = subprocess.Popen(f"samtools sort {rf} > {outf}",
                                shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    process.wait()

    process = subprocess.Popen(f"samtools index {outf}",
                                shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    process.wait()

    irn = bam_dir + f'/r{i}/reads_{i}.fa'
    rn.append(irn)
    process = subprocess.Popen(f"samtools fasta -0 {irn} {outf}",
                                shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    process.wait()
    command_output = process.stdout.read().decode('utf-8')
    print(command_output)

# cat two reads
all_r = bam_dir + '/sim_reads_tmp.fa'
process = subprocess.Popen(f"cat {rn[0]} {rn[1]} > {all_r}",
                                shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
process.wait()

get_reads(bam_dir)
print('ALL COMPLETED! Time cost:',round(time.time()-st,3), 'sec.')


