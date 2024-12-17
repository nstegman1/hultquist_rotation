import os
import glob 

base_path = '/projects/b1042/HultquistLab/Daphne/Novogene_NVUS2024071006/Raw_Data/'
zipped_path = '/projects/p32549/RNA_seq_analysis_ht1080/ht1080_raw_fastqs/'
general_out = '/projects/p32549/RNA_seq_analysis_ht1080/'

# Split R1 files into R1_ht1080.txt and R2_ht1080.txt (spearate, sort, and write out filenames.)

zipped_filenames = os.listdir(zipped_path)
r1_filenames = []
r2_filenames = []

for i in zipped_filenames:
    splitted = i.split('_')
    read_num = splitted[2]
    print(read_num)
    if read_num[0] == '1':
        print('I am R1 ',i)
        r1_filenames.append(i)
    else:
        print('I am R2',i)
        r2_filenames.append(i)

r1_filenames = sorted(r1_filenames)
r2_filenames = sorted(r2_filenames)

print('are they the same length?')
print(len(r1_filenames), len(r2_filenames))

print('do they match?:')
r1_out = open(general_out+'R1_ht1080.txt','w')
r2_out = open(general_out+'R2_ht1080.txt','w')
for i in range(0, len(r1_filenames)):
    print(r1_filenames[i], r2_filenames[i])
    r1_out.write(r1_filenames[i]+'\n')
    r2_out.write(r2_filenames[i]+'\n')

r1_out.close()
r2_out.close()