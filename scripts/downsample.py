import random

random.seed(100)
dic = {}

with open('zika_clean_dates.fa', 'r') as file:
    for row in file:
        row = row.strip()
        if row.startswith('>'):
            seq_name, date = row[1:].split('_')
            year = date.split('-')[0]
            dic.setdefault(year, []).append(seq_name)

# please be evenly distributed pleasepleaseplease
for year, names in dic.items():
    print(f'{year}: {len(names)}')

# Randomly sample sequences per year
save = []
for year, names in dic.items():
    save.extend(random.sample(names, min(15, len(names))))

# Second pass: Collect selected sequences
final = {}

with open('zika_clean_dates.fa', 'r') as file:
    seq_name = None
    seq_full = None 
    for row in file:
        row = row.strip()
        if row.startswith('>'):
            seq_name = row[1:].split('_')[0]
            seq_full = row[1:]
            if seq_name in save:
                final[seq_full] = ""
        elif seq_full in final:
            final[seq_full] += row

# Write downsampled sequences
with open('zika_downsample.fa', 'w') as file:
    for name, seq in final.items():
        file.write(f'>{name}\n{seq}\n')

