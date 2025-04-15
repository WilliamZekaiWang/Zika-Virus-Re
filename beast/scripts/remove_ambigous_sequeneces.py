import csv

dic = {}
    
ambig = {
    "R",  # A or G
    "Y",  # C or T
    "S",  # G or C
    "W",  # A or T
    "K",  # G or T
    "M",  # A or C
    "B",  # C, G, or T
    "D",  # A, G, or T
    "H",  # A, C, or T
    "V",  # A, C, or G
    "N"   # Any base (A, T, C, G)
}

# might have broken the code by trying to fix the problem that the dates are dropped because of the way I split the seq title
with open('zika_dates.fa', 'r') as file:
    for row in file:
        row = row.strip()
        if row.startswith('>'):
            seq_name = row[1:].split()[0]  # Extract just the accession number
            seq_full = row[1:]
            dic[seq_full] = []
        elif seq_full:
            count = 0
            for char in row:
                if char.upper() in ambig:
                    count += 1
                if count > 10:
                    dic.pop(seq_full)
                    seq_name = None
                    seq_full = None
                    break
            if seq_name:
                dic[seq_full].append(row)

with open('zika_clean_dates.fa', 'w') as file:
    for name, seq_lines in dic.items():
        file.write(f'>{name}\n')
        file.write('\n'.join(seq_lines) + '\n')
