"""
script to drop sepcific sequences in fasta file
"""
dic = {}

to_drop = {
'KX198134.2',
'MT505350.1',
'MT505349.1',
'KY766069.1'
}


seq_name = None

with open('zika_downsample.fa', 'r') as file:
    for row in file:
        row = row.strip()
        if row.startswith('>'):
            seq_name = row[1:].split('_')[0]
            full = row[1:]
            if seq_name in to_drop:
                seq_name = None
                full = None
            else:
                dic[full] = []
        elif seq_name:
            dic[full].append(row)

with open('zika_no_out.fa', 'w') as file:
    for name, seq_lines in dic.items():
        file.write(f'>{name}\n')
        file.write('\n'.join(seq_lines) + '\n')

