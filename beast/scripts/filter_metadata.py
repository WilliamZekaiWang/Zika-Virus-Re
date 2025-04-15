import csv

# get headers to keep
keep = []
with open('zika_aln.fa', 'r') as file:
    for row in file:
        if row[0] == '>':
            keep.append(row[1:].strip().split('_')[0]) # because I added decimal years lol  

save = []
with open('zika_metadata.csv', 'r') as file:
    file = csv.DictReader(file)
    for row in file:
        if row['header'] in keep:    
            save.append(row)

print(len(keep))
print(len(save))

# write
with open('zika_metadata_final.csv', 'w', newline='') as file:
    fieldnames = save[0].keys() if save else []  # Get the fieldnames from the first row
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    
    writer.writeheader() 
    writer.writerows(save) 
