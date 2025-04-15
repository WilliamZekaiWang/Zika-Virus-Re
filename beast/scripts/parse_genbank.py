from Bio import SeqIO
import argparse
import csv
import sys

# command line interface
description = """
Extract metadata from a Genbank file and write the results 
to a CSV file.  Write sequences to a separate FASTA file.
"""
parser = argparse.ArgumentParser(description=description)
parser.add_argument("infile", type=argparse.FileType('r'),
                    help="Path to Genbank file to parse.")
parser.add_argument("csvfile", type=argparse.FileType('w'),
                    help="Path to write CSV file with metadata.")
parser.add_argument("seqfile", type=argparse.FileType('w'),
                    help="Path to write FASTA file with sequences.")
args = parser.parse_args()

# prepare output
fields = ['accession', 'Sample City', 'Risk factor', 'Viral load', 
          'CD4 count', 'sample tissue', 'country', 'coldate']
writer = csv.DictWriter(args.csvfile, fieldnames=fields)
writer.writeheader()

records = SeqIO.parse(args.infile, "genbank")
for record in records:
    # extract metadata from comment field
    try:
        mdata = record.annotations['structured_comment']['HIVDataBaseData']
        mdata['accession'] = record.id
        record.description = ""

        # extract source metadata
        source = [feat for feat in record.features if feat.type=='source']
        if source:
            quals = source[0].qualifiers
            mdata['country'] = quals.get('geo_loc_name', [None])[0]
            mdata['coldate'] = quals.get('collection_date', [None])[0]

        writer.writerow(mdata)  # write metadata
        SeqIO.write(record, args.seqfile, 'fasta')  # write sequence
    except:
        continue
