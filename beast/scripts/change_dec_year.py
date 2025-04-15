from Bio import SeqIO
from datetime import datetime
import sys

def date_to_decimal(date_str):
    date = datetime.strptime(date_str, "%Y-%m-%d")
    year = date.year
    start_of_year = datetime(year, 1, 1)
    end_of_year = datetime(year + 1, 1, 1)
    decimal_year = year + (date - start_of_year).total_seconds() / (end_of_year - start_of_year).total_seconds()
    decimal_year -= 1947  # earliest date
    return f"{decimal_year:.6f}"

def process_fasta(input_fasta, output_fasta):
    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            parts = record.id.rsplit("_", 1)
            if len(parts) == 2:
                decimal_year = date_to_decimal(parts[1])
                if decimal_year:
                    record.id = f"{parts[0]}_{decimal_year}"
                    record.description = ""
            SeqIO.write(record, out_handle, "fasta")

# credit to my friend Ray Liu at University of Waterloo for helping me make this
if __name__ == "__main__":
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    
    process_fasta(input_fasta, output_fasta)
