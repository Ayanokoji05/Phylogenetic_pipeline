#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import os

# Set your email for NCBI Entrez API access
Entrez.email = "your_email@example.com"

# Path setup (relative to script location)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
CDS_DIR = os.path.join(DATA_DIR, "cds_filtered")
GENE_LIST_FILE = os.path.join(SCRIPT_DIR, "genes.list")

MIN_CDS_LENGTH = 300  # Minimum CDS length for QC

def read_gene_list(filename=GENE_LIST_FILE):
    with open(filename) as f:
        return [line.strip() for line in f if line.strip()]

def search_and_fetch(gene, species_limit=500):
    query = f"{gene}[Gene]"
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=species_limit)
    record = Entrez.read(handle)
    return record['IdList']

def has_sequence(seq_id):
    """Check if the GenBank record has sequence content using XML mode."""
    try:
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="xml")
        records = Entrez.read(handle)
        if records and 'GBSeq_sequence' in records[0]:
            seq = records[0]['GBSeq_sequence']
            return seq and len(seq) > 0
    except Exception as e:
        print(f"[WARN] XML check failed for {seq_id}: {e}")
    return False

def extract_cds_from_ids(id_list, gene):
    good_cds = []
    for seq_id in id_list:
        if not has_sequence(seq_id):
            print(f"[SKIP] {seq_id} has no sequence content")
            continue
        try:
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            for feature in record.features:
                if feature.type == "CDS" and "gene" in feature.qualifiers:
                    if gene.lower() in [x.lower() for x in feature.qualifiers["gene"]]:
                        seq = feature.extract(record.seq)
                        if len(seq) >= MIN_CDS_LENGTH and seq.startswith("ATG") and len(seq) % 3 == 0:
                            prot = seq.translate(to_stop=True)
                            if "*" not in prot:
                                good_cds.append((record.id, seq))
        except Exception as e:
            print(f"[ERROR] Failed for {seq_id}: {e}")
    return good_cds

def save_cds(gene, cds_list):
    os.makedirs(CDS_DIR, exist_ok=True)
    out_file = os.path.join(CDS_DIR, f"{gene}_cds.fa")
    with open(out_file, "w") as f:
        for sid, seq in cds_list:
            f.write(f">{sid}|{gene}\n{seq}\n")

def main():
    genes = read_gene_list()
    for gene in genes:
        print(f"🔍 Processing gene: {gene}")
        ids = search_and_fetch(gene)
        cds = extract_cds_from_ids(ids, gene)
        print(f"✅ {len(cds)} valid CDS for {gene}")
        if cds:
            save_cds(gene, cds)

if __name__ == "__main__":
    main()