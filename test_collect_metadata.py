#!/usr/bin/env python3

import os
import json
import time
from io import StringIO, BytesIO
from Bio import Entrez, Phylo, SeqIO

Entrez.email = "pratushkumar.pusti@niser.ac.in"

TREE_FILE = "annotated_tree.nwk"
METADATA_DIR = "metadatas"

def detect_id_type(accession):
    """Try to determine whether the accession belongs to nucleotide or protein DB."""
    try:
        handle = Entrez.esummary(db="nucleotide", id=accession)
        Entrez.read(handle)
        handle.close()
        return "nucleotide"
    except:
        try:
            handle = Entrez.esummary(db="protein", id=accession)
            Entrez.read(handle)
            handle.close()
            return "protein"
        except:
            return None  # Unknown accession type

def fetch_metadata(accession):
    db = detect_id_type(accession)
    if not db:
        print(f"⚠️ Unknown accession type for {accession}")
        return None

    print(f"🔍 Fetching metadata for {accession} from {db}")
    
    # First try gbwithparts (text)
    try:
        handle = Entrez.efetch(db=db, id=accession, rettype="gbwithparts", retmode="text", timeout=20)
        response = handle.read()
        handle.close()
        record = SeqIO.read(StringIO(response), "genbank")
    except Exception as e1:
        print(f"⚠️ gbwithparts failed for {accession}: {e1}")
        # Try XML as fallback
        try:
            print(f"⚠️ Retrying {accession} with XML mode...")
            handle = Entrez.efetch(db=db, id=accession, rettype="gb", retmode="xml", timeout=20)
            response = handle.read()
            handle.close()
            record = SeqIO.read(BytesIO(response), "genbank")
        except Exception as e2:
            print(f"❌ Both fetch methods failed for {accession}: {e2}")
            return None

    # Extract gene symbol
    gene_symbol = None
    for feature in record.features:
        if feature.type == "CDS" and "gene" in feature.qualifiers:
            gene_symbol = feature.qualifiers["gene"][0]
            break

    metadata = {
        "Accession": accession,
        "Organism": record.annotations.get("organism", "Unknown"),
        "Gene": gene_symbol if gene_symbol else "Unknown",
        "Length": len(record.seq),
        "Molecule Type": record.annotations.get("molecule_type", "Unknown"),
        "Description": record.description,
        "Keywords": record.annotations.get("keywords", []),
        "References": [ref.title for ref in record.annotations.get("references", [])],
        "Taxonomy": record.annotations.get("taxonomy", []),
        "DB Source": record.annotations.get("db_source", "Unknown"),
        "Version": record.annotations.get("sequence_version", "Unknown"),
        "Sequence": str(record.seq)
    }

    return metadata

def save_metadata(accession, metadata, outdir):
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, f"{accession}.json")
    with open(out_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"✅ Saved metadata to: {out_path}")

def extract_unique_accessions(tree_file):
    tree = Phylo.read(tree_file, "newick")
    accessions = set()
    for leaf in tree.get_terminals():
        label_parts = leaf.name.split("|")
        accession = label_parts[0].strip()
        accessions.add(accession)
    return sorted(accessions)

def main():
    print(f"🌳 Reading tree from: {TREE_FILE}")
    accessions = extract_unique_accessions(TREE_FILE)
    print(f"📦 Found {len(accessions)} unique accessions")

    for accession in accessions:
        metadata_file = os.path.join(METADATA_DIR, f"{accession}.json")
        if os.path.exists(metadata_file):
            print(f"📁 Metadata already exists for {accession} — skipping.")
            continue

        metadata = fetch_metadata(accession)
        if metadata:
            save_metadata(accession, metadata, METADATA_DIR)
        time.sleep(0.5)  # Respect NCBI rate limits

if __name__ == "__main__":
    main()