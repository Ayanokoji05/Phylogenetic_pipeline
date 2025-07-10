#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import time
from typing import List, Optional

# Always provide your email when using NCBI services
Entrez.email = "pratushkumar.pusti@niser.ac.in"
Entrez.api_key = None  # You can add an API key here if you have one
MAX_RETRIES = 3  # Maximum number of retries for API calls
DELAY = 1  # Delay between API calls in seconds

def search_nucleotide(query: str, max_results: int = 100) -> List[str]:
    """Search NCBI nucleotide database with error handling and retries.
    
    Args:
        query: Search term for NCBI nucleotide database
        max_results: Maximum number of results to return
        
    Returns:
        List of nucleotide accession IDs
    """
    for attempt in range(MAX_RETRIES):
        try:
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results)
            record = Entrez.read(handle)
            handle.close()
            return record["IdList"]
        except Exception as e:
            if attempt == MAX_RETRIES - 1:
                raise RuntimeError(f"Failed to search NCBI after {MAX_RETRIES} attempts: {str(e)}")
            time.sleep(DELAY * (attempt + 1))
    return []

def fetch_sequences(id_list: List[str], min_length: int = 500, max_length: int = 10000, 
                   batch_size: int = 100) -> List[SeqIO.SeqRecord]:
    """Fetch sequences and filter based on quality criteria with batch processing.
    
    Args:
        id_list: List of NCBI accession IDs
        min_length: Minimum sequence length to keep
        max_length: Maximum sequence length to keep
        batch_size: Number of sequences to fetch at once
        
    Returns:
        List of high-quality SeqRecord objects
    """
    high_quality_genes = []
    
    # Process IDs in batches to avoid overwhelming the server
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i + batch_size]
        for attempt in range(MAX_RETRIES):
            try:
                handle = Entrez.efetch(db="nucleotide", id=batch, rettype="gb", retmode="text")
                for record in SeqIO.parse(handle, "genbank"):
                    if (min_length <= len(record.seq) <= max_length and
                        has_complete_cds(record) and
                        not_putative(record) and
                        is_not_pseudogene(record)):
                        high_quality_genes.append(record)
                handle.close()
                break  # Success, exit retry loop
            except Exception as e:
                if attempt == MAX_RETRIES - 1:
                    print(f"Warning: Failed to fetch batch {i//batch_size + 1}: {str(e)}")
                    continue
                time.sleep(DELAY * (attempt + 1))
        time.sleep(DELAY)  # Be kind to NCBI servers
    
    return high_quality_genes

def has_complete_cds(record: SeqIO.SeqRecord) -> bool:
    """Check if record has a complete CDS feature.
    
    Args:
        record: SeqRecord object to check
        
    Returns:
        True if record has at least one complete CDS, False otherwise
    """
    for feature in record.features:
        if feature.type == "CDS":
            # Check for both partial flag and presence of start/stop codons
            if (not feature.location.is_partial and 
                "transl_start" in feature.qualifiers and
                "transl_table" in feature.qualifiers):
                return True
    return False

def not_putative(record: SeqIO.SeqRecord) -> bool:
    """Filter out putative, hypothetical, or unnamed proteins.
    
    Args:
        record: SeqRecord object to check
        
    Returns:
        True if the record doesn't contain putative/hypothetical proteins, False otherwise
    """
    for feature in record.features:
        if feature.type == "CDS":
            product = feature.qualifiers.get("product", [""])[0].lower()
            note = feature.qualifiers.get("note", [""])[0].lower()
            if ("putative" in product or 
                "hypothetical" in product or
                "unnamed" in product or
                "putative" in note or
                "hypothetical" in note):
                return False
    return True

def is_not_pseudogene(record: SeqIO.SeqRecord) -> bool:
    """Filter out pseudogenes.
    
    Args:
        record: SeqRecord object to check
        
    Returns:
        True if the record doesn't contain pseudogenes, False otherwise
    """
    for feature in record.features:
        if feature.type == "gene":
            gene_type = feature.qualifiers.get("gene_biotype", [""])[0].lower()
            if gene_type == "pseudogene":
                return False
    return True

if __name__ == "__main__":
    # Example usage
    try:
        query = "human[orgn] AND COX1[gene]"
        print(f"Searching for: {query}")
        ids = search_nucleotide(query, max_results=50)
        print(f"Found {len(ids)} results")
        
        if ids:
            sequences = fetch_sequences(ids)
            print(f"Retrieved {len(sequences)} high-quality sequences")
            
            # Save to file
            if sequences:
                with open("high_quality_sequences.fasta", "w") as f:
                    SeqIO.write(sequences, f, "fasta")
                print("Saved sequences to high_quality_sequences.fasta")
    except Exception as e:
        print(f"Error: {str(e)}")