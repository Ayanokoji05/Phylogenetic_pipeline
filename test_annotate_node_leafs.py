#!/usr/bin/env python3

from Bio import Entrez, Phylo, SeqIO
import os
import time
from io import StringIO

Entrez.email = "pratushkumar.pusti@niser.ac.in"

def detect_id_type(accession):
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
            return "nucleotide"  # default fallback

def get_gene_info_from_id(accession):
    accession = accession.strip()
    db = detect_id_type(accession)

    print(f"🔍 Fetching {accession} from {db}")

    try:
        handle = Entrez.efetch(db=db, id=accession, rettype="gb", retmode="text", timeout=10)
        response = handle.read()
        handle.close()

        if not response.strip():
            print(f"⚠️ Empty response for {accession}")
            return {'organism': 'Unknown', 'gene_symbol': accession}

        record = SeqIO.read(StringIO(response), "genbank")
    except Exception as e:
        print(f"❌ Error fetching {accession}: {e}")
        return {'organism': 'Unknown', 'gene_symbol': accession}

    organism = record.annotations.get("organism", "Unknown")
    gene_symbol = None
    for feature in record.features:
        if feature.type == "CDS" and "gene" in feature.qualifiers:
            gene_symbol = feature.qualifiers["gene"][0]
            break

    if not gene_symbol:
        print(f"⚠️ No gene found in {accession}, checking linked...")
        try:
            handle = Entrez.elink(dbfrom=db, db="nuccore" if db == "protein" else "protein", id=accession, timeout=10)
            link_results = Entrez.read(handle)
            handle.close()

            if link_results and link_results[0].get("LinkSetDb"):
                links = link_results[0]["LinkSetDb"][0].get("Link", [])
                if links:
                    linked_id = links[0]["Id"]
                    handle = Entrez.efetch(db="nuccore", id=linked_id, rettype="gb", retmode="text", timeout=10)
                    linked_response = handle.read()
                    handle.close()

                    if not linked_response.strip():
                        print(f"⚠️ Empty linked record for {accession}")
                        return {'organism': organism, 'gene_symbol': accession}

                    linked_record = SeqIO.read(StringIO(linked_response), "genbank")

                    for feature in linked_record.features:
                        if feature.type == "CDS" and "gene" in feature.qualifiers:
                            gene_symbol = feature.qualifiers["gene"][0]
                            break
        except Exception as e:
            print(f"⚠️ Linked record fetch failed: {e}")

    return {
        'organism': organism,
        'gene_symbol': gene_symbol if gene_symbol else accession
    }


def annotate_tree_with_gene_info(tree_file, output_file):
    tree = Phylo.read(tree_file, 'newick')

    for i, leaf in enumerate(tree.get_terminals()):
        if not leaf.name:
            print(f"❌ [{i}] Leaf has no name — skipping.")
            continue

        print(f"▶️ [{i}] Processing {leaf.name}")
        accession = leaf.name.split("|")[0]  # Extract accession without trailing info
        gene_info = get_gene_info_from_id(accession)

        # ✅ Include gene symbol in label
        new_name = f"{accession}|{gene_info['gene_symbol']}|{gene_info['organism']}"
        leaf.name = new_name


# Example usage
if __name__ == "__main__":
    annotate_tree_with_gene_info("raxml_output.raxml.support", "annotated_tree.nwk")