#!/usr/bin/env python3

#!/usr/bin/env python3

#!/usr/bin/env python3

import os
import json
import subprocess
from collections import defaultdict
from Bio import Phylo

TREE_FILE = "annotated_tree.nwk"
METADATA_DIR = "metadatas"
OUTPUT_LATEX = "metadata_report.tex"

def read_json_metadata(filepath):
    with open(filepath, "r") as f:
        return json.load(f)

def collect_metadata(tree_file, metadata_dir):
    tree = Phylo.read(tree_file, "newick")
    grouped = defaultdict(list)

    for leaf in tree.get_terminals():
        accession = leaf.name.split('|')[0].strip()
        metadata_path = os.path.join(metadata_dir, f"{accession}.json")

        if not os.path.exists(metadata_path):
            print(f"⚠️ Missing metadata for {accession}")
            continue

        try:
            data = read_json_metadata(metadata_path)
        except Exception as e:
            print(f"❌ Error reading {accession}: {e}")
            continue

        gene = data.get("Gene", "Unknown")
        grouped[gene].append(data)

    return grouped

def safe_latex(s):
    return str(s).replace('\\', r'\\').replace('_', r'\_').replace('%', r'\%')\
        .replace('&', r'\&').replace('#', r'\#').replace('$', r'\$')\
        .replace('{', r'\{').replace('}', r'\}').replace('^', r'\^{}')\
        .replace('~', r'\~{}')

def write_latex(grouped_data, output_file):
    with open(output_file, "w") as f:
        # Preamble with standard packages
        f.write(r"\documentclass[10pt]{article}" + "\n")
        f.write(r"\usepackage[T1]{fontenc}" + "\n")
        f.write(r"\usepackage{lmodern}" + "\n")
        f.write(r"\usepackage{longtable}" + "\n")
        f.write(r"\usepackage[colorlinks=true,linkcolor=blue]{hyperref}" + "\n")
        f.write(r"\usepackage{array}" + "\n")
        f.write(r"\usepackage{tocloft}" + "\n")
        f.write(r"\usepackage{microtype}" + "\n")
        f.write(r"\usepackage[margin=1in]{geometry}" + "\n")
        f.write(r"\renewcommand{\cftsecleader}{\cftdotfill{\cftdotsep}}" + "\n")
        
        # Document metadata
        f.write(r"\title{Metadata Report}" + "\n")
        f.write(r"\date{\today}" + "\n")
        f.write(r"\begin{document}" + "\n")
        f.write(r"\maketitle" + "\n")
        
        # Table of Contents (automatically populated by sections)
        f.write(r"\tableofcontents" + "\n")
        f.write(r"\clearpage" + "\n")

        # Metadata sections with automatic TOC entries
        for gene, entries in sorted(grouped_data.items()):
            safe_gene = safe_latex(gene)
            f.write(fr"\section{{Gene: {safe_gene}}}" + "\n")
            
            for entry in entries:
                f.write(r"{\footnotesize" + "\n")
                f.write(r"\begin{longtable}{>{\raggedright\arraybackslash}p{4.5cm} >{\raggedright\arraybackslash}p{11.5cm}}" + "\n")
                f.write(r"\textbf{Field} & \textbf{Value} \\" + "\n")
                f.write(r"\hline" + "\n")

                # Standard fields
                for key in ["Accession", "Organism", "Length", "Molecule Type", "Description"]:
                    val = entry.get(key, "Unknown")
                    f.write(f"{safe_latex(key)} & {safe_latex(val)} \\\\" + "\n")

                # Optional fields
                if entry.get("Keywords"):
                    kws = ", ".join(entry["Keywords"])
                    f.write(f"Keywords & {safe_latex(kws)} \\\\" + "\n")

                if entry.get("Taxonomy"):
                    tax = ", ".join(entry["Taxonomy"])
                    f.write(f"Taxonomy & {safe_latex(tax)} \\\\" + "\n")

                f.write(r"\end{longtable}" + "\n")
                f.write(r"}" + "\n\n")
                f.write(r"\vspace{1em}" + "\n")

        f.write(r"\end{document}" + "\n")
    print(f"✅ LaTeX report written to {output_file}")

def compile_pdf(latex_file):
    try:
        print(f"🛠️ Compiling LaTeX with pdflatex...")
        # First compilation (generates TOC)
        subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", latex_file], 
            capture_output=True
        )
        # Second compilation (resolves references)
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", latex_file], 
            capture_output=True, 
            text=True
        )
        
        if result.returncode != 0:
            print("❌ LaTeX compilation failed")
            print(result.stderr)
            return False
        
        return True
        
    except FileNotFoundError:
        print("❌ pdflatex not found. Please install LaTeX:")
        print("   sudo apt-get install texlive-latex-base texlive-latex-recommended")
        return False

def main():
    print("🚀 Starting metadata report generation...")
    
    if not os.path.exists(TREE_FILE):
        print(f"❌ Tree file not found: {TREE_FILE}")
        return
    
    if not os.path.exists(METADATA_DIR):
        print(f"❌ Metadata directory not found: {METADATA_DIR}")
        return
    
    grouped_metadata = collect_metadata(TREE_FILE, METADATA_DIR)
    
    if not grouped_metadata:
        print("❌ No metadata found. Check your tree file and metadata directory.")
        return
    
    print(f"📊 Found {len(grouped_metadata)} gene groups with {sum(len(entries) for entries in grouped_metadata.values())} total entries")
    
    write_latex(grouped_metadata, OUTPUT_LATEX)
    success = compile_pdf(OUTPUT_LATEX)
    
    if success:
        print("✅ Report generation completed successfully!")
    else:
        print("❌ Report generation failed during PDF compilation.")

if __name__ == "__main__":
    main()