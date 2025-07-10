#!/usr/bin/env python3

import os
from Bio import Phylo
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import subprocess

def tree_to_forest(tree):
    """Convert a Bio.Phylo tree to LaTeX forest format recursively."""
    def _to_forest(clade):
        if clade.is_terminal():
            node_str = f"[{clade.name}]"
        else:
            children = [_to_forest(child) for child in clade.clades]
            node_str = f"[{clade.name or ''} {' '.join(children)}]"
        return node_str
    
    return _to_forest(tree.root)

def create_tree_pdf(input_tree="annotated_tree.nwk", output_pdf="tree_figure.pdf"):
    """
    Create a PDF of a phylogenetic tree using LaTeX forest rendering.
    
    Args:
        input_tree (str): Path to Newick format tree file (.nwk)
        output_pdf (str): Path for output PDF file
    """
    # 1. Read the Newick tree file
    tree = Phylo.read(input_tree, "newick")
    
    # 2. Generate LaTeX forest code
    forest_code = tree_to_forest(tree)
    
    # 3. Create a temporary LaTeX file
    latex_content = r"""
\documentclass{article}
\usepackage[paperwidth=17in,paperheight=11in,margin=0in]{geometry}
\usepackage{graphicx}
\usepackage{tikz}
\usepackage{forest}

\begin{document}
\thispagestyle{empty}
\begin{figure}[p]
    \centering
    \begin{forest}
        for tree={
            grow'=0,
            child anchor=west,
            parent anchor=east,
            anchor=west,
            calign=first,
            edge path={
                \noexpand\path [draw, \forestoption{edge}]
                (!u.south west) +(7.5pt,0) |- (.child anchor)\forestoption{edge label};
            },
            font=\sffamily
        }
""" + forest_code + r"""
    \end{forest}
    \caption{Phylogenetic Tree Visualization}
    \label{fig:tree}
\end{figure}
\end{document}
"""
    
    # 4. Save LaTeX to temporary file
    temp_tex = "temp_tree.tex"
    with open(temp_tex, "w") as f:
        f.write(latex_content)
    
    # 5. Compile LaTeX to PDF
    try:
        subprocess.run(["pdflatex", "-interaction=nonstopmode", temp_tex], check=True)
        
        # 6. Rename the output
        if os.path.exists("temp_tree.pdf"):
            os.rename("temp_tree.pdf", output_pdf)
            print(f"Successfully generated {output_pdf}")
        else:
            print("Error: PDF generation failed")
    except subprocess.CalledProcessError as e:
        print(f"LaTeX compilation failed: {e}")
    except FileNotFoundError:
        print("Error: pdflatex not found. Please install a LaTeX distribution.")
    
    # Clean up auxiliary files
    for ext in [".aux", ".log"]:
        if os.path.exists("temp_tree" + ext):
            os.remove("temp_tree" + ext)

if __name__ == "__main__":
    create_tree_pdf()