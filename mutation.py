import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import Draw
import random
import io
import numpy as np
import requests

# Streamlit Page Title
st.title("BIOPHYLO: Phylogenetic Tree & Molecular Docking Tool")

# Sequence Input Section
st.header("1️⃣ Enter Sequence Alignment (FASTA Format)")
seq_input = st.text_area("Paste your sequence alignment here (FASTA format)")

# Protein & Ligand Input Sections
st.header("2️⃣ Enter Protein List")
protein_input = st.text_area("Enter protein names (one per line)")
proteins = [p.strip() for p in protein_input.split('\n') if p.strip()]

st.header("3️⃣ Enter Ligand List")
ligand_input = st.text_area("Enter ligand names (one per line)")
ligands = [l.strip() for l in ligand_input.split('\n') if l.strip()]


# Function to fetch SMILES notation from PubChem
@st.cache_data
def get_smiles(ligand_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ligand_name}/property/IsomericSMILES/TXT"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    return None


# Auto-generate Molecular Docking Results
def generate_docking_results(proteins, ligands):
    if not proteins or not ligands:
        return None
    data = {
        "Protein_ID": random.choices(proteins, k=10),
        "Ligand_ID": random.choices(ligands, k=10),
        "Docking_Score": [random.uniform(-10, -5) for _ in range(10)]
    }
    df = pd.DataFrame(data)
    return df


# Function to generate phylogenetic tree from sequence alignment
def generate_phylogenetic_tree(seq_input, docking_df):
    if seq_input:
        alignment = AlignIO.read(io.StringIO(seq_input), "fasta")
        calculator = DistanceCalculator("identity")
        dm = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)

        fig, ax = plt.subplots()
        Phylo.draw(tree, axes=ax)

        # Assign colors based on docking scores
        color_map = {}
        if docking_df is not None:
            min_score = min(docking_df["Docking_Score"])
            max_score = max(docking_df["Docking_Score"])
            norm = plt.Normalize(vmin=min_score, vmax=max_score)
            cmap = plt.cm.coolwarm

            for _, row in docking_df.iterrows():
                score = row["Docking_Score"]
                color_map[row["Protein_ID"]] = cmap(norm(score))

            best_score = min(docking_df["Docking_Score"])
            best_entry = docking_df.loc[docking_df["Docking_Score"] == best_score]
            st.write(f"### Best Docking Score: {best_score:.2f}")
            st.dataframe(best_entry)

        # Color nodes in the tree based on docking scores
        for clade in tree.get_terminals():
            if clade.name in color_map:
                clade.color = color_map[clade.name]

        st.pyplot(fig)
    else:
        st.warning("Please enter a valid sequence alignment in FASTA format.")


# Function to visualize ligand structures
def visualize_ligand(ligand_name):
    smiles = get_smiles(ligand_name)
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol)
            st.image(img, caption=f"Ligand: {ligand_name} ({smiles})")
        else:
            st.warning(f"Invalid ligand structure for {ligand_name}")
    else:
        st.warning(f"SMILES not found for {ligand_name}")


# Main execution logic
if proteins and ligands:
    docking_df = generate_docking_results(proteins, ligands)
    if docking_df is not None:
        st.write("### Docking Scores:")
        st.dataframe(docking_df)

        # Generate phylogenetic tree with docking results only if sequence input is provided
        if seq_input:
            generate_phylogenetic_tree(seq_input, docking_df)
        else:
            st.warning("Please enter a sequence alignment to generate the phylogenetic tree.")

        st.write("### Ligand Structures:")
        for ligand in ligands:
            visualize_ligand(ligand)
    else:
        st.warning("Please enter both protein and ligand lists to generate docking results.")
elif seq_input:
    # Only generate phylogenetic tree without docking results if no proteins/ligands provided
    generate_phylogenetic_tree(seq_input, None)