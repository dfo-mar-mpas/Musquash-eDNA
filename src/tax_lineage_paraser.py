
import pandas as pd
import sys

def load_ncbi_taxonomy(nodes_path, names_path):
    print("Loading NCBI taxonomy into memory... (1-2 minutes)")
    
    # taxid -> (parent_id, rank)
    nodes = {}
    with open(nodes_path, 'r') as f:
        for line in f:
            parts = line.split('|')
            tax_id = parts[0].strip()
            parent_id = parts[1].strip()
            rank = parts[2].strip()
            nodes[tax_id] = (parent_id, rank)

    # taxid -> scientific_name
    names = {}
    with open(names_path, 'r') as f:
        for line in f:
            if 'scientific name' in line:
                parts = line.split('|')
                tax_id = parts[0].strip()
                name = parts[1].strip()
                names[tax_id] = name
                
    return nodes, names

def get_lineage_string(tax_id, nodes, names):
    if pd.isna(tax_id):
        return ""
    
    # Convert tax_id to string to match dictionary keys
    current_id = str(int(tax_id))
    
    # NCBI uses 'superkingdom' for Bacteria/Archaea
    target_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    lineage = {rank: "" for rank in target_ranks}
    
    # Traverse up the tree to the root
    while current_id != '1' and current_id in nodes:
        parent_id, rank = nodes[current_id]
        if rank in lineage:
            lineage[rank] = names.get(current_id, f"unknown_{current_id}")
        current_id = parent_id
        
    # Join with commas in the specific order requested
    return ",".join([lineage[r] for r in target_ranks])

def process_csv(input_csv, output_csv, nodes_path, names_path, taxid_column='taxid'):
    # Load Taxonomy
    nodes, names = load_ncbi_taxonomy(nodes_path, names_path)
    
    # Read CSV
    print(f"Reading {input_csv}...")
    df = pd.read_csv(input_csv)
    
    if taxid_column not in df.columns:
        print(f"Error: Column '{taxid_column}' not found. Available columns: {list(df.columns)}")
        return

    # Apply lineage mapping
    print("Mapping lineages...")
    df['taxlineage'] = df[taxid_column].apply(lambda x: get_lineage_string(x, nodes, names))
    
    # Save to new CSV
    df.to_csv(output_csv, index=False)
    print(f"Success! File saved as: {output_csv}")

# --- Configuration ---
WORKING_DIR = '/mnt/e/projects/Blast_visualizer/euk/'
INPUT_FILE = f'{WORKING_DIR}/cosmic_euk_20260114_213237_coordinates.csv'
OUTPUT_FILE = f'{WORKING_DIR}/euk_annotated_lineages.csv'
NAMES_DMP_PATH = '/mnt/e/projects/databases/tax_dump/ncbi_taxdump'
NODES_DMP = f'{NAMES_DMP_PATH}/nodes.dmp'
NAMES_DMP = f'{NAMES_DMP_PATH}/names.dmp'
TAXID_COL = 'taxonomy_id'  # Change this to match your CSV's column name

if __name__ == "__main__":
    process_csv(INPUT_FILE, OUTPUT_FILE, NODES_DMP, NAMES_DMP, TAXID_COL)