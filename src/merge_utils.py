# #!/usr/bin/env python3
# """
# Utilities for merging taxonomies and filtering based on expected lists.
# """

import csv
import os

def load_taxonomy_dict(file_path):
    """
    Load a QIIME taxonomy file into a dictionary.
    Returns: {feature_id: {'taxon': str, 'conf': str}}
    """
    tax_dict = {}
    if not os.path.exists(file_path):
        return tax_dict
        
    with open(file_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None) # Skip header
        for row in reader:
            if len(row) >= 2:
                conf = row[2] if len(row) > 2 else "0.0"
                tax_dict[row[0]] = {'taxon': row[1], 'conf': conf}
    return tax_dict

def clean_rank_name(rank_part):
    """
    Extract scientific name from QIIME rank string.
    Handles 's__Homo sapiens (human)' -> 'Homo sapiens'
    """
    if "__" in rank_part:
        name = rank_part.split("__", 1)[1]
    else:
        name = rank_part
        
    if " (" in name and name.endswith(")"):
        name = name.rpartition(" (")[0]
        
    return name.strip()

def _lineage_contains_expected(lineage_string, target_names):
    """
    Checks if any name in the QIIME lineage string exists in the target_names set.
    """
    if not lineage_string or lineage_string == "Unclassified":
        return False
    
    ranks = lineage_string.split(';')
    for rank in ranks:
        if "s__" in rank:
            clean_name = clean_rank_name(rank)
            if clean_name and clean_name in target_names:
                return True
    return False

def merge_and_filter_taxonomy(alpha_file, beta_file, lookup_file, output_file):
    """
    Merges two taxonomies (Alpha and Beta) based on precedence logic:
    
    1. If Alpha assignment is in expected taxa -> Keep Alpha
    2. Else if Beta assignment is in expected taxa -> Keep Beta
    3. Else if Alpha is Unclassified but Beta is Classified -> Keep Beta (Rescue)
    4. Else -> Keep Alpha (Default)
    
    Finally, marks for removal if the chosen taxonomy is not in the expected list.
    """
    
    # 1. Load Data
    alpha_data = load_taxonomy_dict(alpha_file)
    beta_data = load_taxonomy_dict(beta_file)
    
    target_names = set()
    if os.path.exists(lookup_file):
        with open(lookup_file, 'r') as f:
            for line in f:
                name = line.strip()
                if name: target_names.add(name)
    else:
        print(f"Warning: Lookup file not found: {lookup_file}")

    all_ids = set(alpha_data.keys()) | set(beta_data.keys())
    print(f"Merging {len(all_ids)} features...")
    
    counts = {
        'kept_alpha_expected': 0,
        'kept_beta_expected': 0,
        'rescued_beta': 0,
        'default_alpha': 0,
        'total_marked': 0
    }
    
    with open(output_file, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow(['Feature ID', 'Taxon', 'Confidence'])
        
        for fid in all_ids:
            # Get entries, defaulting to Unclassified if missing
            alpha_entry = alpha_data.get(fid, {'taxon': 'Unclassified', 'conf': '0.0'})
            beta_entry = beta_data.get(fid, {'taxon': 'Unclassified', 'conf': '0.0'})
            
            alpha_tax = alpha_entry['taxon']
            beta_tax = beta_entry['taxon']
            
            # Pre-calculate conditions
            alpha_is_expected = _lineage_contains_expected(alpha_tax, target_names)
            beta_is_expected = _lineage_contains_expected(beta_tax, target_names)
            alpha_is_unclassified = (alpha_tax == "Unclassified")
            beta_is_classified = (beta_tax != "Unclassified")

            # --- SELECTION LOGIC ---
            
            # Rule 1: If Alpha is in expected taxa -> Keep Alpha
            if alpha_is_expected:
                final_tax = alpha_tax
                final_conf = alpha_entry['conf']
                counts['kept_alpha_expected'] += 1
                
            # Rule 2: If not, look at Beta; if in expected taxa -> Keep Beta
            elif beta_is_expected:
                final_tax = beta_tax
                final_conf = beta_entry['conf']
                counts['kept_beta_expected'] += 1
                
            # Rule 3: If Alpha is unclassified but Beta has assignment -> Keep Beta
            # (Note: We only reach here if Beta was NOT expected, otherwise Rule 2 matches)
            elif alpha_is_unclassified and beta_is_classified:
                final_tax = beta_tax
                final_conf = beta_entry['conf']
                counts['rescued_beta'] += 1
                
            # Rule 4: Otherwise -> Keep Alpha
            else:
                final_tax = alpha_tax
                final_conf = alpha_entry['conf']
                counts['default_alpha'] += 1
            
            # --- FILTERING LOGIC ---
            
            should_remove = False
            
            if final_tax == "Unclassified":
                should_remove = True
            elif alpha_is_expected and final_tax == alpha_tax:
                should_remove = False # Rule 1
            elif beta_is_expected and final_tax == beta_tax:
                should_remove = False # Rule 2
            elif not _lineage_contains_expected(final_tax, target_names):
                should_remove = True # Rule 3 or 4 resulted in unexpected taxon
            
            if should_remove:
                final_tax = final_tax + "(remove)"
                counts['total_marked'] += 1
            
            writer.writerow([fid, final_tax, final_conf])
            
    print(f"Stats: Alpha Expected (Priority 1): {counts['kept_alpha_expected']}")
    print(f"Stats: Beta Expected (Priority 2): {counts['kept_beta_expected']}")
    print(f"Stats: Beta Rescued (Priority 3): {counts['rescued_beta']}")
    print(f"Stats: Alpha Default (Priority 4): {counts['default_alpha']}")
    print(f"Final Filter: {counts['total_marked']} features marked for removal.")
    
    return output_file


def generate_taxonomy_summary(alpha_file, beta_file, lookup_file, summary_output_file):
    """
    Compares Alpha and Beta taxonomies against an expected species list 
    and against each other, exporting ONLY the summary statistics CSV.
    Automatically creates the target directory chain if it is missing.
    """
    # 1. Load Data
    alpha_data = load_taxonomy_dict(alpha_file)
    beta_data = load_taxonomy_dict(beta_file)
    
    target_names = set()
    if os.path.exists(lookup_file):
        with open(lookup_file, 'r') as f:
            for line in f:
                name = line.strip()
                if name: target_names.add(name)
    else:
        print(f"Warning: Lookup file not found: {lookup_file}")

    all_ids = set(alpha_data.keys()) | set(beta_data.keys())
    
    summary_counts = {
        'alpha': {'unclassified': 0, 'correct_match': 0, 'mismatch_expected': 0, 'unexpected': 0},
        'beta': {'unclassified': 0, 'correct_match': 0, 'mismatch_expected': 0, 'unexpected': 0}
    }
    
    # 2. Process and Categorize
    for fid in all_ids:
        alpha_entry = alpha_data.get(fid, {'taxon': 'Unclassified', 'conf': '0.0'})
        beta_entry = beta_data.get(fid, {'taxon': 'Unclassified', 'conf': '0.0'})
        
        alpha_tax = alpha_entry['taxon']
        beta_tax = beta_entry['taxon']
        
        alpha_is_expected = _lineage_contains_expected(alpha_tax, target_names)
        beta_is_expected = _lineage_contains_expected(beta_tax, target_names)
        
        # Categorize Alpha
        if alpha_tax == "Unclassified":
            summary_counts['alpha']['unclassified'] += 1
        elif alpha_is_expected:
            if alpha_tax == beta_tax:
                summary_counts['alpha']['correct_match'] += 1
            else:
                summary_counts['alpha']['mismatch_expected'] += 1
        else:
            summary_counts['alpha']['unexpected'] += 1

        # Categorize Beta
        if beta_tax == "Unclassified":
            summary_counts['beta']['unclassified'] += 1
        elif beta_is_expected:
            if alpha_tax == beta_tax:
                summary_counts['beta']['correct_match'] += 1
            else:
                summary_counts['beta']['mismatch_expected'] += 1
        else:
            summary_counts['beta']['unexpected'] += 1

    # 3. Secure Directory & Write Summary CSV File
    parent_dir = os.path.dirname(summary_output_file)
    if parent_dir: 
        os.makedirs(parent_dir, exist_ok=True)

    dataset_name = "_".join(os.path.basename(summary_output_file).split("-")[1:3])
    
    
    with open(summary_output_file, 'w', newline='') as f_sum:
        sum_writer = csv.writer(f_sum, delimiter=',')
        sum_writer.writerow([dataset_name, 'Alpha Taxonomy', 'Beta Taxonomy'])
        
        sum_writer.writerow([
            'Correct', #'correctly classified and matching the other alpha/beta taxonomy', 
            summary_counts['alpha']['correct_match'] + summary_counts['alpha']['mismatch_expected'], 
            summary_counts['beta']['correct_match']+ summary_counts['beta']['mismatch_expected']
        ])

        sum_writer.writerow([
            'Offtarget/False positive', #'classified but not present in the expected species list', 
            summary_counts['alpha']['unexpected'], 
            summary_counts['beta']['unexpected']
        ])
        sum_writer.writerow([
            'Unclassified', 
            summary_counts['alpha']['unclassified'], 
            summary_counts['beta']['unclassified']
        ])

    print(f"Summary metrics successfully compiled and exported to: {summary_output_file}")
    return summary_output_file