#!/usr/bin/env python3
"""
Convert RDP results to QIIME-compatible taxonomy file.
Uses taxonomy_utils for all taxonomy interactions.
"""

import os
# Import logic from shared module
from taxonomy_utils import (
    NCBITaxonomyLookup,
    get_hierarchical_common_names,
    get_common_names_from_ncbi,
    format_qiime_taxonomy
)

def extract_species_from_rdp(line, threshold=0.8):
    """Parse RDP line to find best assignment."""
    parts = line.strip().split('\t')
    best_tax, best_conf = "unassigned", 0.0
    full_path = []
    
    # RDP columns start at index 8 and come in triplets (Taxon, Rank, Conf)
    i = 8
    while i < len(parts) - 2:
        name, _, conf_str = parts[i], parts[i+1], parts[i+2]
        conf = float(conf_str)
        if name and not name.startswith("undef_"):
            if conf >= threshold:
                full_path.append(name)
                best_tax, best_conf = name, conf
        i += 3
    return best_tax, best_conf, "; ".join(full_path)

def rdp_to_qiime_standardized(input_file, output_file, names_dmp_path=None, 
                            confidence_threshold=0.9, use_local_db=True,
                            include_common_names=True, hierarchical_names=True):
    
    # 1. Initialize shared Lookup Class
    lookup = NCBITaxonomyLookup(names_dmp_path, use_local_db)
    
    # Data holders
    rows = [] # (fid, taxid, rdp_tax, conf)
    unique_taxids = set()
    
    # 2. First Pass: Parse RDP & Find TaxIDs
    print("Parsing RDP file...")
    with open(input_file, 'r') as f:
        for line in f:
            fid = line.split('\t')[0]
            species, conf, rdp_tax = extract_species_from_rdp(line, confidence_threshold)
            
            taxid = None
            if species != "unassigned":
                # Use shared lookup
                taxid = lookup.search_taxid(species)
                if taxid: unique_taxids.add(taxid)
            
            rows.append((fid, taxid, rdp_tax, conf))
            
    # 3. Batch Lookup Common Names (Shared Logic)
    common_map = {}
    if include_common_names and unique_taxids:
        print(f"Fetching common names for {len(unique_taxids)} taxa...")
        if hierarchical_names:
            common_map = get_hierarchical_common_names(unique_taxids, names_dmp_path)
        else:
            common_map = get_common_names_from_ncbi(unique_taxids, names_dmp_path)

    # 4. Second Pass: Write Output
    print("Writing output...")
    with open(output_file, 'w') as out:
        out.write("Feature ID\tTaxon\tConfidence\n")
        
        for fid, taxid, rdp_tax, conf in rows:
            final_str = "Unclassified"
            
            # If we found a valid NCBI TaxID, standardize the lineage
            if taxid:
                # Use shared lookup to get [Kingdom, Phylum...] list
                names_list = lookup.get_lineage_list(taxid)
                # Use shared formatter
                final_str = format_qiime_taxonomy(names_list, common_map.get(taxid))
            elif conf >= confidence_threshold:
                # Fallback to RDP string if confident but no NCBI match
                final_str = rdp_tax
                
            out.write(f"{fid}\t{final_str}\t{conf:.6f}\n")

def batch_process_rdp_files(input_dir, output_dir, years, **kwargs):
    os.makedirs(output_dir, exist_ok=True)
    for year in years:
        infile = f"{input_dir}/Musq-{year}-COI-RDP-results"
        outfile = f"{output_dir}/Musq-{year}-COI-RDP-standardized.qiime"
        if os.path.exists(infile):
            rdp_to_qiime_standardized(infile, outfile, **kwargs)
            print(f"Finished {year}")

if __name__ == "__main__":
    # Example trigger
    pass