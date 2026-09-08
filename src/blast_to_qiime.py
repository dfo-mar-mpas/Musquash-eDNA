#!/usr/bin/env python3
"""
Convert BLAST results to QIIME-compatible taxonomy file.
Uses taxonomy_utils for core logic.
"""

import os
from collections import defaultdict
# Import everything from our new robust module
from taxonomy_utils import (
    get_hierarchical_common_names, 
    get_common_names_from_ncbi,
    run_taxonkit_batch,
    format_qiime_taxonomy
)

def parse_fasta_ids(fasta_file):
    """Extract IDs from FASTA."""
    ids = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                ids.add(line[1:].strip().split()[0])
    return ids

def extract_taxids_from_blast(blast_file, use_top_scoring=False):
    """Extract taxids from BLAST results."""
    taxids = set()
    feature_data = defaultdict(list)
    
    with open(blast_file, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 6: continue
            
            feat_id = fields[0]
            # Standard metrics extraction...
            metrics = {
                'pident': float(fields[2]) if len(fields)>2 else 0,
                'evalue': float(fields[4]) if len(fields)>4 else 1,
                'bitscore': float(fields[11]) if len(fields)>11 else 0
            }
            
            # 1. Try TaxID column (MIDORI2)
            try:
                taxid = int(fields[6])
                if taxid == 0: raise ValueError
                taxids.add(taxid)
                metrics['taxid'] = taxid
                feature_data[feat_id].append(metrics if use_top_scoring else taxid)
                continue
            except ValueError: pass
            
            # 2. Try Lineage String (nt_euk)
            try:
                lineage_part = fields[1].split('###')[1].split(';')[-1] # "...###...;s_Species_12345"
                taxid = int(lineage_part.split('_')[-1])
                taxids.add(taxid)
                metrics['taxid'] = taxid
                feature_data[feat_id].append(metrics if use_top_scoring else taxid)
            except (IndexError, ValueError): pass
            
    return taxids, feature_data

def generate_qiime_lines(feature_data, taxid_to_lineage, query_ids, use_top_scoring, scoring_method):
    lines = ["Feature ID\tTaxon\tConfidence"]
    processed = set()
    
    for fid, data in feature_data.items():
        # --- SCORING LOGIC ---
        if use_top_scoring:
            # Sort by chosen metric
            key = lambda x: x['bitscore'] # Default
            if scoring_method == 'pident': key = lambda x: x['pident']
            elif scoring_method == 'evalue_inv': key = lambda x: 1/(x['evalue']+1e-300)
            
            valid_hits = [x for x in data if x['taxid'] in taxid_to_lineage]
            if not valid_hits:
                taxon, conf = "Unclassified", 0.0
            else:
                top = sorted(valid_hits, key=key, reverse=True)[0]
                top_tax = taxid_to_lineage[top['taxid']]
                count = sum(1 for x in valid_hits if taxid_to_lineage.get(x['taxid']) == top_tax)
                taxon, conf = top_tax, count/len(valid_hits)
        else:
            # Consensus
            valid_taxes = [taxid_to_lineage[t] for t in data if t in taxid_to_lineage]
            if not valid_taxes:
                taxon, conf = "Unclassified", 0.0
            else:
                counts = defaultdict(int)
                for t in valid_taxes: counts[t] += 1
                taxon = max(counts, key=counts.get)
                conf = counts[taxon] / len(valid_taxes)
        
        if taxon != "Unclassified":
            lines.append(f"{fid}\t{taxon}\t{conf:.6f}")
            processed.add(fid)
            
    for qid in query_ids:
        if qid not in processed: lines.append(f"{qid}\tUnclassified\t0.000000")
        
    return lines

def convert_blast_to_qiime(blast_path, fasta_path, output_path, **kwargs):
    # Unpack kwargs
    names_dmp = kwargs.get('names_dmp_path')
    include_common = kwargs.get('include_common_names', False)
    hierarchical = kwargs.get('use_hierarchical_names', True)
    
    # 1. Parse
    q_ids = parse_fasta_ids(fasta_path)
    taxids, feat_data = extract_taxids_from_blast(blast_path, kwargs.get('use_top_scoring', False))
    
    # 2. Get Raw Lineages (List of names)
    # Using the shared Batch function!
    print("Getting lineages via TaxonKit...")
    raw_lineages = run_taxonkit_batch(taxids, names_dmp)
    
    # 3. Get Common Names
    common_names = {}
    if include_common:
        if hierarchical:
            common_names = get_hierarchical_common_names(taxids, names_dmp)
        else:
            common_names = get_common_names_from_ncbi(taxids, names_dmp)
            
    # 4. Format Lineages
    # Using the shared Format function!
    final_lineages = {}
    for tid, names in raw_lineages.items():
        cname = common_names.get(tid)
        final_lineages[tid] = format_qiime_taxonomy(names, cname)
        
    # 5. Generate Output
    lines = generate_qiime_lines(feat_data, final_lineages, q_ids, 
                               kwargs.get('use_top_scoring'), kwargs.get('scoring_method'))
    
    with open(output_path, 'w') as f: f.write('\n'.join(lines))
    print(f"Done: {output_path}")

def convert_blast_to_qiime_dir(input_dir, blast, fasta, output, **kwargs):
    convert_blast_to_qiime(os.path.join(input_dir, blast), 
                          os.path.join(input_dir, fasta), 
                          os.path.join(input_dir, output), **kwargs)