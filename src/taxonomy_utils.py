#!/usr/bin/env python3
"""
Shared utilities for taxonomy processing.
Contains common logic for:
1. NCBI Name Lookups (names.dmp)
2. QIIME String Formatting
3. TaxonKit CLI wrapping (Batch & Single)
4. Unified Taxonomy Lookup (Local DB -> API -> CLI)
"""

import os
import subprocess
import tempfile
import requests
import time
from collections import defaultdict

# ==============================================================================
# 1. STRING FORMATTING
# ==============================================================================

def format_qiime_taxonomy(names_list, common_name=None):
    """
    Construct a QIIME-compatible lineage string from a list of names.
    """
    prefixes = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    qiime_parts = []
    
    length = min(len(names_list), len(prefixes))
    
    for i in range(len(prefixes)):
        if i < length and names_list[i] and names_list[i].strip():
            name = names_list[i].strip()
            if i == len(prefixes) - 1 and common_name:
                qiime_parts.append(f"{prefixes[i]}{name} ({common_name})")
            else:
                qiime_parts.append(f"{prefixes[i]}{name}")
        else:
            qiime_parts.append(prefixes[i])
            
    return ';'.join(qiime_parts)

# ==============================================================================
# 2. COMMON NAME LOOKUPS
# ==============================================================================

def get_common_names_from_ncbi(taxids, names_dmp_path=None):
    """Get common names for taxids directly from NCBI names.dmp file."""
    if not names_dmp_path:
        names_dmp_path = _find_file('names.dmp')
        
    if not names_dmp_path or not os.path.exists(names_dmp_path):
        return {}
    
    taxid_to_common = {}
    try:
        with open(names_dmp_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                fields = line.strip().split('\t|\t')
                if len(fields) >= 4:
                    try:
                        taxid = int(fields[0])
                        if taxid in taxids:
                            name_class = fields[3].rstrip('\t|').strip()
                            if name_class in ['genbank common name', 'common name', 'blast name']:
                                # Prioritize genbank common name
                                if name_class == 'genbank common name':
                                    taxid_to_common[taxid] = fields[1].strip()
                                elif taxid not in taxid_to_common:
                                    taxid_to_common[taxid] = fields[1].strip()
                    except (ValueError, IndexError):
                        continue
    except Exception as e:
        print(f"Warning reading names.dmp: {e}")
    return taxid_to_common


def get_hierarchical_common_names(taxids, names_dmp_path=None, nodes_dmp_path=None):
    """Get common names by looking up the taxonomic hierarchy."""
    if not names_dmp_path:
        names_dmp_path = _find_file('names.dmp')
    if not nodes_dmp_path and names_dmp_path:
        nodes_dmp_path = names_dmp_path.replace('names.dmp', 'nodes.dmp')
        
    if not names_dmp_path or not os.path.exists(names_dmp_path): return {}
    
    # 1. Load Hierarchy
    taxid_to_parent = {}
    if nodes_dmp_path and os.path.exists(nodes_dmp_path):
        try:
            with open(nodes_dmp_path, 'r') as f:
                for line in f:
                    fields = line.strip().split('\t|\t')
                    if len(fields) >= 2:
                        taxid_to_parent[int(fields[0])] = int(fields[1])
        except Exception: pass
        
    # 2. Load All Common Names
    all_common = {}
    try:
        with open(names_dmp_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                fields = line.strip().split('\t|\t')
                if len(fields) >= 4:
                    name_class = fields[3].rstrip('\t|')
                    if name_class in ['genbank common name', 'common name', 'blast name']:
                        tid = int(fields[0])
                        if tid not in all_common or name_class == 'genbank common name':
                            all_common[tid] = fields[1]
    except Exception: pass
    
    # 3. Traverse
    results = {}
    for tid in taxids:
        curr = tid
        visited = set()
        while curr and curr not in visited:
            visited.add(curr)
            if curr in all_common:
                results[tid] = all_common[curr]
                break
            curr = taxid_to_parent.get(curr)
            if curr == taxid_to_parent.get(curr): break 
    return results

def _find_file(filename):
    opts = [filename, f'~/.taxonkit/{filename}', f'/usr/local/share/taxonkit/{filename}']
    for p in opts:
        ep = os.path.expanduser(p)
        if os.path.exists(ep): return ep
    return None

# ==============================================================================
# 3. TAXONKIT WRAPPER
# ==============================================================================

def run_taxonkit_batch(taxids, data_dir=None):
    """Run taxonkit for a list of taxids."""
    if not taxids: return {}
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp_in:
        for t in taxids: tmp_in.write(f"{t}\n")
        tmp_in_path = tmp_in.name
        
    try:
        cmd = ['taxonkit', 'lineage']
        if data_dir: cmd.extend(['--data-dir', os.path.dirname(data_dir)])
        cmd.append(tmp_in_path)
        
        res1 = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp_mid:
            tmp_mid.write(res1.stdout)
            tmp_mid_path = tmp_mid.name
            
        fmt = "{K};{p};{c};{o};{f};{g};{s}"
        cmd2 = ['taxonkit', 'reformat', '-F', '-f', fmt, tmp_mid_path]
        res2 = subprocess.run(cmd2, capture_output=True, text=True, check=True)
        
        results = {}
        for line in res2.stdout.strip().split('\n'):
            if not line.strip(): continue
            parts = line.split('\t')
            if len(parts) >= 3:
                try:
                    tid = int(parts[0])
                    lineage_str = parts[2] if len(parts) > 2 else parts[1]
                    results[tid] = lineage_str.split(';')
                except ValueError: pass
        
        if os.path.exists(tmp_mid_path): os.unlink(tmp_mid_path)
        return results
        
    except subprocess.CalledProcessError as e:
        print(f"TaxonKit Error: {e.stderr}")
        return {}
    finally:
        if os.path.exists(tmp_in_path): os.unlink(tmp_in_path)

# ==============================================================================
# 4. UNIFIED TAXONOMY LOOKUP CLASS (UPDATED)
# ==============================================================================

class NCBITaxonomyLookup:
    """
    Unified interface for taxonomy lookups.
    Handles 'Genus_species' -> 'Genus species' conversion automatically.
    """
    def __init__(self, names_dmp_path=None, use_local_db=True):
        self.names_dmp_path = names_dmp_path
        self.nodes_dmp_path = names_dmp_path.replace('names.dmp', 'nodes.dmp') if names_dmp_path else None
        
        self.ncbi = None
        if use_local_db:
            try:
                from ete3 import NCBITaxa
                self.ncbi = NCBITaxa()
            except Exception:
                print("Warning: ETE3 not available or DB init failed. Using API/CLI.")
    
    def search_taxid(self, name):
        """Find taxid by name (exact or fuzzy)."""
        clean_name = name.strip()
        
        # Define variations to try (Original, Space replaced, Underscore replaced)
        variations = [clean_name]
        if '_' in clean_name:
            variations.append(clean_name.replace('_', ' ')) # Homo_sapiens -> Homo sapiens
        # if ' ' in clean_name:
        #     variations.append(clean_name.replace(' ', '_')) # Homo sapiens -> Homo_sapiens
            
        # 1. Try ETE3 (Local)
        if self.ncbi:
            try:
                for var in variations:
                    n2t = self.ncbi.get_name_translator([var])
                    if var in n2t: return n2t[var][0]
            except Exception: pass
            
        # 2. Try API (Online)
        try:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            for var in variations:
                params = {'db': 'taxonomy', 'term': var, 'retmode': 'xml', 'retmax': 1}
                res = requests.get(url, params=params)
                if '<Id>' in res.text:
                    return int(res.text.split('<Id>')[1].split('</Id>')[0])
                time.sleep(0.1) # Be nice to NCBI
        except Exception: pass
        
        return None

    def get_lineage_list(self, taxid):
        """Returns list [Kingdom, Phylum, ..., Species] for a taxid."""
        if not taxid: return []
        
        # 1. Try ETE3
        if self.ncbi:
            try:
                lineage = self.ncbi.get_lineage(taxid)
                names = self.ncbi.get_taxid_translator(lineage)
                ranks = self.ncbi.get_rank(lineage)
                std_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
                r_to_n = {ranks[t]: names[t] for t in lineage}
                if 'kingdom' in r_to_n and 'superkingdom' not in r_to_n: 
                    r_to_n['superkingdom'] = r_to_n['kingdom']
                return [r_to_n.get(r, '') for r in std_ranks]
            except: pass
            
        # 2. Fallback to TaxonKit
        batch_res = run_taxonkit_batch([taxid], self.names_dmp_path)
        return batch_res.get(taxid, [])
    
    def get_common_name(self, species_name, hierarchical=False):
        """
        Get common name for a species by name string.
        
        Args:
            species_name: Species name (e.g., 'Homo sapiens' or 'Homo_sapiens')
            hierarchical: If True, walks up the taxonomy tree to find a common name
                        if the species itself has none.
        
        Returns:
            Common name string, or None if not found.
        """
        # Step 1: Resolve name -> taxid
        taxid = self.search_taxid(species_name)
        if not taxid:
            print(f"Warning: Could not find taxid for '{species_name}'")
            return None

        # Step 2: Try ETE3 first (fastest, local)
        if self.ncbi:
            try:
                lineage = self.ncbi.get_lineage(taxid) if hierarchical else [taxid]
                for tid in lineage[::-1]:  # Walk from species -> root
                    names = self.ncbi.get_common_names([tid]) if hasattr(self.ncbi, 'get_common_names') else {}
                    # ETE3 stores common names in its SQLite DB
                    result = self.ncbi.get_taxid_translator([tid])
                    # Fall through to names.dmp lookup for common names
                # ETE3 doesn't expose common names directly, so fall through
            except Exception:
                pass

        # Step 3: Use names.dmp lookups
        taxids_set = {taxid}
        if hierarchical:
            lookup_fn = get_hierarchical_common_names
            result = lookup_fn(
                taxids_set,
                names_dmp_path=self.names_dmp_path,
                nodes_dmp_path=self.nodes_dmp_path
            )
        else:
            lookup_fn = get_common_names_from_ncbi
            result = lookup_fn(
                taxids_set,
                names_dmp_path=self.names_dmp_path
            )

        if taxid in result:
            return result[taxid]

        # Step 4: Fallback to NCBI API
        try:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            params = {'db': 'taxonomy', 'id': taxid, 'retmode': 'json'}
            res = requests.get(url, params=params).json()
            summary = res.get('result', {}).get(str(taxid), {})
            common = summary.get('commonname', '').strip()
            if common:
                return common
        except Exception:
            pass

        return None
    
    def get_common_names_batch(self, species_list, hierarchical=False):
        """
        Look up common names for a list of species in one pass.
        Much faster than calling get_common_name() 70 times.
        
        Returns: dict {species_name: common_name_or_None}
        """
        results = {sp: None for sp in species_list}

        # Step 1: Resolve all names -> taxids in bulk via ETE3
        name_to_taxid = {}
        if self.ncbi:
            try:
                clean = {sp: sp.replace('_', ' ') for sp in species_list}
                translated = self.ncbi.get_name_translator(list(clean.values()))
                for sp, clean_name in clean.items():
                    if clean_name in translated:
                        name_to_taxid[sp] = translated[clean_name][0]
            except Exception:
                pass

        # API fallback only for species ETE3 missed
        missing = [sp for sp in species_list if sp not in name_to_taxid]
        for sp in missing:
            print(f"API fallback for sp taxid: {sp}")
            taxid = self.search_taxid(sp)  # still per-species but only for misses
            if taxid:
                name_to_taxid[sp] = taxid

        if not name_to_taxid:
            return results

        # Step 2: Single pass through names.dmp for ALL taxids at once
        taxid_set = set(name_to_taxid.values())
        if hierarchical:
            print(f"getting hierarchial commonname: {taxid_set}")
            common_by_taxid = get_hierarchical_common_names(
                taxid_set,
                names_dmp_path=self.names_dmp_path,
                nodes_dmp_path=self.nodes_dmp_path
            )
        else:
            print(f"getting  name from ncbi: {taxid_set}")
            common_by_taxid = get_common_names_from_ncbi(
                taxid_set,
                names_dmp_path=self.names_dmp_path
            )

        # Step 3: API fallback only for taxids still missing a common name
        still_missing_taxids = taxid_set - set(common_by_taxid.keys())
        if still_missing_taxids:
            try:
                ids_str = ','.join(str(t) for t in still_missing_taxids)
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                res = requests.get(url, params={'db': 'taxonomy', 'id': ids_str, 'retmode': 'json'}).json()
                for taxid in still_missing_taxids:
                    print(f"still missing id for: {taxid}")
                    common = res.get('result', {}).get(str(taxid), {}).get('commonname', '').strip()
                    if common:
                        common_by_taxid[taxid] = common
            except Exception:
                pass

        # Step 4: Map back to species names
        for sp, taxid in name_to_taxid.items():
            results[sp] = common_by_taxid.get(taxid)

        return results
    
    def get_ranks_batch(self, species_list, rank='class'):
        """
        Get a specific taxonomic rank for a list of species names.
        
        Args:
            species_list: List of species names (e.g., ['Homo sapiens', 'Panthera_leo'])
            rank: Taxonomic rank to retrieve. One of:
                'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'
        
        Returns:
            dict: {species_name: rank_value_or_None}
        """
        STANDARD_RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        if rank not in STANDARD_RANKS:
            raise ValueError(f"Invalid rank '{rank}'. Must be one of: {STANDARD_RANKS}")

        results = {sp: None for sp in species_list}

        # Step 1: Resolve all species names -> taxids in bulk via ETE3
        name_to_taxid = {}
        if self.ncbi:
            try:
                clean = {sp: sp.replace('_', ' ') for sp in species_list}
                translated = self.ncbi.get_name_translator(list(clean.values()))
                for sp, clean_name in clean.items():
                    if clean_name in translated:
                        name_to_taxid[sp] = translated[clean_name][0]
            except Exception:
                pass

        # API fallback for any species ETE3 missed
        missing = [sp for sp in species_list if sp not in name_to_taxid]
        for sp in missing:
            taxid = self.search_taxid(sp)
            if taxid:
                name_to_taxid[sp] = taxid

        if not name_to_taxid:
            return results

        # Step 2: Resolve lineages in bulk
        rank_idx = STANDARD_RANKS.index(rank)
        taxid_to_rank_value = {}

        if self.ncbi:
            # ETE3 path: resolve all lineages at once
            all_taxids = list(set(name_to_taxid.values()))
            try:
                for taxid in all_taxids:
                    lineage = self.ncbi.get_lineage(taxid)
                    ranks = self.ncbi.get_rank(lineage)
                    names = self.ncbi.get_taxid_translator(lineage)

                    # Build rank -> name map
                    r_to_n = {ranks[t]: names[t] for t in lineage}
                    if 'kingdom' in r_to_n and 'superkingdom' not in r_to_n:
                        r_to_n['superkingdom'] = r_to_n['kingdom']

                    taxid_to_rank_value[taxid] = r_to_n.get(rank)
            except Exception:
                pass

        # TaxonKit fallback for taxids ETE3 couldn't resolve
        missing_taxids = [t for t in set(name_to_taxid.values()) if t not in taxid_to_rank_value]
        if missing_taxids:
            batch_res = run_taxonkit_batch(missing_taxids, self.names_dmp_path)
            for taxid, lineage_list in batch_res.items():
                value = lineage_list[rank_idx] if rank_idx < len(lineage_list) else None
                taxid_to_rank_value[taxid] = value or None  # Treat empty string as None

        # Step 3: Map taxids back to original species names
        for sp, taxid in name_to_taxid.items():
            results[sp] = taxid_to_rank_value.get(taxid)

        return results
