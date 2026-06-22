#!/usr/bin/env python3
from __future__ import annotations

import multiprocessing as mp
import os
import sys
import logging
import pandas as pd
from typing import Iterable, List, Optional, Tuple
from collections import defaultdict
import pysam

# Global BAM handle and config for worker processes
taxa_dict = {}
lineage_cache = {}  # Cache for reference taxid to qualified taxids mapping
ref_to_extract_taxid = defaultdict(list)  # Mapping from reference taxid to set of qualified taxids to extract

def load_criteria_from_log(gottcha_log: str) -> Tuple[float, float, int, str]:
    """Parse match parameters from the gottcha log file."""
    try:
        with open(gottcha_log) as f:
            for line in f:
                if "Min Match Identity" in line:
                    match_identity = float(line.split(" ")[-1].strip())
                elif "Min Match Fraction" in line:
                    match_fraction = float(line.split(" ")[-1].strip())
                elif "Min Match Length" in line:
                    match_length = int(line.split(" ")[-1].strip())
                elif "SNI-score" in line:
                    sniScore = line.split(" ")[-1].strip()
            return (match_identity, match_fraction, match_length, sniScore)
    except Exception as e:
        logging.error(f"Error parsing gottcha log file {gottcha_log}: {e}")
        sys.exit(1)


def parse_taxids(taxid_arg: str, 
                 res_df: pd.DataFrame, 
                 full_tsv_fn: str,
                 sni_score_cutoff: float,
                 sni_score_species: float,
                 sni_score_strain: float) -> Tuple[dict, list]:
    """Parse taxids from command line arg or file"""

    global ref_to_extract_taxid
    taxa_df = pd.DataFrame()

    # Load the full report TSV if available, otherwise use the res_df from memory
    if res_df.shape[1] > 0:
        taxa_df = res_df
        logging.info(f"Successfully loaded result summary with {len(taxa_df)} taxonomic entries")
    else:
        try:
            logging.info(f"Reading result summary file {full_tsv_fn}...")
            taxa_df = pd.read_csv(full_tsv_fn,
                                  sep='\t',
                                  quoting=3,
                                  on_bad_lines='skip',
                                  dtype={'NOTE': str, 'TAXID': str, 'PARENT_TAXID': str})

            logging.info(f"Successfully loaded result summary with {len(taxa_df)} entries")
        except Exception as e:
            logging.error(f"Error reading result summary file: {e}")
            sys.exit(1)

    # get full lineages of each strain-level taxid in the report
    df = taxa_df[['LEVEL','NAME', 'SNI_SCORE', 'TAXID', 'PARENT_TAXID']]
    df_lineages = df[df['LEVEL']=='strain'].copy().reset_index(drop=True)
    df_lineages = df_lineages.rename(columns={'TAXID':'strain', 'PARENT_TAXID':'species'})

    RANKS = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]

    for idx, RANK in enumerate(RANKS):
        rank_df = df[df['LEVEL']==RANK].reset_index(drop=True).copy()
        rank_df[RANK] = rank_df['TAXID']
        df_lineages = df_lineages.merge(rank_df[[RANK,'PARENT_TAXID']], on=RANK, how='left')

        if idx == len(RANKS)-1:
            next_rank = 'root'
        else:
            next_rank = RANKS[idx+1]
        df_lineages = df_lineages.rename(columns={'PARENT_TAXID':next_rank})

    logging.debug(f"Initial lineages dataframe:\n{df_lineages.head()}")

    # expand the input taxa for extracction to the taxid in the restuls, and store to qualified_taxa
    qualified_idx = pd.Series([True]*len(taxa_df))

    if taxid_arg and taxid_arg != 'all':
        if taxid_arg.startswith('@'):
            # Read taxids from file
            filename = taxid_arg[1:]  # Remove @ prefix
            try:
                with open(filename) as f:
                    taxa_list = [x.strip() for x in f.readlines() if x.strip() and not x.startswith('#')]
            except IOError as e:
                logging.error(f"Error reading taxid file {filename}: {e}")
                sys.exit(1)
        else:
            # Parse comma-separated list
            taxa_list = [x.strip() for x in taxid_arg.split(',')]

        if taxa_list:
            qualified_idx &= (taxa_df['TAXID'].isin(taxa_list) | taxa_df['NAME'].isin(taxa_list))

    logging.info(f"Found {qualified_idx.sum()} qualified taxa after filtering by taxid/name")

    # Filter out entries with "Filtered out" notes
    if qualified_idx is not None:
        qualified_taxa = taxa_df.loc[qualified_idx, ['LEVEL', 'TAXID']].copy()
    else:
        qualified_taxa = taxa_df[['LEVEL', 'TAXID']].copy()

    # Further filter by SNI score thresholds
    for level, extract_taxid in qualified_taxa.itertuples(index=False):
        # set sni cutoff
        if level == 'strain':
            sni_score_cutoff = sni_score_strain
        elif level == 'species':
            sni_score_cutoff = sni_score_species
        
        logging.debug(f"Processing - level: {level}; taxid: {extract_taxid}; sni_cutoff: {sni_score_cutoff}")

        idx = (df_lineages[level]==extract_taxid) & (df_lineages['SNI_SCORE']>=sni_score_cutoff)
        for ref_taxid in df_lineages[idx]['strain']:
            ref_to_extract_taxid[ref_taxid].append(extract_taxid)

    if len(ref_to_extract_taxid) == 0:
        logging.warning("No qualified taxa for extraction. No reads can be extracted.")
        sys.exit(0)

    return taxa_df.set_index('TAXID')[['LEVEL','NAME']].to_dict(orient='index'), ref_to_extract_taxid


def _init_worker(bam_path: str,
                 taxa_dict: dict,
                 format: str,
                 min_frac: float,
                 min_idt: float,
                 min_alen: int,
                 min_mapq: int = 0,
                 htslib_threads: int = 1,
                 include_secondary: bool = False,
                 include_supplementary: bool = False,
                 include_duplicates: bool = False,
                 include_qcfail: bool = False) -> None:
    """Initializer for each worker process: open BAM once and stash filters."""
    global _BAM, _CFG
    _BAM = pysam.AlignmentFile(bam_path, "rb", threads=htslib_threads)
    _CFG = {
        "min_mapq": min_mapq,
        "taxa_dict": taxa_dict,
        "min_frac": min_frac,
        "min_idt": min_idt,
        "min_alen": min_alen,
        "include_secondary": include_secondary,
        "include_supplementary": include_supplementary,
        "include_duplicates": include_duplicates,
        "include_qcfail": include_qcfail,
        "format": format
    }


def _iter_tasks(references: List[str],
                max_per_taxon: int,
                acc_list: list, 
                acc_list_action: str
                ) -> Iterable[Tuple[str, int, int]]:
    
    global ref_to_extract_taxid

    for ref in references:
        # Extract taxid from reference
        try:
            acc, rstart, rend, ref_taxid, _ = ref.split('|')
        except ValueError:
            logging.debug(f"Malformed reference: {ref}")
            continue  # Skip malformed references

        if ref_taxid not in ref_to_extract_taxid:
            continue  # Skip references that don't belong to any qualified taxon

        # Skip if accession is in the exclusion list (if applicable)
        aoi_flag = False
        if acc_list:
            if (acc in acc_list) or (ref in acc_list):
                aoi_flag = True
                if acc_list_action == 'filter_out':
                    continue
            else:
                if acc_list_action == 'filter_in':
                    continue

        if len(ref_to_extract_taxid[ref_taxid]) > 0:
            yield (ref, max_per_taxon, aoi_flag)


def _extract_worker(task: Tuple[str, int, bool]) -> dict:
    """
    Process one (rname, start0, end0) chunk.

    Returns:
      (rname, start0, end0, numreads, covbases, mismatches_total,
       consensus_diff, mean_depth)
    """
    global _BAM, _CFG, ref_to_extract_taxid
    assert _BAM is not None, "Worker BAM handle not initialized"
    
    taxon_seqs = defaultdict(list)  # Dictionary to hold sequences for each taxid
    ref, max_per_taxon, aoi_flag = task
    
    try:
        ref_taxid = ref.split('|')[3]
    except ValueError:
        logging.debug(f"Malformed reference: {ref}")

    bam = _BAM

    min_mapq = _CFG["min_mapq"]
    taxa_dict = _CFG["taxa_dict"]
    format = _CFG["format"]
    min_frac = _CFG["min_frac"]
    min_idt = _CFG["min_idt"]
    min_alen = _CFG["min_alen"]
    inc_sec = _CFG["include_secondary"]
    inc_sup = _CFG["include_supplementary"]
    inc_dup = _CFG["include_duplicates"]
    inc_qcf = _CFG["include_qcfail"]

    # Iterate reads overlapping this region.
    for aln in bam.fetch(ref):
        # Basic filters
        if aln.is_unmapped:
            continue
        if (not inc_sec) and aln.is_secondary:
            continue
        if (not inc_sup) and aln.is_supplementary:
            continue
        if (not inc_dup) and aln.is_duplicate:
            continue
        if (not inc_qcf) and aln.is_qcfail:
            continue
        if aln.mapping_quality < min_mapq:
            continue

        if min_idt > 0.0 and aln.has_tag('NM'):
            mm_idt = aln.get_tag('NM') / aln.alen
            if min_idt > (1-mm_idt):
                continue

        if min_frac > 0.0:
            if (aln.alen / aln.query_length) < min_frac and (aln.alen / bam.get_reference_length(ref)) < min_frac:
                continue

        if min_alen > 0 and aln.alen < min_alen:
            continue

        # Process the matching taxa
        for extract_taxid in ref_to_extract_taxid[ref_taxid]:

            # Only collect up to max_per_taxon sequences per taxon
            if (max_per_taxon==0) or (len(taxon_seqs[extract_taxid]) < max_per_taxon):
                # Create FASTA entry with taxonomy information
                level = taxa_dict[extract_taxid]['LEVEL']
                name = taxa_dict[extract_taxid]['NAME'].replace(' ', '_')  # Replace spaces with underscores for FASTA headers
                rname = aln.query_name
                region = (aln.reference_start+1, aln.reference_end)
                mapping_idt = 1 - (aln.get_tag('NM') / aln.alen) if aln.has_tag('NM') else 0
                mapping_frac = max((aln.alen / aln.query_length), (aln.alen / bam.get_reference_length(ref))) if bam.get_reference_length(ref) > 0 else 0

                # determine if the read is the first or second mate
                mate = ''
                if aln.is_paired:
                    mate = '.1' if aln.is_read1 else '.2'

                if format == 'fasta':
                    fasta_entry = f">{rname}{mate}|{ref}:{region[0]}..{region[1]} LEVEL={level} NAME={name} TAXID={extract_taxid} AOI={aoi_flag} MG={aln.alen} MI={mapping_idt:.2f} MF={mapping_frac:.2f}\n{aln.query_sequence}\n"
                else:
                    fasta_entry = f"@{rname}{mate}|{ref}:{region[0]}..{region[1]} LEVEL={level} NAME={name} TAXID={extract_taxid} AOI={aoi_flag} MG={aln.alen} MI={mapping_idt:.2f} MF={mapping_frac:.2f}\n{aln.query_sequence}\n+\n{aln.query_qualities_str}\n"
                
                taxon_seqs[extract_taxid].append(fasta_entry)

    return taxon_seqs


def extract_sequences_by_taxonomy(bam_path: str,
                                  taxa_dict: dict,
                                  ref_to_extract_taxid: dict,
                                  o,
                                  numthreads: int,
                                  matchFraction: float,
                                  matchIdentity: float,
                                  matchLength: int,
                                  max_per_taxon: int,
                                  acc_list: list,
                                  acc_list_action: str,
                                  format: str = 'fasta'):
    """
    Extract sequences mapping to taxa from the full taxonomy report.

    For each taxon in the full report, extract up to max_per_taxon sequences.

    Parameters:
        bam_path (str): Path to the BAM file
        taxa_dict (dict): Dictionary containing taxonomy information
        qualified_taxids (list): List of qualified taxonomic IDs to extract
        o (file): Output file handle for the extracted sequences
        numthreads (int): Number of threads to use for processing
        matchFraction (float): Minimum fraction required for a valid match
        matchIdentity (float): Minimum identity required for a valid match
        matchLength (int): Minimum length required for a valid match
        max_per_taxon (int): Maximum number of sequences to extract per taxon; 0 is unlimited.
        format (str): Output format ('fasta' or 'fastq')

    Returns:
        tuple: (taxon_count, seq_count) - Number of taxa and total sequences extracted
    """

    if not os.path.exists(bam_path):
        logging.fatal(f"ERROR: BAM not found: {bam_path}")
        return 2

    # Open BAM in main process to validate index and obtain reference lengths
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            if not bam.has_index():
                logging.fatal("ERROR: BAM index (.bai) not found or not readable. Pysam requires an index.")
                return 2
            references = list(bam.references)
    except Exception as e:
        logging.fatal(f"ERROR: Failed to open BAM: {e}")
        return 2

    # Generate tasks for worker processes
    tasks = _iter_tasks(references, max_per_taxon, acc_list, acc_list_action)
    chunk_results = {}
    # Merge results from this chunk
    all_taxon_seqs = {}

    # (default is one-based if neither flag used; argparse sets one_based True by default)
    # endpos will be end0 in both conventions; interpretation differs.

    pool = mp.Pool(
        processes=numthreads,
        initializer=_init_worker,
        initargs=(
            bam_path,
            taxa_dict,
            format,
            matchFraction,
            matchIdentity,
            matchLength,
        ),
    )

    logging.info(f"Starting extraction for qualified taxa...")

    try:
        mapper = pool.imap_unordered
        for result in mapper(_extract_worker, tasks):
            for key, value in result.items():
                chunk_results.setdefault(key, []).extend(value)
    
        for taxid, seqs in chunk_results.items():
            logging.info(f"Extracted {len(seqs)} sequences for taxid {taxid}")

            if taxid not in all_taxon_seqs:
                all_taxon_seqs[taxid] = []

            # Add sequences, respecting the max_per_taxon limit
            if max_per_taxon > 0:
                remaining = max_per_taxon - len(all_taxon_seqs[taxid])

                if remaining > 0:
                    all_taxon_seqs[taxid].extend(seqs[:remaining])
            else:
                all_taxon_seqs[taxid].extend(seqs)
    finally:
        pool.close()
        pool.join()

    # Write sequences to output file
    logging.info("Writing sequences to output file...")

    total_seqs = 0
    taxon_count = 0

    for taxid, seqs in all_taxon_seqs.items():
        if seqs:  # If we got any sequences for this taxon
            taxon_count += 1
            total_seqs += len(seqs)
            
            output_seqs = "".join(seqs)
            o.write(output_seqs)

    return taxon_count, total_seqs