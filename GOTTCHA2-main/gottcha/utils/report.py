from . import taxonomy
import pandas as pd
import sys

def generate_report_file(rep_df: pd.DataFrame, o: str, fullreport_o: str, fmt: str="tsv") -> bool:
    """
    Generate taxonomy profiling result files in TSV or CSV format.

    Creates two files: a summary file with only qualified results, and
    a full report with all results including filtered entries.

    Parameters:
        rep_df (pandas.DataFrame): Taxonomy results DataFrame
        o (file): Output file handle for the summary results
        fullreport_o (str): Path for the full report file
        fmt (str): Output format, either 'tsv' or 'csv'

    Returns:
        bool: True if successful
    """

    # Fields for full mode
    cols = [# summary 1-11
            'LEVEL', 'NAME', 'TAXID', 'READ_COUNT', 'TOTAL_BP_MAPPED',
            'SNI_SCORE', 'COVERED_SIG_LEN', 'BEST_SIG_COV', 'DEPTH', 'REL_ABUNDANCE_GC',
            'REL_ABUNDANCE',
            # parents 12-13
            'PARENT_NAME', 'PARENT_TAXID',
            # read stats 14-18
            'AOI_READ_COUNT', 'TOTAL_READ_LEN', 'TOTAL_BP_MISMATCH', 'TOTAL_BP_INDEL', 'READ_WT_SNI',
            # signature stats 19-26
            'CONSENSUS_SEQ_SNI', 'SNI_CI95_LH', 'SIG_COV', 'MAPPED_SIG_LEN', 'TOTAL_SIG_LEN', 
            'COVERED_SIG_DEPTH', 'COVERED_MAPPED_SIG_COV', 'ZSCORE',
            # abundance 27-29
            'GENOMIC_CONTENT_EST', 'ABUNDANCE', 'REL_ABUNDANCE_DEPTH',
            # ref genome 30-32
            'SIG_LEVEL', 'GENOME_COUNT', 'GENOME_SIZE', 'NOTE']

    # replace SIG_LEVEL back to their original ranks
    major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7, "strain":8}
    major_ranks = {v:k for k,v in major_ranks.items()}
    rep_df['SIG_LEVEL'] = rep_df['SIG_LEVEL'].map(major_ranks)

    # get qualified taxa
    non_qualified_idx = rep_df['NOTE'].str.contains('Filtered out', na=False) | rep_df['NOTE'].str.contains('Not shown', na=False)
    qualified_df = rep_df.loc[~non_qualified_idx, cols[:11]]  # first 11 columns are summary

    sep = ',' if fmt=='csv' else '\t'

    # save full report
    rep_df[cols].to_csv(fullreport_o, index=False, sep=sep, float_format='%.6f', quoting=2 if fmt=='csv' else 0)

    # save summary
    qualified_df.to_csv(o, index=False, sep=sep, float_format='%.6f', quoting=2 if fmt=='csv' else 0)

    return True


def generate_biom_file(res_df: pd.DataFrame, o: str, tg_rank: str, sampleid: str) -> bool:
    """
    Generate a BIOM format file from taxonomy results.

    Creates a BIOM (Biological Observation Matrix) formatted file for
    compatibility with downstream microbiome analysis tools.

    Parameters:
        res_df (pandas.DataFrame): Taxonomy results DataFrame
        o (file): Output file handle
        tg_rank (str): Target taxonomic rank to include in the output
        sampleid (str): Sample identifier

    Returns:
        bool: True if successful

    Raises:
        SystemExit: If the biom library version is incompatible
    """
    import numpy as np
    import biom
    from biom.table import Table
    if biom.__version__ < '2.1.7':
        sys.exit("[ERROR] Biom library requires v2.1.7 or above.\n")

    target_df = pd.DataFrame()
    target_idx = (res_df['LEVEL']==tg_rank)
    target_df = res_df.loc[target_idx, ['ABUNDANCE','TAXID']]
    target_df['LINEAGE'] = target_df['TAXID'].apply(lambda x: taxonomy.taxid2lineage(x, True, True)).str.split('|')

    sample_ids = [sampleid]
    data = np.array(target_df['ABUNDANCE']).reshape(len(target_df), 1)
    observ_ids = target_df['TAXID']
    observ_metadata = [{'taxonomy': x} for x in target_df['LINEAGE'].tolist()]
    biom_table = Table(data, observ_ids, sample_ids, observ_metadata, table_id='GOTTCHA2')
    biom_table.to_json('GOTTCHA2', direct_io=o)

    return True


def generate_lineage_file(target_df: pd.DataFrame, o: str) -> bool:
    """
    Generate a lineage file showing taxonomic paths with abundances.

    Creates a tab-delimited file with abundance values followed by
    the complete taxonomic lineage for each taxon.

    Parameters:
        target_df (pandas.DataFrame): DataFrame containing abundance and taxids
        o (str): Output file path

    Returns:
        bool: True if successful
    """
    lineage_df = target_df['TAXID'].apply(lambda x: taxonomy.taxid2lineage(x, True, True)).str.split('|', expand=True)
    result = pd.concat([target_df['ABUNDANCE'], lineage_df], axis=1, sort=False)
    result.to_csv(o, index=False, header=False, sep='\t', float_format='%.4f')

    return True


def generate_mpa_file(target_df: pd.DataFrame, o: str) -> bool:
    """
    Generate a lineage file showing taxonomic lineage with abundances in MPA format.

    Creates a tab-delimited file with abundance values followed by
    the complete taxonomic lineage for each taxon.

    Parameters:
        target_df (pandas.DataFrame): DataFrame containing abundance and taxids
        o (str): Output file path

    Returns:
        bool: True if successful
    """

    lineage_df = target_df['TAXID'].apply(lambda x: taxonomy.taxid2lineage(x, all_major_rank=True, print_strain=False, space2underscore=True, sep=";"))
    result = pd.concat([target_df[['TAXID', 'REL_ABUNDANCE', 'REL_ABUNDANCE_GC', 'READ_COUNT', 'SIG_COV']], lineage_df], axis=1, sort=False)
    result.to_csv(o, index=False, header=True, sep='\t', float_format='%.4f')

    return True
