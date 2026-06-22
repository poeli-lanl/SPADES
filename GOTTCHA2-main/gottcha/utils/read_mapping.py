from pathlib import Path
import re
import subprocess
import logging
from typing import List, Tuple
import pandas as pd
import logging
import pandas as pd

def minimap2(reads: List, db: str, threads: int, mm_options: str, presetx: str, samfile: Path, logfile: Path) -> Tuple[int, str, int, bool]:
    """
    Map reads to the reference database using minimap2.

    Builds and executes a command to run minimap2 for read mapping, with parameters
    adjusted based on input settings. Filters the SAM output to keep only relevant
    alignments.

    Parameters:
        reads (List): List of input read file paths
        db (str): Path to the minimap2 index of the reference database
        threads (int): Number of threads to use
        mm_options (str): Minimap2 options for read mapping
        presetx (str): Minimap2 preset mode ('sr', 'map-pb', or 'map-ont')
        samfile (Path): Output SAM file path
        logfile (Path): Log file path
        nanopore (bool): Whether to use Nanopore-specific settings

    Returns:
        Tuple[int, str, int, bool]: (
            exitcode (int): Exit code from the mapping process,
            cmd (str): Command that was executed,
            input_read_count (int): Number of input reads,
            multi_part_index_flag (bool): Flag indicating if a multi-part index was used
        )
    """
    input_file = " ".join(reads)
    mapped_re = re.compile(r"mapped (\d+) sequences")
    multi_part_index_flag = False
    input_read_count = 0

    # Minimap2 options for short reads: the options here is essentailly the -x 'sr' equivalent with some modifications on scoring
    sr_opts = f"-x sr {mm_options} -a -N20 --eqx --secondary=no --sam-hit-only"

    if presetx != 'sr':
        sr_opts = f"-x {presetx} -N20 --secondary=no --sam-hit-only -a"

    mm2_cmd    = f"minimap2 {sr_opts} -t{threads} {db} {input_file}"
    filter_cmd = f"sed '/^@/d'"  # filter out header lines

    # proc = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, universal_newlines=True, bufsize=1)

    with samfile.open("w", encoding="utf-8") as out_f:
        mm2 = subprocess.Popen(
            mm2_cmd,
            shell=True,
            stdout=subprocess.PIPE,      # -> sed
            stderr=subprocess.PIPE,      # <- read THIS (minimap2 only)
            text=True,
            bufsize=1,
        )

        sed = subprocess.Popen(
            filter_cmd,
            shell=True,
            stdin=mm2.stdout,
            stdout=out_f,
            stderr=subprocess.PIPE,      # sed stderr (optional)
            text=True,
            bufsize=1,
        )

        mm2.stdout.close()  # allow mm2 to get SIGPIPE if sed exits

        with logfile.open("a", encoding="utf-8") as f:
            # Stream / parse minimap2 stderr
            for line in mm2.stderr:
                if "For a multi-part index" in line:
                    multi_part_index_flag = True
                
                m = mapped_re.search(line)
                if m:
                    logging.debug(line)
                    input_read_count += int(m.group(1))
                f.write(line)

        mm2.stderr.close()
        sed.stderr.close()

        rc_mm = mm2.wait()
        rc_sed = sed.wait()

    return rc_mm, mm2_cmd, input_read_count, multi_part_index_flag


def post_processing_sam(samfile: Path, samfile_temp: Path) -> Tuple[bool, int, int]:
    """
    Removing multiple hits from the SAM file by keeping only the best alignment for each read.

    Parameters:
        samfile (str): Path to the SAM file
        samfile_temp (str): Path to the temporary SAM file with only the best alignments

    Returns:
        Tuple[bool, int, int]: (
            multiple_hits_removed (bool): False if no multiple hits were found, True if multiple hits were removed,
            total_alignments (int): Total number of alignments in the SAM file,
            top_score_hits (int): Number of top score hits after filtering
        )
    """
    logging.info(f'Loading the SAM file...')

    df = pd.read_csv(samfile,
                     sep='\t',
                     header=None,
                     usecols=[0, 1, 13],
                     names=['QNAME', 'FLAG', 'AS'],
                     converters={
                         'AS': lambda x: x.replace('AS:i:', '')
                     },
                     dtype={'QNAME': 'str', 'FLAG': 'uint16'}
    )

    aln_count = len(df)
    logging.info(f'Total alignments in SAM file: {aln_count}')

    df[['AS']] = df[['AS']].astype('int16')

    logging.info(f'Filtering non-primary hits...')
    # for each row, if the flag bitwise AND with 256 (not primary alignment) or 2048 (supplementary), then remove them from the df
    df = df[~(df['FLAG'] & (256|2048)).astype(bool)]
    logging.info(f'After removing non-primary hits: {len(df)}')

    logging.info(f'Identifying top score hits...')
    # if FLAG bitwise AND with 128 (second in pair), append '/2' to the QNAME
    idx = (df['FLAG'] & 128).astype(bool)
    df.loc[idx, 'QNAME'] = df.loc[idx, 'QNAME'] + '/2'

    # get the index with the best alignment score for each read
    idxmax = df.groupby('QNAME')['AS'].idxmax()
    logging.info(f'Total top score hits: {len(idxmax)}')

    if len(idxmax) == aln_count:
        logging.info(f'No multiple hits found. Keeping the original SAM file.')
        return False, aln_count, aln_count
    else:
        # Create a set of indices for faster lookup
        idxmax_set = set(idxmax.values)
        del idxmax

        logging.info(f'Writing top score hits...')
        with samfile_temp.open("w", encoding="utf-8") as fout, samfile.open("r", encoding="utf-8") as fin:
            lines_to_write = []
            for idx, line in enumerate(fin):
                if idx in idxmax_set:
                    lines_to_write.append(line)
                    if len(lines_to_write) >= 1000:
                        fout.writelines(lines_to_write)
                        lines_to_write.clear()
                        logging.debug(f'Written {idx} alignments...')

            if lines_to_write:
                fout.writelines(lines_to_write)
        logging.info(f'{len(idxmax_set)} hits written to {samfile_temp}.')

        return True, aln_count, len(idxmax_set)
