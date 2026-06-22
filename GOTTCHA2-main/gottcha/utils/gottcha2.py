#!/usr/bin/env python3

__author__    = "Po-E (Paul) Li, Bioscience Division, Los Alamos National Laboratory"
__credits__   = ["Po-E Li", "Anna Chernikov", "Jason Gans", "Tracey Freites", "Patrick Chain"]
__version__   = "2.2.3"

import argparse as ap
import sys, os, time, subprocess
import pandas as pd
import numpy as np
import gc
from re import search, findall
from multiprocessing import Pool, set_start_method
from itertools import chain
import math
import logging
from types import SimpleNamespace

try:
    # Try relative import first (for package usage)
    from . import taxonomy
    from . import ont_utils
except ImportError:
    # Fall back to direct import (for script usage)
    import taxonomy
    import ont_utils

def parse_args(ver, args):
    """
    Parse and validate command line arguments for GOTTCHA2.

    This function sets up the argument parser, defines all possible command-line
    options, parses the provided arguments, and performs validation to ensure
    the configuration is valid and complete.

    Parameters:
        ver (str): Version string to display in help messages
        args (list): Command line arguments to parse

    Returns:
        argparse.Namespace: Object containing all validated arguments

    Raises:
        SystemExit: If validation fails or --version is specified
    """
    p = ap.ArgumentParser( prog='gottcha2.py', description="""Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) is an
            annotation-independent and signature-based metagenomic taxonomic profiling tool
            that has significantly smaller FDR than other profiling tools. This program
            is a wrapper to map input reads to pre-computed signature databases using minimap2
            and/or to profile mapped reads in SAM format. (VERSION: %s)""" % ver)

    eg = p.add_mutually_exclusive_group( required=True )

    eg.add_argument( '-i','--input', metavar='[FASTQ]', nargs='+', type=str,
                    help="Input FASTQ/FASTA file(s). Use space to separate multiple input files.")

    eg.add_argument( '-s','--sam', metavar='[SAMFILE]', nargs=1, type=str,
                    help="Specify the input SAM file. Use '-' for standard input. ")

    p.add_argument( '-d','--database', metavar='[GOTTCHA2_db]', type=str, default=None,
                    help="The path and prefix of the GOTTCHA2 database.")

    p.add_argument( '-l','--dbLevel', metavar='[LEVEL]', type=str, default='',
                    choices=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'],
                    help="""Specify the taxonomic level of the input database. You can choose one rank from "superkingdom", "phylum", "class", "order", "family", "genus", "species" and "strain". The value will be auto-detected if the input database ended with levels (e.g. GOTTCHA_db.species).""")

    p.add_argument( '-ti','--taxInfo', metavar='[PATH]', type=str, default='',
                    help="""Specify the path to the taxonomy information directory or file. The program will attempt to locate a matching .tax.tsv file for the specified database. If it cannot find one, it will use the ‘taxonomy_db’ directory located in the same directory as the executable by default.""")

    p.add_argument( '-np','--nanopore', action="store_true",
                    help="Indicate that the input reads are from Oxford Nanopore sequencing platform. This option enables read splitting and error rate set to 0.03 if not specified.")

    p.add_argument( '-e', '--extract', metavar='TAXON[,TAXON2,...]', type=str, default=None,
                    help=(
                        "Extract mapped reads for specific taxa to a FASTA or FASTQ file.\n"
                        "You can specify taxa in one of the following ways:\n"
                        "  - Comma-separated list of taxon IDs:  e.g., -e '1234,5678'\n"
                        "  - File containing a list of taxon IDs (one per line):  e.g., -e '@taxids.txt'\n"
                        "  - File with read limits and format: e.g., -e '@taxids.txt:1000:fasta'\n"
                        "    This limits the number of reads extracted per taxon to <NUMBER> and outputs in <FORMAT> (fasta or fastq).\n"
                        "  Use 'all' to extract all matching taxa/reads.\n"
                        "[default: None]"
                    )
    )

    p.add_argument('-ef', '--extractFullRef', action='store_true',
                    help=(
                        "Extract up to 20 sequences per reference from the SAM file and save them to a FASTA file. "
                        "Equivalent to using: -e 'all:20:fasta'."
                    )
    )

    p.add_argument('-eo', '--extractOnly', action='store_true',
                    help='While --extract is specified, this option will only extract the reads and not perform any further processing of the SAM file.'
    )

    p.add_argument( '-fm','--format', metavar='[STR]', type=str, default='tsv',
                    choices=['tsv','csv','biom'],
                    help='Format of the results; available options include tsv, csv or biom. [default: tsv]')

    p.add_argument( '-r','--relAbu', metavar='[FIELD]', type=str, default='DEPTH',
                    choices=['DEPTH','READ_COUNT','GENOMIC_CONTENT_EST'],
                    help='The field will be used to calculate relative abundance. You can specify one of the following fields: "DEPTH", "READ_COUNT", "GENOMIC_CONTENT_EST". [default: DEPTH]')

    p.add_argument( '-t','--threads', metavar='<INT>', type=int, default=1,
                    help="Number of threads [default: 1]")

    p.add_argument( '-o','--outdir', metavar='[DIR]', type=str, default='.',
                    help="Output directory [default: .]")

    p.add_argument( '-p','--prefix', metavar='<STR>', type=str, required=False,
                    help="Prefix of the output file [default: <INPUT_FILE_PREFIX>]")

    p.add_argument( '-xm','--presetx', metavar='<STR>', type=str, required=False, default='sr',
                    choices=['sr','map-pb','map-ont'],
                    help="The preset option (-x) for minimap2. Default value 'sr' for short reads. [default: sr]")

    p.add_argument( '--m2options', metavar='<STR>', type=str, required=False, default='auto',
                    help="The minimap2 mapping options for short reads. Do not use this option unless you know what you are doing. [default: 'auto']")

    p.add_argument( '-mi','--matchIdentity', metavar='<FLOAT>', type=float,
                    help="Minimum identity (0.0-1.0) required for a valid match. [default: 0.95 for short reads, 0.9 for nanopore reads]")

    p.add_argument( '-mf','--matchFraction', metavar='<FLOAT>', type=float, default=0.99,
                    help="Minimum fraction (0.0-1.0) of the read or signature fragment required to be considered a valid match. [default: 0.99]")

    p.add_argument( '-mg','--matchLength', metavar='<INT>', type=int, default=100,
                    help="Minimum length of the alignment required to be considered a valid match. [default: 100]")

    p.add_argument( '-ss','--sniScore', metavar='<FLOAT>[,<FLOAT>,<FLOAT>]', type=str, default='0.9,0.95,0.99',
                    help="Signature nucleotide identity (SNI) score thresholds for taxonomic aggregation: other levels (first), species level (first value), and strain level (second value); if only one value is provided, all three levels use that value. [default: 0.9,0.95,0.99]")

    p.add_argument( '-Mc','--minCov', metavar='<FLOAT>', type=float, default=0,
                    help="Minimum signature coverage to be considered valid in abundance calculation. [default: 0]")

    p.add_argument( '-Mr','--minReads', metavar='<INT>', type=int, default=0,
                    help="Minimum number of reads to be considered valid in abundance calculation. [default: 0]")

    p.add_argument( '-Ml','--minLen', metavar='<INT>', type=int, default=0,
                    help="Minimum signature length to be considered valid in abundance calculation. [default: 0]")

    p.add_argument( '-Mz','--maxZscore', metavar='<FLOAT>', type=float, default=0,
                    help="Maximum estimated z-score for the depths of the mapped region. Set to 0 to disable. [default: 0]")

    p.add_argument( '-nc','--noCutoff', action="store_true",
                    help="Remove all cutoffs applied during the taxonomic profiling stage (alignment thresholds will remain applied). This option is equivalent to use [-Mc 0 -Mr 0 -Ml 0 -Mz 0 -ss 0,0,0]")

    p.add_argument( '-a','--accList', metavar='[FILE]', required=False, type=str,
                    help="A file of list with accession-of-interest (e.g. plasmid accessions).")

    p.add_argument( '-aa','--accListAction', choices=['filter_out', 'filter_in', 'report_only'], default='report_only', type=str,
                    help=("Action for aligned reads mapping to the accession list. "
                          "'filter_out': discard reads matching accession-of-interest in the list. "
                          "'filter_in': output only reads matching accession-of-interest in the list. "
                          "'report_only': do not filter; report reads matching accession-of-interest in the list (AOI_READ_COUNT). "
                          "[default: report_only]"))

    p.add_argument( '-rm','--removeMultipleHits', choices=['yes', 'no', 'auto'], default='auto', type=str,
                    help="The multiple hit removal step is automatically enabled for sequence input files and disabled for SAM files. Users can explicitly control this behavior by specifying 'yes' or 'no' to force the step to be enabled or disabled. [default: auto]")

    p.add_argument( '-er','--errorRate', metavar='<FLOAT>', type=float,
                    help="Estimated error rate for sequencing data. [default: 0.005]")

    p.add_argument( '-c','--stdout', action="store_true",
                    help="Write on standard output.")

    p.add_argument( '--mpa', action="store_true",
                    help="Generate output in MetaPhlAn format.")

    eg.add_argument( '-v','--version', action="store_true",
                    help="Print version number.")

    p.add_argument( '--silent', action="store_true",
                    help="Disable all messages.")

    p.add_argument( '--verbose', action="store_true",
                    help="Provide verbose messages.")

    p.add_argument( '--debug', action="store_true",
                    help="Debug mode. Provide verbose running messages and keep all temporary files.")

    args_parsed = p.parse_args(args)

    """
    Checking options
    """
    if args_parsed.version:
        print( ver )
        sys.exit(0)

    if args_parsed.extract and args_parsed.extractFullRef:
        p.error( '--extract and --extractFullRef are incompatible options.' )

    if not args_parsed.database:
        p.error( '--database option is missing.' )

    if args_parsed.input and args_parsed.sam:
        p.error( '--input and --same are incompatible options.' )

    if args_parsed.database:
        #assign default path for database name
        if "/" not in args_parsed.database and not os.path.isfile( args_parsed.database + ".mmi" ):
            bin_dir = os.path.dirname(os.path.realpath(__file__))
            args_parsed.database = bin_dir + "/database/" + args_parsed.database

    if args_parsed.database and args_parsed.database.endswith(".mmi"):
        args_parsed.database.replace('.mmi','')

    if args_parsed.database and args_parsed.input:
        if not os.path.isfile( args_parsed.database + ".mmi" ):
            p.error( 'Database index %s.mmi not found.' % args_parsed.database )

    if args_parsed.input:
        validated_inputs = []
        for path in args_parsed.input:
            if path == '-':
                p.error('--input does not support reading from stdin ("-"). Please provide a file path.')
            if not os.path.isfile(path):
                p.error(f'Input file {path} not found.')
            validated_inputs.append(SimpleNamespace(name=os.path.abspath(path)))
        args_parsed.input = validated_inputs

    if args_parsed.sam:
        sam_path = args_parsed.sam[0]
        if sam_path != '-' and not os.path.isfile(sam_path):
            p.error(f'SAM file {sam_path} not found.')
        sam_path_name = sam_path if sam_path == '-' else os.path.abspath(sam_path)
        args_parsed.sam = [SimpleNamespace(name=sam_path_name)]

    if args_parsed.accList:
        if not os.path.isfile(args_parsed.accList):
            p.error(f'Accession exclusion list {args_parsed.accList} not found.')
        args_parsed.accList = os.path.abspath(args_parsed.accList)

    if args_parsed.nanopore and args_parsed.input and len(args_parsed.input) != 1:
        p.error( '--nanopore option requires a single input read file.' )

    if not args_parsed.prefix:
        if args_parsed.input:
            name = search(r'([^\/\.]+)\..*$', args_parsed.input[0].name )
            args_parsed.prefix = name.group(1)
        elif args_parsed.sam:
            name = search(r'([^\/]+).\w+.\w+$', args_parsed.sam[0].name )
            args_parsed.prefix = name.group(1)
        else:
            args_parsed.prefix = "GOTTCHA_"

    if not args_parsed.dbLevel:
        if args_parsed.database:
            major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7, "strain":8}
            parts = args_parsed.database.split('.')
            for part in parts:
                if part in major_ranks:
                    args_parsed.dbLevel = part
                    break
        elif args_parsed.sam:
            name = search(r'\.gottcha_(\w+).sam$', args_parsed.sam[0].name )
            try:
                args_parsed.dbLevel = name.group(1)
            except:
                pass

        if not args_parsed.dbLevel:
            p.error( '--dbLevel is missing and cannot be auto-detected.' )

    if args_parsed.removeMultipleHits == 'auto':
        if args_parsed.input:
            args_parsed.removeMultipleHits = "yes"
        else:
            args_parsed.removeMultipleHits = "no"

    if args_parsed.extractFullRef:
        args_parsed.extract = 'all:20:fasta'

    if args_parsed.noCutoff:
        args_parsed.sniScore = '0,0,0'

    if args_parsed.m2options == 'auto':
        args_parsed.m2options = '-s60'

    if args_parsed.matchIdentity:
        if args_parsed.matchIdentity < 0 or args_parsed.matchIdentity > 1:
            p.error( '--matchIdentity must be between 0 and 1.' )

    if args_parsed.matchIdentity is None:
         if args_parsed.nanopore:
            args_parsed.matchIdentity = 0.9
            args_parsed.matchFraction = 0.9
         else:
            args_parsed.matchIdentity = 0.95

    if not args_parsed.errorRate:
        if args_parsed.nanopore:
            args_parsed.errorRate = 0.03
        else:
            args_parsed.errorRate = 0.005

    return args_parsed

def dependency_check(cmd):
    """
    Verify that external dependencies are available in the system.

    Attempts to execute the specified command with --help and checks if it runs
    successfully. Exits the program if the command is not found or fails.

    Parameters:
        cmd (str): Command to check

    Returns:
        None

    Raises:
        SystemExit: If the command is not found or fails
    """
    try:
        subprocess.check_call([cmd, "--help"], stdout=subprocess.DEVNULL)
    except Exception as e:
        sys.stderr.write(f"[ERROR] {cmd}: {e}\n")
        sys.exit(1)

def merge_ranges(ranges):
    """
    Merge overlapping or consecutive genomic ranges.

    Takes a list of (start, end) tuples representing genomic ranges and merges
    any ranges that overlap or are directly adjacent (consecutive positions).
    This is used to calculate accurate linear coverage for mapped reads.

    Parameters:
        ranges (list): List of (start, end) tuples representing genomic ranges

    Returns:
        list: List of merged (start, end) tuples with no overlaps

    Example:
        >>> merge_ranges([(1, 5), (4, 8), (10, 12)])
        [(1, 8), (10, 12)]
    """
    # Sort ranges by start position
    sorted_ranges = sorted(ranges, key=lambda x: x[0])
    merged = []
    for current in sorted_ranges:
        if not merged:
            merged.append(current)
        else:
            last = merged[-1]
            # Check if current range overlaps or is adjacent to the last range
            if current[0] <= last[1] + 1:
                # Merge the ranges
                merged[-1] = (last[0], max(last[1], current[1]))
            else:
                merged.append(current)
    return merged

def worker(filename, chunkStart, chunkSize, matchFraction, matchIdentity, matchLength, split_read_flag=False):
    """
    Process a chunk of a SAM file to extract mapping information.

    This function is intended to be run in parallel to process different chunks
    of a SAM file. It parses lines within the specified chunk and builds a dictionary
    of reference sequences with their mapped regions, base counts, read counts, etc.

    Parameters:
        filename (str): Path to the SAM file to process
        chunkStart (int): Byte position in the file where to start reading
        chunkSize (int): Number of bytes to read from the start position
        matchFraction (float): Minimum fraction required for a valid match
        matchIdentity (float): Minimum identity required for a valid match
        split_read_flag (bool): Flag indicating whether reads are splitted (optional)

    Returns:
        dict: Dictionary with reference sequences as keys and mapping statistics as values
        int: Number of lines processed in this chunk
        int: Number of invalid matches found
        int: Number of accession#s of interest found
    """
    # processing alignments in SAM format
    f = open( filename )
    f.seek(chunkStart)
    lines = f.read(chunkSize).splitlines()
    res={}
    invalid_match_count=0

    for line in lines:
        k, r, nm, nid, rd, rs, rq, flag, cigr, read_len, pri_aln_flag, valid_match_flag, sr_chunk_flag = parse(line, matchFraction, matchIdentity, matchLength)
        # parsed values from SAM line
        # only k, r, n, pri_aln_flag, valid_flag are used
        # k: reference name
        # r: (start, end) of mapped region
        # nm: number of mismatches
        # nid: number of indels
        # read_len: read length
        # pri_aln_flag: whether this is a primary alignment
        # valid_match_flag: whether this alignment meets match criteria
        # valid_acc_flag: whether this alignment maps to a valid accession
        # sr_chunk_flag: whether this is a chunked read

        if not valid_match_flag:
            invalid_match_count += 1
            continue

        if pri_aln_flag:
            if k in res:
                res[k]['REGIONS'] = merge_ranges(res[k]['REGIONS']+[r])
                res[k]["MB"] += r[1] - r[0] + 1
                res[k]["MR"] += 0 if split_read_flag and (not sr_chunk_flag) else 1
                res[k]["NM"] += nm
                res[k]["ID"] += nid
                res[k]["RL"] += read_len
            else:
                res[k]={}
                res[k]["REGIONS"] = [r]
                res[k]["MB"] = r[1] - r[0] + 1
                res[k]["MR"] = 0 if split_read_flag and (not sr_chunk_flag) else 1
                res[k]["NM"] = nm
                res[k]["ID"] = nid
                res[k]["RL"] = read_len

    return (res, len(lines), invalid_match_count)

def parse(line, matchFraction, matchIdentity, matchLength):
    """
    Parse a line from a SAM file and extract relevant mapping information.

    Parses alignment details from a SAM format line, including reference ID,
    match position, mismatches, sequence quality, and flags. Determines if
    the alignment is a valid match based on matchFraction criteria.

    Parameters:
        line (str): A line from a SAM file
        matchFraction (float): Minimum fraction required for a valid match
        matchIdentity (float): Minimum identity required for a valid match

    Returns:
        tuple: (
            ref (str): Reference identifier,
            (start, end) (tuple): Mapped region coordinates,
            mismatches (int): Number of mismatches,
            read_name (str): Name of the read,
            read_seq (str): Read sequence,
            read_qual (str): Read quality string,
            flag (str): SAM flag,
            cigar (str): CIGAR string,
            primary_alignment_flag (bool): Whether this is a primary alignment,
            valid_match_flag (bool): Whether this alignment meets match criteria,
            sr_chunk_flag (bool): Whether this is a chunked read
        )

    Example:
        SAM format example:
        read1   0   ABC|1|100|GCF_12345|    11  0   5S10M3S *   0   0   GGGGGCCCCCCCCCCGGG  HHHHHHHHHHHHHHHHHH  NM:i:0  MD:Z:10 AS:i:10 XS:i:0
        read2   16  ABC|1|100|GCF_12345|    11  0   3S10M5S *   0   0   GGGCCCCCCCCCCGGGGG  HHHHHHHHHHHHHHHHHH  NM:i:0  MD:Z:10 AS:i:10 XS:i:0
    """
    temp = line.split('\t')
    name = temp[0]
    cigr = temp[5]
    mapped_len = 0

    # parse NM tag for mismatch length
    mismatch_len = search(r'NM:i:(\d+)', line)
    mismatch_len = int(mismatch_len.group(1)) if mismatch_len else 0

    # parse the CIGAR string for match length, search all r'(\d+)M' and sum all matches
    if "M" in cigr:
        mapped_len = sum(int(num) for num in findall(r'(\d+)M', cigr))
    else:
        mapped_len = sum(int(num) for num in findall(r'(\d+)=', cigr))
        mapped_len += mismatch_len

    # count for deletion and insertion length
    ins_len = 0
    del_len = 0
    if 'I' in cigr:
        ins_len = sum(int(num) for num in findall(r'(\d+)I', cigr)) if 'I' in cigr else 0
    if 'D' in cigr:
        del_len = sum(int(num) for num in findall(r'(\d+)D', cigr)) if 'D' in cigr else 0
    indel_len = ins_len + del_len

    start = int(temp[3])
    end   = start + mapped_len + del_len - 1
    read_len = len(temp[9])
    # determine if this is a first seen chunked read
    sr_chunk_flag = True if line.endswith('ZC:i:1') else False

    ref = temp[2].rstrip('|')
    ref = ref[: -2 if ref.endswith(".0") else None ]

    (acc, rstart, rend, taxid) = ref.split('|')

    # this is a workaround for the case when the reference name is not in the format of "accession|start|end|taxid"
    rlen = int(rend)-int(rstart)+1

    # check if this is a primary alignment 256=secondary, 2048=supplementary
    primary_alignment_flag=False if int(temp[1]) & 2304 else True

    # the alignment region should cover at least matchFraction(proportion) of the read or signature fragment
    valid_match_flag = True

    # get the max matching proportion of the read or the signature fragment
    match_prop = max(mapped_len/rlen, mapped_len/read_len)
    # get the match identity of the mapped region
    match_idt = ((mapped_len-(indel_len + mismatch_len)) / mapped_len) if mapped_len > 0 else 0

    if (match_prop >= matchFraction) and (match_idt >= matchIdentity) and (mapped_len >= matchLength):
        valid_match_flag = True
    else:
        valid_match_flag = False

    return (ref,
            (start, end),
            mismatch_len,
            indel_len,
            name,
            temp[9],
            temp[10],
            temp[1],
            temp[5],
            read_len,
            primary_alignment_flag,
            valid_match_flag,
            sr_chunk_flag)

def time_spend(start):
    """
    Calculate and format elapsed time since a given start time.

    Parameters:
        start (float): Starting time in seconds (as returned by time.time())

    Returns:
        str: Formatted time string in HH:MM:SS format
    """
    done = time.time()
    elapsed = done - start
    return time.strftime( "%H:%M:%S", time.gmtime(elapsed) )

def chunkify(fname, size=1*1024*1024):
    """
    Split a file into chunks for parallel processing.

    Divides a file into chunks of approximately the specified size, ensuring
    that all alignments for a single read are kept in the same chunk.
    This is critical for accurate processing of multi-mapped reads.

    Parameters:
        fname (str): Path to the file to be chunked
        size (int): Approximate chunk size in bytes (default: 1MB)

    Yields:
        tuple: (chunkStart, chunkSize) where:
            - chunkStart (int): Byte position to start reading
            - chunkSize (int): Number of bytes to read from that position
    """
    fileEnd = os.path.getsize(fname)
    with open(fname, "rb") as f:
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(size, 1)
            f.readline()
            # put all alignments of a read in the same chunck
            line = f.readline().decode('ascii')
            tmp = line.split('\t')
            if chunkEnd <= fileEnd and line:
                last = f.tell()
                line = f.readline().decode('ascii')
                while line.startswith(tmp[0]):
                    last = f.tell()
                    line = f.readline().decode('ascii')
                f.seek(last)
            # current position
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                break

def process_sam_file(sam_fn, numthreads, matchFraction, matchIdentity, matchLength, split_read_flag=False):
    """
    Process a SAM file using parallel execution to extract mapping information.

    Divides the SAM file into chunks, processes each chunk in parallel using a thread pool,
    and then merges the results. Computes the linear coverage for each reference sequence.

    Parameters:
        sam_fn (str): Path to the SAM file
        numthreads (int): Number of parallel processes to use
        matchFraction (float): Minimum fraction required for a valid match
        matchIdentity (float): Minimum identity required for a valid match
        split_read_flag (bool): Flag indicating whether to split reads (optional)

    Returns:
        tuple: (
            result (dict): Dictionary with references as keys and mapping statistics as values,
            mapped_reads (int): Total number of reads that mapped
            tol_alignment_count (int): Total number of alignments processed
            tol_invalid_match_count (int): Total number of invalid matches found
            tol_reportable_acc_count (int): Total number of accession#s of interest found
        )
    """
    result = taxonomy._autoVivification()
    mapped_reads = 0

    print_message(f" - Processing with {numthreads} cpus...", argvs.silent, begin_t, logfile)
    pool = Pool(processes=numthreads)
    jobs = []
    results = []
    tol_invalid_match_count = 0
    tol_alignment_count = 0

    for chunkStart,chunkSize in chunkify(sam_fn):
        _args = (sam_fn,
                 chunkStart,
                 chunkSize,
                 matchFraction,
                 matchIdentity,
                 matchLength,
                 split_read_flag)
        jobs.append( pool.apply_async(worker, _args) )

    #wait for all jobs to finish
    tol_jobs = len(jobs)
    cnt=0
    for job in jobs:
        results.append( job.get() )
        cnt+=1
        if argvs.debug:
            logging.debug( f"[DEBUG] Progress: {cnt}/{tol_jobs} ({cnt/tol_jobs*100:.1f}) chunks done.")

    #clean up
    pool.close()

    print_message(f" - Merging {tol_jobs} jobs...", argvs.silent, begin_t, logfile)
    for res_tuples in results:
        (res, alignment_count, invalid_match_count) = res_tuples
        tol_alignment_count += alignment_count
        tol_invalid_match_count += invalid_match_count

        for k in res:
            if k in result:
                result[k]['REGIONS'] = merge_ranges(result[k]['REGIONS']+res[k]['REGIONS'])
                result[k]["MB"] += res[k]["MB"]
                result[k]["MR"] += res[k]["MR"]
                result[k]["NM"] += res[k]["NM"]
                result[k]["ID"] += res[k]["ID"]
                result[k]["RL"] += res[k]["RL"]
            else:
                result[k]={}
                result[k].update(res[k])

    # convert mapped regions to covered signature length
    refs = result.keys()
    for k in list(refs):
        if not result[k]["MR"]:
            del result[k]
        else:
            result[k]["SC"] = sum(end - start + 1 for start, end in result[k]['REGIONS'])
            del result[k]['REGIONS']
            mapped_reads += result[k]["MR"]

    for k in result:
        logging.debug(f'Processed {k}: {result[k]["MR"]} reads, {result[k]["SC"]} covbases, {result[k]["NM"]} mismatches, {result[k]["MB"]} mapped_bases, {result[k]["ID"]} indels, {result[k]["RL"]} read_length')

    return result, mapped_reads, tol_alignment_count, tol_invalid_match_count

def extract_sequences_by_taxonomy(sam_fn,
                                  taxa_dict,
                                  qualified_taxids,
                                  o,
                                  numthreads,
                                  matchFraction,
                                  matchIdentity,
                                  matchLength,
                                  max_per_taxon,
                                  acc_list,
                                  acc_list_action,
                                  format='fasta'):
    """
    Extract sequences mapping to taxa from the full taxonomy report.

    For each taxon in the full report, extract up to max_per_taxon sequences.

    Parameters:
        sam_fn (str): Path to the SAM file
        full_tsv_fn (str): Path to the full taxonomy report file
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

    # Calculate total file size for progress reporting
    file_size = os.path.getsize(sam_fn)
    chunk_positions = list(chunkify(sam_fn))
    total_chunks = len(chunk_positions)

    print_message(f"Processing SAM file ({file_size/1024/1024:.1f} MB) in {total_chunks} chunks with {numthreads} threads",
                argvs.silent, begin_t, logfile)

    # Process in batches to show progress
    batch_size = max(1, total_chunks // 10)  # Show progress in ~10% increments

    # Process the SAM file in batches
    all_taxon_seqs = {}
    processed_chunks = 0

    for i in range(0, total_chunks, batch_size):
        batch_chunks = chunk_positions[i:i+batch_size]

        pool = Pool(processes=numthreads)
        jobs = []

        # Submit jobs for this batch
        for chunkStart, chunkSize in batch_chunks:
            jobs.append(pool.apply_async(
                OptimizedFastaWorker,
                (sam_fn, chunkStart, chunkSize, taxa_dict, qualified_taxids, matchFraction, matchIdentity, matchLength, max_per_taxon, acc_list, acc_list_action, format)
            ))

        # Process results as they complete
        for job in jobs:
            chunk_results = job.get()
            processed_chunks += 1

            # Report progress
            progress = processed_chunks / total_chunks * 100
            print_message(f"Progress: {progress:.1f}% ({processed_chunks}/{total_chunks} chunks)",
                        argvs.silent, begin_t, logfile)

            # Merge results from this chunk
            for taxid, seqs in chunk_results.items():
                logging.info(f"Processing {len(seqs)} sequences for taxid {taxid}")

                if taxid not in all_taxon_seqs:
                    all_taxon_seqs[taxid] = []

                # Add sequences, respecting the max_per_taxon limit
                if max_per_taxon > 0:
                    remaining = max_per_taxon - len(all_taxon_seqs[taxid])

                    if remaining > 0:
                        all_taxon_seqs[taxid].extend(seqs[:remaining])
                else:
                    all_taxon_seqs[taxid].extend(seqs)

        # Clean up this batch's pool
        pool.close()
        pool.join()

    # Write sequences to output file
    print_message("Writing sequences to output file...",
                argvs.silent, begin_t, logfile)

    total_seqs = 0
    taxon_count = 0

    for taxid, seqs in all_taxon_seqs.items():
        if seqs:  # If we got any sequences for this taxon
            taxon_count += 1
            for seq in seqs:  # Write up to max_per_taxon
                o.write(seq)
                total_seqs += 1

    return taxon_count, total_seqs

def OptimizedFastaWorker(filename, chunkStart, chunkSize, taxa_dict, qualified_taxids, matchFraction, matchIdentity, matchLength, max_per_taxon, acc_list, acc_list_action, format):
    """
    Worker function that processes a chunk of the SAM file and extracts sequences for all taxa.

    Parameters:
        filename (str): Path to the SAM file
        chunkStart (int): Starting position in the file
        chunkSize (int): Size of the chunk to process
        taxa_dict (dict): Dictionary mapping taxids to their level and name
        qualified_taxids (list): List of taxids to check against
        matchFraction (float): Minimum fraction required for a valid match
        matchIdentity (float): Minimum identity required for a valid match
        max_per_taxon (int): Maximum sequences to extract per taxon
        acc_list (list, optional): List of accession-of-interest
        acc_list_action (str, optional): Action to take with the accession list (e.g., "filter_out")
        format (str): Output format ('fasta' or 'fastq')

    Returns:
        dict: Dictionary mapping taxids to lists of their FASTA sequences
    """
    taxon_seqs = {}  # Dictionary to hold sequences for each taxid
    processed_lines = 0

    # Process the SAM file chunk
    f = open(filename)
    f.seek(chunkStart)
    lines = f.read(chunkSize).splitlines()

    # Create a more efficient lookup for cached lineages
    lineage_cache = {}

    for line in lines:
        processed_lines += 1

        try:
            ref, region, nm, nid, rname, rseq, rq, flag, cigr, read_len, pri_aln_flag, valid_match_flag, sr_chunk_flag = parse(line, matchFraction, matchIdentity, matchLength)

            if not (pri_aln_flag and valid_match_flag):
                continue

            # Extract taxid from reference
            try:
                acc, rstart, rend, ref_taxid = ref.split('|')
            except ValueError:
                logging.debug(f"Malformed reference: {ref}")
                continue  # Skip malformed references

            # Skip if accession is in the exclusion list (if applicable)
            aoi_flag = False
            if acc_list:
                if acc in acc_list:
                    aoi_flag = True
                    if acc_list_action == 'filter_out':
                        continue
                else:
                    if acc_list_action == 'filter_in':
                        continue

            # Check if we already know what qualified taxa this reference belongs to
            if ref_taxid in lineage_cache:
                matching_taxids = lineage_cache[ref_taxid]
            else:
                # If not, find all qualified taxa this reference belongs to
                matching_taxids = []
                ref_lineage = None

                for q_taxid in qualified_taxids:
                    # Avoid recomputing the lineage for each taxid check
                    if ref_lineage is None:
                        ref_lineage = taxonomy.taxid2fullLineage(ref_taxid, space2underscore=False)

                    if f"|{q_taxid}|" in ref_lineage:
                        matching_taxids.append(q_taxid)

                # Cache the result
                lineage_cache[ref_taxid] = matching_taxids

            # Process the matching taxa
            for taxid in matching_taxids:
                # Handle reverse complement if needed
                rc_seq = None
                if int(flag) & 16:
                    if rc_seq is None:  # Only compute once if needed
                        g = findall(r'\d+\w', cigr)
                        cigr = "".join(list(reversed(g)))
                        rc_seq = seqReverseComplement(rseq)
                        rq = rq[::-1]
                    seq_to_use = rc_seq
                else:
                    seq_to_use = rseq

                # Initialize list for this taxid if needed
                if taxid not in taxon_seqs:
                    taxon_seqs[taxid] = []

                # Only collect up to max_per_taxon sequences per taxon
                if (max_per_taxon==0) or (len(taxon_seqs[taxid]) < max_per_taxon):
                    # Create FASTA entry with taxonomy information
                    level = taxa_dict[taxid]['level']
                    name = taxa_dict[taxid]['name']

                    # determine if the read is the first or second mate
                    mate = ''
                    if (int(flag) & 64) | (int(flag) & 128):
                        if int(flag) & 64:
                            mate = '.1'
                        elif int(flag) & 128:
                            mate = '.2'

                    mapping_len = region[1] - region[0] + 1
                    ref_len = int(rend) - int(rstart) + 1
                    mapping_idt = 1 - (nm + nid) / mapping_len
                    mapping_frac = max(mapping_len/read_len, mapping_len/ref_len)

                    if format == 'fasta':
                        fasta_entry = f">{rname}{mate}|{ref}:{region[0]}..{region[1]} LEVEL={level} NAME={name} TAXID={taxid} AOI={aoi_flag} MG={mapping_len} MI={mapping_idt:.2f} MF={mapping_frac:.2f}\n{seq_to_use}\n"
                    else:
                        fasta_entry = f"@{rname}{mate}|{ref}:{region[0]}..{region[1]} LEVEL={level} NAME={name} TAXID={taxid} AOI={aoi_flag} MG={mapping_len} MI={mapping_idt:.2f} MF={mapping_frac:.2f}\n{seq_to_use}\n+\n{rq}\n"
                    taxon_seqs[taxid].append(fasta_entry)
        except Exception as e:
            # Skip problematic lines
            logging.debug(f"Error processing line: {line}\n{e}")
            continue

    return taxon_seqs

def seqReverseComplement(seq):
    """
    Generate the reverse complement of a DNA sequence.

    Creates a mapping dictionary for complementary bases and applies it to the
    reversed sequence. Handles both uppercase and lowercase nucleotides.

    Parameters:
        seq (str): DNA sequence string

    Returns:
        str: Reverse complemented DNA sequence
    """
    seq1 = 'ACGTURYSWKMBDHVNTGCAAYRSWMKVHDBNacgturyswkmbdhvntgcaayrswmkvhdbn'
    seq_dict = { seq1[i]:seq1[i+16] for i in range(64) if i < 16 or 32<=i<48 }
    return "".join([seq_dict[base] for base in reversed(seq)])

def group_refs_to_strains(r, acc_list, acc_list_action):
    """
    Group reference mapping results by strains and calculate strain-level statistics.

    Converts the mapping results dictionary to a pandas DataFrame and groups by
    taxonomic identifier. Calculates various statistics including total mapped bases,
    read counts, coverage, and depth of coverage.

    Parameters:
        r (dict): Dictionary with reference sequences as keys and mapping statistics
                 as values (output from process_sam_file)
        acc_list (list, optional): List of accession-of-interest
        acc_list_action (str, optional): Action to take with the accession list (e.g., "report_only")
    Returns:
        pandas.DataFrame: DataFrame with strain-level statistics
    """
    # covert mapping info to df
    r_df = pd.DataFrame.from_dict(r, orient='index').reset_index()
    r_df.rename(columns={"index": "RNAME"}, inplace=True)
    # retrieve sig fragment info
    r_df['RNAME'] = r_df['RNAME'].str.rstrip('|')
    r_df[['ACC','RSTART','REND','TAXID']] = r_df['RNAME'].str.split('|', expand=True)

    # add AOI read count
    r_df['RR'] = 0
    aoi_read_count = 0

    if acc_list:
        idx = (r_df['ACC'].isin(acc_list) | r_df['RNAME'].isin(acc_list))
        r_df.loc[idx, 'RR'] = r_df.loc[idx, 'MR'] # report the read count for the accession#s of interest
        aoi_read_count = r_df.loc[idx, 'MR'].sum()

        if acc_list_action == 'filter_out':
            r_df = r_df.loc[~idx] # set mapped bases, read count, mismatch and covered sig len to 0 for the accession#s of interest
        elif acc_list_action == 'filter_in':
            r_df = r_df[idx].reset_index(drop=True)

        # if after applying the accession list filter, there is no valid mapping left, exit the program
        if len(r_df) == 0:
            logging.info(f"No valid mappings after applying accession list filter. Exiting.")
            sys.exit(0)

    r_df['RSTART'] = r_df['RSTART'].astype(int)
    r_df['REND'] = r_df['REND'].astype(int)
    r_df['SLEN'] = r_df['REND']-r_df['RSTART']+1 # length of the signature fragment

    # group by strain
    str_df = r_df.groupby(['TAXID']).agg({
        'MB':'sum', # of mapped bases
        'MR':'sum', # of mapped reads
        'NM':'sum', # of mismatches
        'ID':'sum', # of indels
        'SC':'sum', # covered signature length
        'SLEN':'sum', # length of this signature fragments (mapped)
        'RL':'sum', # length of the reads
        'RR':'sum'  # reportable read count
    }).reset_index()
    # total length of signatures
    str_df['TS'] = str_df['TAXID'].map(df_stats['TotalLength'])
    str_df['bLC'] = str_df['SC']/str_df['TS'] # bLC:  best linear coverage of a strain
    str_df['RD'] = str_df['MB']/str_df['TS'] # roll-up DoC
    str_df['NOTE'] = str_df['TAXID'].map(df_stats['Note']).fillna('') # note for the strain

    # rename columns
    str_df.rename(columns={
        "MB":   "TOTAL_BP_MAPPED",
        "MR":   "READ_COUNT",
        "NM":   "TOTAL_BP_MISMATCH",
        "ID":   "TOTAL_BP_INDEL",
        "RL":   "TOTAL_READ_LEN",
        "SC":   "COVERED_SIG_LEN",
        "SLEN": "MAPPED_SIG_LEN", # length of the mapped signature fragments (entire fragment)
        "TS":   "TOTAL_SIG_LEN",
        "RD":   "DEPTH",
        "bLC":  "BEST_SIG_COV",
        "RR":   "AOI_READ_COUNT"
    }, inplace=True)

    # check if TOTAL_SIG_LEN is 0, report the TAXID and exit
    # this should not happen if the database and corresponding stats file are correct
    if str_df['TOTAL_SIG_LEN'].eq(0).any():
        logging.fatal(f"Error: total signature length is ZERO for some mapped strains. Please check your database.")
        sys.exit(1)

    # get genome size
    str_df['SIG_LEVEL'] = str_df['TAXID'].map(df_stats['DB_level'])
    str_df['GENOME_SIZE'] = str_df['TAXID'].map(df_stats['GenomeSize'])
    str_df['GENOME_COUNT'] = 1

    # infer total genome contents
    str_df['GENOMIC_CONTENT_EST'] = str_df['TOTAL_BP_MAPPED']/str_df['TOTAL_SIG_LEN']*str_df['GENOME_SIZE']

    # estimate z-score
    str_df['ZSCORE'] = str_df.apply(lambda x: pile_lvl_zscore(x.TOTAL_BP_MAPPED, x.TOTAL_SIG_LEN, x.COVERED_SIG_LEN), axis=1)

    return str_df, aoi_read_count

def aggregate_taxonomy(str_df, abu_col, tg_rank, mc, mr, ml, mz, sni_score_species, sni_score_strain, sni_score_cutoff, error_rate):
    """
    Aggregate strain-level results to higher taxonomic ranks.

    Starting from strain-level mapping data, this function rolls up statistics to
    higher taxonomic ranks (species, genus, family, etc.). It applies the specified
    cutoff criteria to filter results and marks entries that fall below these thresholds.

    The process of aggregating taxonomic data is done in a bottom-up manner, starting from the strain level and moving up to the superkingdom level.
        1. First identify the taxon name and taxid at each major rank for each strain.
        2. Then, identify the strains that meet the cutoff criteria.
        3. For each rank starting from species, aggregate the qualify strains by summing up the relevant statistics (e.g., total mapped bases, read counts, etc.) to each rank.
        4. For unqualified strains, a note is added to indicate the reason for exclusion.
        3. Finally, the aggregated data is stored in a DataFrame(rep_df), which is returned as the output of the function.

    Parameters:
        r (dict): Dictionary with reference sequences as keys and mapping stats as values
        abu_col (str): Column name to use for abundance calculations
        tg_rank (str): Target taxonomic rank
        mc (float): Minimum linear coverage threshold
        mr (int): Minimum read count threshold
        ml (int): Minimum covered signature length threshold
        mz (float): Maximum Z-score threshold (0 to disable)
        sni_score_cutoff (float): SNI-score cutoff for all levels
        sni_score_species (float): SNI-score cutoff for species level
        sni_score_strain (float): SNI-score cutoff for strain level
        error_rate (float): Error rate for SNI-score inference
        acc_list (list, optional): List of accessions to filter
        acc_list_action (str, optional): Action to take with the accession list (e.g., "report_only")

    Returns:
        pandas.DataFrame: DataFrame with rolled-up taxonomy at all ranks
    """

    major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7,"strain":8}

    # produce columns for the final report at each ranks
    rep_df = pd.DataFrame()

    # add taxonomic lineage info
    ranks = list(major_ranks.keys())[::-1]

    def get_taxid_lineage(taxid):
        """get taxid lineage with {rank}_names and {rank}_taxids"""
        lineage = taxonomy.taxid2lineageDICT(taxid).values()
        return [d['name'] for d in lineage]+[d['taxid'] for d in lineage]

    try:
        cols = [f'{r}_name' for r in ranks]+[f'{r}_taxid' for r in ranks]
        str_df[cols] = str_df['TAXID'].map(get_taxid_lineage).to_list()
    except Exception as e:
        logging.error(f"Error processing rank {ranks}: {e}. Please verify that your taxonomy file matches the expected database.")
        sys.exit(1)

    # decide top signature level, convert the rank to the corresponding number
    str_df['SIG_LEVEL'] = str_df['SIG_LEVEL'].map(major_ranks)

    # infer the SNI-score for each strain
    str_df["SIG_COV"] = str_df["COVERED_SIG_LEN"]/str_df["TOTAL_SIG_LEN"]
    str_df = infer_sni_score(str_df, error_rate)

    total_abundance = str_df[abu_col].sum()

    # iterate through ranks to get index and value
    for idx, rank in enumerate(ranks):
        str_df['LEVEL'] = rank
        str_df[['LVL_NAME', 'LVL_TAXID']] = str_df[[f'{rank}_name', f'{rank}_taxid']]

        if rank=='superkingdom':
            str_df[['PARENT_NAME', 'PARENT_TAXID']] = ['root', '1']
        else:
            str_df[['PARENT_NAME', 'PARENT_TAXID']] = str_df[[f'{ranks[idx+1]}_name', f'{ranks[idx+1]}_taxid']]

        # rollup strains that make cutoffs
        lvl_df = None
        if rank == 'strain':
            lvl_df = str_df.copy()
        else:
            lvl_df = str_df.groupby('LVL_NAME').agg({
                'LEVEL':'first',
                'LVL_TAXID':'first',
                'PARENT_NAME':'first',
                'PARENT_TAXID':'first',
                'TOTAL_BP_MAPPED': 'sum',
                'READ_COUNT': 'sum',
                'TOTAL_BP_MISMATCH': 'sum',
                'TOTAL_BP_INDEL': 'sum',
                'TOTAL_READ_LEN': 'sum',
                'COVERED_SIG_LEN': 'sum',
                'MAPPED_SIG_LEN': 'sum',
                'TOTAL_SIG_LEN': 'sum',
                'DEPTH': 'sum',
                'AOI_READ_COUNT': 'sum',
                'BEST_SIG_COV': 'max',
                'ZSCORE': 'min',
                'GENOMIC_CONTENT_EST': 'sum',
                'SIG_LEVEL': 'max',
                'GENOME_COUNT': 'count',
                'GENOME_SIZE': 'sum',
                'SNI_SCORE': 'max',
                'NOTE': lambda x: '; '.join(list(x.unique()))
            })

            # find the index of the row with max SNI_SCORE in each group
            # pull out the low/high bounds from those rows
            idx = str_df.groupby('LVL_NAME')['SNI_SCORE'].idxmax()
            score_bounds = ( str_df
                        .loc[idx, ['LVL_NAME', 'SNI_NAIVE', 'SNI_CI95_LH']]
                        .set_index('LVL_NAME') )
            lvl_df = lvl_df.join(score_bounds).reset_index()

        # calculate the relative abundance of each taxon
        # the abundance and the relative abundance is calculated based on the specified column (abu_col)
        # Other depth-based and adjusted-genomic-content-based abundance values are also included
        lvl_df['ABUNDANCE'] = lvl_df[abu_col]
        lvl_df['REL_ABUNDANCE'] = lvl_df[abu_col]/total_abundance

        lvl_df['ABUNDANCE_DEPTH'] = lvl_df['DEPTH']
        if lvl_df['DEPTH'].sum() > 0:
            lvl_df['REL_ABUNDANCE_DEPTH'] = lvl_df['ABUNDANCE_DEPTH']/lvl_df['ABUNDANCE_DEPTH'].sum()
        else:
            lvl_df['REL_ABUNDANCE_DEPTH'] = 0
        lvl_df['ABUNDANCE_GC'] = lvl_df['GENOMIC_CONTENT_EST']
        lvl_df['REL_ABUNDANCE_GC'] = lvl_df['GENOMIC_CONTENT_EST']/lvl_df['GENOMIC_CONTENT_EST'].sum()

        # if 'NOTE' is not empty, add '; ' to the end of the string
        lvl_df['NOTE'] = lvl_df['NOTE'].apply(lambda x: f'{x}; ' if x else x)

        # A note is added if a taxa has a rank with higher resolution than the target rank (signature-level), suggesting potential bias.
        # However, this note is not added if the taxa has a corresponding signature to support the observation.

        # add not shown reason
        idx = lvl_df['SIG_LEVEL'] < major_ranks[rank]
        lvl_df.loc[idx, 'NOTE'] += f"Not shown ({rank}-result could be biased); "

        # add SCORE reason
        if rank == 'strain':
            filtered = (lvl_df['SNI_SCORE'] < sni_score_strain)
            lvl_df.loc[filtered, 'NOTE'] += "Filtered out (strain SNI_SCORE > " + lvl_df.loc[filtered, 'SNI_SCORE'].astype(str) + "); "
        elif rank == 'species':
            filtered = (lvl_df['SNI_SCORE'] < sni_score_species)
            lvl_df.loc[filtered, 'NOTE'] += "Filtered out (species SNI_SCORE > " + lvl_df.loc[filtered, 'SNI_SCORE'].astype(str) + "); "
        else:
            filtered = (lvl_df['SNI_SCORE'] < sni_score_cutoff)
            lvl_df.loc[filtered, 'NOTE'] += "Filtered out (SNI_SCORE cutoff > " + lvl_df.loc[filtered, 'SNI_SCORE'].astype(str) + "); "


        # add filtered reason
        filtered = (lvl_df['BEST_SIG_COV'] < mc)
        lvl_df.loc[filtered, 'NOTE'] += "Filtered out (minCov > " + lvl_df.loc[filtered, 'BEST_SIG_COV'].astype(str) + "); "
        filtered = (lvl_df['READ_COUNT'] < mr)
        lvl_df.loc[filtered, 'NOTE'] += "Filtered out (minReads > " + lvl_df.loc[filtered, 'READ_COUNT'].astype(str) + "); "
        filtered = (lvl_df['COVERED_SIG_LEN'] < ml)
        lvl_df.loc[filtered, 'NOTE'] += "Filtered out (minLen > " + lvl_df.loc[filtered, 'COVERED_SIG_LEN'].astype(str) + "); "

        if mz > 0:
            filtered = (lvl_df['ZSCORE'] > mz)
            lvl_df.loc[filtered, 'NOTE'] += "Filtered out (maxZscore < " + lvl_df.loc[filtered, 'ZSCORE'].astype(str) + "); "

        # concart ranks-dataframe to the report-dataframe
        rep_df = pd.concat([lvl_df.sort_values('ABUNDANCE', ascending=False), rep_df], ignore_index=True)

    # add additional columns
    rep_df = rep_df.assign(
        SIG_COV                = rep_df["COVERED_SIG_LEN"]/rep_df["TOTAL_SIG_LEN"],
        READ_IDT               = (rep_df["TOTAL_BP_MAPPED"]-rep_df["TOTAL_BP_MISMATCH"])/rep_df["TOTAL_READ_LEN"],
        COVERED_MAPPED_SIG_COV = rep_df["COVERED_SIG_LEN"]/rep_df["MAPPED_SIG_LEN"],
        COVERED_SIG_DEPTH      = rep_df["TOTAL_BP_MAPPED"]/rep_df["COVERED_SIG_LEN"],
    )

    rep_df.drop(columns=['TAXID'], inplace=True)
    rep_df.rename(columns={"LVL_NAME": "NAME", "LVL_TAXID": "TAXID"}, inplace=True)

    logging.debug(f'rep_df: {rep_df}')

    return rep_df


def infer_sni_score(df, error_rate):
    """
    Estimate the Average Nucleotide Identity (SNI-score) together with 95% confidence intervals:
    - widens the interval when only a fraction of the signature space is actually covered ( SIG_COV )
    - project mismatches onto those unique positions
    - automatically becomes narrower as more signature bases are covered
    """

    df = df.copy()

    # from scipy.stats import norm
    # z = norm.ppf(0.5 + conf/2)  # ≈1.96
    z = 1.959963984540054

    # use only unique covered signature bases
    n = df["COVERED_SIG_LEN"]
    cov = df["SIG_COV"].clip(lower=1e-12)  # avoid n_eff = 0

    # remove the expected 0.5 % sequencing-error penalty
    m_rate = (df["TOTAL_BP_MISMATCH"]/df["TOTAL_BP_MAPPED"])
    m_rate_adj = m_rate - error_rate
    m_rate_adj = m_rate_adj.clip(lower=1e-12)  # avoid negative values

    # observed SCORE
    p_hat = 1 - m_rate_adj
    p = 1 - m_rate

    # cov is the coverage of the signature space, used to widen the confidence interval when only a fraction of the signature is covered
    n_eff  = n * cov

    z2     = z ** 2
    denom  = 1 + z2 / n_eff
    center = (p_hat + z2 / (2 * n_eff)) / denom
    hw     = (z * np.sqrt(
                (p_hat * (1 - p_hat)) / n_eff + z2 / (4 * n_eff ** 2)
             ) / denom)

    # observed-identity CI
    id_low, id_high = center - hw, center + hw

    # convert to true SCORE by adding the sequencing-error penalty
    score_naive = np.minimum(1, p)
    score_low  = np.clip(id_low, 0, 1)
    score_high = np.clip(id_high, 0, 1)

    score_ci95 = "[" + score_low.round(6).astype(str) + "-" + score_high.round(6).astype(str) + "]"

    df = df.assign(
        SNI_NAIVE    = score_naive.round(6),
        SNI_SCORE    = center.round(6),
        SNI_CI95_LH  = score_ci95
    )

    return df


def pile_lvl_zscore(tol_bp, tol_sig_len, linear_len):
    """
    Calculate Z-score for the depth of coverage of mapped regions.

    This determines how unusual the coverage depth is compared to expected depth
    based on a statistical model. Higher Z-scores may indicate biased mapping.

    Parameters:
        tol_bp (int): Total number of mapped bases
        tol_sig_len (int): Total length of the signature
        linear_len (int): Linear length (de-duplicated) covered by mappings

    Returns:
        float: Z-score for the depth distribution (or 0 if calculation fails)
    """
    try:
        avg_doc = tol_bp/tol_sig_len
        lin_doc = tol_bp/linear_len
        v = (linear_len*(lin_doc-avg_doc)**2 + (tol_sig_len-linear_len)*(avg_doc)**2)/tol_sig_len
        sd = math.sqrt(v)
        if sd == 0.0:
            return 0
        else:
            return (lin_doc-avg_doc)/sd
    except:
        return 0

def generate_taxonomy_file(rep_df, o, fullreport_o, fmt="tsv"):
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
    cols = ['LEVEL', 'NAME', 'TAXID', 'READ_COUNT', 'TOTAL_BP_MAPPED',
            'SNI_SCORE', 'COVERED_SIG_LEN', 'BEST_SIG_COV', 'DEPTH', 'REL_ABUNDANCE_GC', 'REL_ABUNDANCE', # summary
            'PARENT_NAME', 'PARENT_TAXID', # parants
            'AOI_READ_COUNT', 'TOTAL_READ_LEN', 'READ_IDT', 'TOTAL_BP_MISMATCH', 'TOTAL_BP_INDEL', 'SNI_NAIVE', 'SNI_CI95_LH', # read stats
            'SIG_COV', 'MAPPED_SIG_LEN', 'TOTAL_SIG_LEN', 'COVERED_SIG_DEPTH', 'COVERED_MAPPED_SIG_COV', 'ZSCORE', # signature stats
            'GENOMIC_CONTENT_EST', 'ABUNDANCE', 'REL_ABUNDANCE_DEPTH', # abundance
            'SIG_LEVEL', 'GENOME_COUNT', 'GENOME_SIZE', 'NOTE' # ref genome
            ]


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

def generate_biom_file(res_df, o, tg_rank, sampleid):
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

def generate_lineage_file(target_df, o):
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

def generate_mpa_file(target_df, o):
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

def readMapping(reads, db, threads, mm_options, presetx, samfile, logfile):
    """
    Map reads to the reference database using minimap2.

    Builds and executes a command to run minimap2 for read mapping, with parameters
    adjusted based on input settings. Filters the SAM output to keep only relevant
    alignments.

    Parameters:
        reads (list): List of input read file objects
        db (str): Path to the minimap2 database (without .mmi extension)
        threads (int): Number of threads to use
        mm_options (str): Minimap2 options for read mapping
        presetx (str): Minimap2 preset mode ('sr', 'map-pb', or 'map-ont')
        samfile (str): Output SAM file path
        logfile (str): Log file path
        nanopore (bool): Whether to use Nanopore-specific settings

    Returns:
        tuple: (
            exitcode (int): Exit code from the mapping process,
            cmd (str): Command that was executed,
            errs (str): Error output from the command
        )
    """
    input_file = " ".join([x.name for x in reads])

    # Minimap2 options for short reads: the options here is essentailly the -x 'sr' equivalent with some modifications on scoring
    sr_opts = f"-x sr {mm_options} -a -N20 --eqx --secondary=no --sam-hit-only"

    if presetx != 'sr':
        sr_opts = f"-x {presetx} -N20 --secondary=no --sam-hit-only -a"

    bash_cmd   = f"set -o pipefail; set -x;"
    mm2_cmd    = f"minimap2 {sr_opts} -t{threads} {db}.mmi {input_file}"
    filter_cmd = f"sed '/^@/d'"  # filter out header lines
    cmd        = f"{bash_cmd} {mm2_cmd} 2>> {logfile} | {filter_cmd} > {samfile}"

    logging.info(f"Readmapping command: {mm2_cmd}")

    proc = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    outs, errs = proc.communicate()
    exitcode = proc.poll()

    return exitcode, mm2_cmd, errs

def preprocess_nanopore_reads(reads, outdir, prefix, silent):
    """
    Split nanopore reads into shorter chunks before mapping.

    Parameters:
        reads (list): List of input read file objects (expects exactly one)
        outdir (str): Output directory for temporary files
        prefix (str): Prefix for the generated chunk file
        silent (bool): Silence flag for logging

    Returns:
        list: Updated list of read objects pointing to the chunked read file
    """
    os.makedirs(outdir, exist_ok=True)

    if len(reads) != 1:
        print_message("ERROR: Nanopore read processing expects a single input file.", silent, begin_t, logfile, errorout=1)

    input_path = reads[0].name
    output_path = os.path.join(outdir, f"{prefix}.split_reads.fasta.gz")

    print_message("Splitting nanopore reads into chunks...", silent, begin_t, logfile)
    try:
        chunk_count = ont_utils.split_to_fasta(input_path, output_path, split_length=150, step_length=150, drop_tail=True)
    except Exception as e:
        print_message(f"ERROR: Failed to split nanopore reads: {e}", silent, begin_t, logfile, errorout=1)
    else:
        if chunk_count == 0:
            print_message("ERROR: No reads were produced after splitting nanopore reads.", silent, begin_t, logfile, errorout=1)
        print_message(f" - {chunk_count} chunks written to {output_path}", silent, begin_t, logfile)

    return [SimpleNamespace(name=output_path)]

def split_reads_samfile_postprocessing(samfile, samfile_temp):
    """
    Clean up SAM file by removing inconsistent split-read alignments.

    Parameters:
        samfile (str): Path to the input SAM file
        samfile_temp (str): Path to the output cleaned SAM file

    Returns:
        bool: True if successful
    """
    total_chunks = 0

    logging.info(f'Reading split-read alignments from the sam file...')

    df = pd.read_csv(samfile,
                sep='\t',
                header=None,
                usecols=[0, 2],
                names=['QNAME', 'REF'],
                dtype={'QNAME': 'str', 'REF': 'str'}
    )

    df['QNAME_MAIN'] = df['QNAME'].str.split('|').str[0]
    df['REF_TAXID'] = df['REF'].str.split('|').str[-2]

    logging.info(f'Filtering out inconsistent chunks of reads...')
    # get the index with the most frequent taxid for each read
    idxmax = df.assign(_taxid_cnt=df.groupby(["QNAME_MAIN", "REF_TAXID"])["REF_TAXID"].transform("size")) \
                .loc[lambda x: x["_taxid_cnt"] == x.groupby("QNAME_MAIN")["_taxid_cnt"].transform("max")] \
                .index
    # Create a set of indices for faster lookup
    idxmax_set = set(idxmax.values)

    # in the idxmax, get the index of the first occurrence of the chunk (QNAME_MAIN) for each read
    idx1st = df.loc[idxmax].drop_duplicates(subset="QNAME_MAIN", keep="first").index
    idx1st_set = set(idx1st.values)

    total_chunks = len(df)
    del idxmax
    del df

    logging.info(f'Writing qualified hits...')
    with open(samfile_temp, 'w') as fout, open(samfile, 'r') as fin:
        for idx, line in enumerate(fin):
            if not idx%100000:
                logging.debug(f'Processed {idx} lines...')

            if idx in idxmax_set:
                if idx in idx1st_set:
                    fout.write(f"{line.rstrip()}\tZC:i:1\n")
                else:
                    fout.write(line)

    logging.info(f'Done writing {len(idxmax_set)} hits.')

    return True, total_chunks, len(idxmax_set)

def remove_multiple_hits(samfile, samfile_temp):
    """
    Removing multiple hits from the SAM file by keeping only the best alignment for each read.

    Parameters:
        samfile (str): Path to the SAM file
        samfile_temp (str): Path to the temporary SAM file with only the best alignments

    Returns:
        bool: False if no multiple hits were found, True if multiple hits were removed
    """
    logging.info(f'Loading the sam file...')

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

    logging.info(f'After removing non-parmary hits: {len(df)}')
    logging.info(f'Identifying top score hits...')
    # if FLAG bitwise AND with 128 (second in pair), append '/2' to the QNAME
    idx = (df['FLAG'] & 128).astype(bool)
    df.loc[idx, 'QNAME'] = df.loc[idx, 'QNAME'] + '/2'

    # get the index with the best alignment score for each read
    idxmax = df.groupby('QNAME')['AS'].idxmax()
    logging.info(f'Total top score hits: {len(idxmax)}')

    if len(idxmax) == aln_count:
        logging.info(f'No multiple hits found. Keeping the original SAM file.')
        return False
    else:
        # Create a set of indices for faster lookup
        idxmax_set = set(idxmax.values)
        del idxmax

        logging.info(f'Writing top score hits...')
        with open(samfile_temp, 'w') as fout, open(samfile, 'r') as fin:
            for idx, line in enumerate(fin):
                if not idx%100000:
                    logging.debug(f'Processed {idx} lines...')

                if idx in idxmax_set:
                    fout.write(line)
        logging.info(f'Done writing hits.')

        return True

def load_acc_list(filepath):
    """
    Load a list of accession numbers.

    Reads a file containing accession numbers (one per line) and returns them
    as a set for efficient lookup during processing. Empty lines are ignored.

    Parameters:
        filepath (str): Path to the file containing accession numbers

    Returns:
        set: Set of accession numbers. Returns empty set if input file is empty.

    Example:
        acc_list = load_acc_list('plasmid.txt')
    """
    with open(filepath) as f:
        acc_list = f.read().splitlines()

    if len(acc_list) == 0:
        logging.warning(f"Accession list is empty.")
        return set()
    else:
        return set(acc_list)

def loadDatabaseStats(db_stats_file):
    """
    Load database signature statistics from a stats file.

    Reads a tab-delimited stats file containing information about
    taxonomic signatures and their lengths.

    Parameters:
        db_stats_file (str): Path to the database stats file

    Returns:
        pd.Dataframe: df indexed with taxid, contains signature lengths and genome sizes

    Note:
        The input stats file is an 9-column tab-delimited file with:
        1. Rank
        2. Name
        3. Taxid
        4. Superkingdom
        5. NumOfSeq
        6. Max
        7. Min
        8. TotalLength
        9. GenomeSize
       10. Note (optional)
    """

    # Determine the number of columns in the stats file
    header = pd.read_csv(db_stats_file, nrows=0, sep='\t', header=None)
    valid_col_count = len(header.columns)

    usecols = [0, 2, 7, 8]
    names = ['DB_level', 'Taxid', 'TotalLength', 'Note']

    if valid_col_count == 10:
        usecols=[0, 2, 7, 8, 9]
        names=['DB_level', 'Taxid', 'TotalLength', 'GenomeSize', 'Note']

    # Set header to None to support files without headers
    df_stats = pd.read_csv(db_stats_file,
                           low_memory=False,
                           sep='\t',
                           header=None,
                           usecols=usecols,
                           names=names,
                           dtype={'DB_level': str, 'Taxid': str},
                           index_col='Taxid')

    # If 'Note' column is not present, create it with empty strings
    if not 'Note' in df_stats:
        df_stats['Note'] = ''
    if not 'GenomeSize' in df_stats:
        df_stats['GenomeSize'] = 0

    # Remove the row with index 'Taxid' if it exists
    # This is to handle the case when the stats file having the header
    if 'Taxid' in df_stats.index:
        df_stats = df_stats.drop('Taxid')

    # Make sure the format is consistent
    try:
        df_stats['TotalLength'] = pd.to_numeric(df_stats['TotalLength'], errors='raise')
    except ValueError:
        logging.error(f"Error processing stats file. Please check the format of the file.")
        sys.exit(1)

    # This is to handle the case when the stats file does not have GenomeSize column
    # In that case, the 'Note' column will be loaded as 'GenomeSize' and filled with 0
    df_stats['GenomeSize'] = pd.to_numeric(df_stats['GenomeSize'], errors='coerce')
    df_stats['GenomeSize'] = df_stats['GenomeSize'].fillna(0).astype(int)

    return df_stats

def print_message(msg, silent, start, logfile, errorout=0):
    """
    Print and log a timestamped message.

    Writes a message to the log file and optionally to stderr. Can also
    terminate the program with an error message.

    Parameters:
        msg (str): Message to print
        silent (bool): If True, suppress output to stderr
        start (float): Start time for timestamp calculation
        logfile (str): Path to the log file
        errorout (int): If non-zero, exit with error after printing

    Returns:
        None

    Raises:
        SystemExit: If errorout is non-zero
    """
    message = "[%s] %s\n" % (time_spend(start), msg)

    with open( logfile, "a" ) as f:
        f.write( message )
        f.close()

    if errorout:
        sys.exit( message )
    elif not silent:
        sys.stderr.write( message )

def parse_taxids(taxid_arg, res_df, full_tsv_fn):
    """Parse taxids from command line arg or file"""

    taxa_list = []
    qualified_taxa = pd.DataFrame()
    taxa_df = pd.DataFrame()

    if res_df.shape[1] > 0:
        taxa_df = res_df
        print_message(f"Successfully loaded taxonomy profile with {len(taxa_df)} entries",
                    argvs.silent, begin_t, logfile)
    else:
        try:
            print_message(f"Reading taxonomy file {full_tsv_fn}...",
                        argvs.silent, begin_t, logfile)
            taxa_df = pd.read_csv(full_tsv_fn,
                                sep='\t',
                                engine='python',
                                quoting=3,
                                on_bad_lines='skip',
                                dtype={'NOTE': str})

            print_message(f"Successfully loaded taxonomy profile with {len(taxa_df)} entries",
                        argvs.silent, begin_t, logfile)
        except Exception as e:
            print_message(f"Error reading taxonomy file: {e}", argvs.silent, begin_t, logfile, errorout=1)

    # Filter in entries specified by taxid_arg
    filtered_idx = None

    if 'NOTE' in taxa_df.columns:
        filtered_idx = ~taxa_df['NOTE'].str.contains('Filtered out', na=False)

    if taxid_arg and taxid_arg != 'all':
        if taxid_arg.startswith('@'):
            # Read taxids from file
            filename = taxid_arg[1:]  # Remove @ prefix
            try:
                with open(filename) as f:
                    taxa_list = [x.strip() for x in f.readlines() if x.strip() and not x.startswith('#')]
            except IOError as e:
                sys.stderr.write(f"Error reading taxid file {filename}: {e}\n")
                sys.exit(1)
        else:
            # Parse comma-separated list
            taxa_list = [x.strip() for x in taxid_arg.split(',')]

        if taxa_list:
            filtered_idx &= (taxa_df['TAXID'].isin(taxa_list) | taxa_df['NAME'].isin(taxa_list))

    # Filter out entries with "Filtered out" notes
    if filtered_idx is not None:
        qualified_taxa = taxa_df[filtered_idx]
        print_message(f"Found {len(qualified_taxa)} qualified taxa after filtering",
                    argvs.silent, begin_t, logfile)
    else:
        qualified_taxa = taxa_df

    # Ensure these columns exist
    if not all(col in qualified_taxa.columns for col in ['LEVEL', 'NAME', 'TAXID']):
        print_message(f"Required columns missing in taxonomy file. Available columns: {qualified_taxa.columns.tolist()}",
                    argvs.silent, begin_t, logfile, errorout=1)

    # Pre-compute a mapping from reference taxids to qualified taxids
    # This avoids expensive lineage lookups during processing
    print_message("Building taxonomy lookup index...",
                argvs.silent, begin_t, logfile)

    taxa_dict = {}

    # Gather all qualified taxids
    qualified_taxids = []
    for _, row in qualified_taxa[['LEVEL', 'NAME', 'TAXID']].iterrows():
        if pd.notna(row['TAXID']):
            taxid = str(row['TAXID']).strip()
            qualified_taxids.append(taxid)
            taxa_dict[taxid] = {
                'level': str(row['LEVEL']).replace(' ', '_') if pd.notna(row['LEVEL']) else 'unknown',
                'name': str(row['NAME']).replace(' ', '_') if pd.notna(row['NAME']) else 'unknown'
            }

    print_message(f"Starting extraction for {len(qualified_taxids)} taxa...",
                argvs.silent, begin_t, logfile)

    logging.debug(f"Qualified taxa: {qualified_taxa}")
    logging.debug(f"Qualified taxids: {qualified_taxids}")

    return taxa_dict, qualified_taxids

def main(args):
    """
    Main execution function for GOTTCHA2.
    """
    global argvs
    global logfile
    global begin_t
    global df_stats
    global acc_list

    argvs = parse_args( __version__, args )
    begin_t  = time.time()
    sam_fp   = argvs.sam[0] if argvs.sam else ""
    samfile  = f"{argvs.outdir}/{argvs.prefix}.gottcha_{argvs.dbLevel}.sam" if not argvs.sam else sam_fp.name
    logfile  = f"{argvs.outdir}/{argvs.prefix}.gottcha_{argvs.dbLevel}.log"
    set_start_method("fork") # for default multiprocessing method
    acc_list = set()
    split_read_flag = False
    res_df = pd.DataFrame() # aggregated restuls

    logging_level = logging.WARNING

    if argvs.debug:
        logging_level = logging.DEBUG
    elif argvs.silent:
        logging_level = logging.FATAL
    elif argvs.verbose:
        logging_level = logging.INFO

    logging.basicConfig(
        level=logging_level,
        format='%(asctime)s [%(levelname)s] %(module)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M',
    )

    #dependency check
    if sys.version_info < (3,8):
        sys.exit("[ERROR] Python 3.8 or above is required.")

    dependency_check("minimap2")
    # dependency_check("gawk")

    #prepare output object
    argvs.relAbu = argvs.relAbu.upper()
    outfile_full = "%s/%s.full.tsv" % (argvs.outdir, argvs.prefix)
    outfile_lineage = "%s/%s.lineage.tsv" % (argvs.outdir, argvs.prefix)
    outfile_mpa = "%s/%s.mpa.tsv" % (argvs.outdir, argvs.prefix)

    # remove previous log file if exists
    if os.path.isfile(logfile):
        os.remove(logfile)

    out_fp = sys.stdout
    outfile = "STDOUT"

    if not argvs.stdout and not argvs.extractOnly:
        #create output directory if not exists
        if not os.path.exists(argvs.outdir):
            os.makedirs(argvs.outdir)

        outfile = "%s/%s.tsv" % (argvs.outdir, argvs.prefix)
        if argvs.format == "csv":

            outfile = "%s/%s.csv" % (argvs.outdir, argvs.prefix)
        elif argvs.format == "biom":
            outfile = "%s/%s.biom" % (argvs.outdir, argvs.prefix)

        out_fp = open(outfile, 'w')

    # display the command line
    logging.info( ' '.join(sys.argv) )

    print_message(f"GOTTCHA (v{__version__})", argvs.silent, begin_t, logfile)
    print_message(f"Arguments and dependencies checked:", argvs.silent, begin_t, logfile)
    if argvs.input:
        print_message(f"    Input reads        : {[x.name for x in argvs.input]}",     argvs.silent, begin_t, logfile)
    print_message(f"    Input SAM file     : {samfile}",           argvs.silent, begin_t, logfile)
    print_message(f"    Database           : {argvs.database}",    argvs.silent, begin_t, logfile)
    print_message(f"    Database level     : {argvs.dbLevel}",     argvs.silent, begin_t, logfile)
    print_message(f"    Abundance          : {argvs.relAbu}",      argvs.silent, begin_t, logfile)
    print_message(f"    Output path        : {argvs.outdir}",      argvs.silent, begin_t, logfile)
    print_message(f"    Prefix             : {argvs.prefix}",      argvs.silent, begin_t, logfile)
    print_message(f"    Threads            : {argvs.threads}",     argvs.silent, begin_t, logfile)
    print_message(f"    SNI-score (g,s,n)  : {argvs.sniScore}",    argvs.silent, begin_t, logfile)
    if argvs.nanopore:
        print_message(f"    Nanopore mode      : Enabled", argvs.silent, begin_t, logfile)
    if argvs.errorRate >= 0.0:
        print_message(f"    Read error rate    : {argvs.errorRate}", argvs.silent, begin_t, logfile)
    if argvs.accList:
        print_message(f"    Acc-of-int list    : {argvs.accList}", argvs.silent, begin_t, logfile)
    if argvs.accList:
        print_message(f"    Acc-of-int action  : {argvs.accListAction}", argvs.silent, begin_t, logfile)
    if argvs.extract:
        print_message(f"    Extract seqs       : {argvs.extract}", argvs.silent, begin_t, logfile)
    if argvs.minCov > 0:
        print_message(f"    Minimal SIG cov    : {argvs.minCov}", argvs.silent, begin_t, logfile)
    if argvs.minLen > 0:
        print_message(f"    Minimal SIG length : {argvs.minLen}", argvs.silent, begin_t, logfile)
    if argvs.maxZscore > 0:
        print_message(f"    Maximal zScore     : {argvs.maxZscore}", argvs.silent, begin_t, logfile)
    if argvs.minReads > 0:
        print_message(f"    Minimal reads      : {argvs.minReads}", argvs.silent, begin_t, logfile)
    if argvs.matchIdentity > 0:
        print_message(f"    Min match identity : {argvs.matchIdentity}", argvs.silent, begin_t, logfile)
    if argvs.matchFraction > 0:
        print_message(f"    Min match fraction : {argvs.matchFraction}", argvs.silent, begin_t, logfile)
    if argvs.matchLength > 0:
        print_message(f"    Min match length   : {argvs.matchLength}", argvs.silent, begin_t, logfile)

    #load taxonomy
    print_message("Loading taxonomy information...", argvs.silent, begin_t, logfile)
    custom_taxa_tsv = None
    dbpath = None
    if os.path.isdir(argvs.taxInfo):
        dbpath = argvs.taxInfo
    elif os.path.isfile(argvs.taxInfo):
        custom_taxa_tsv = argvs.taxInfo
    elif os.path.isfile( argvs.database + ".tax.tsv" ):
        custom_taxa_tsv = argvs.database + ".tax.tsv"

    logging.info(f"Taxonomy file: {custom_taxa_tsv}")

    taxonomy.loadTaxonomy(dbpath=dbpath,
                    cus_taxonomy_file=custom_taxa_tsv,
                    auto_download=False)
    print_message(f" - {len(taxonomy.taxNames)} taxa loaded.", argvs.silent, begin_t, logfile)

    #load database stats
    print_message("Loading database stats...", argvs.silent, begin_t, logfile)
    if os.path.isfile( argvs.database + ".stats" ):
        df_stats = loadDatabaseStats(argvs.database+".stats")
    else:
        print_message(f"ERROR: {argvs.database+'.stats'} not found.", argvs.silent, begin_t, logfile, errorout=1)

    print_message(f" - {df_stats.shape[0]} entries loaded.", argvs.silent, begin_t, logfile)
    print_message(f" - signatures at {df_stats['DB_level'].unique().tolist()} levels loaded.", argvs.silent, begin_t, logfile)

    if argvs.accList:
        print_message("Loading accession#s of interest list...", argvs.silent, begin_t, logfile)
        acc_list = load_acc_list(argvs.accList)
        print_message(f" - {len(acc_list)} accession#s of interest loaded.", argvs.silent, begin_t, logfile)

    #main process
    if argvs.input:
        # if nanopore option is on, preprocessing reads
        if argvs.nanopore:
            print_message("Checking nanopore read files...", argvs.silent, begin_t, logfile)
            argvs.input = preprocess_nanopore_reads(argvs.input, argvs.outdir, argvs.prefix, argvs.silent)
            split_read_flag = True

        print_message("Running read-mapping...", argvs.silent, begin_t, logfile)
        exitcode, cmd, msg = readMapping( argvs.input, argvs.database, argvs.threads, argvs.m2options, argvs.presetx, samfile, logfile)
        print_message(f"Logfile saved to {logfile}.", argvs.silent, begin_t, logfile)
        print_message(f"COMMAND: {cmd}", argvs.silent, begin_t, logfile)

        if exitcode != 0:
            # if size of the samfile is zero
            sys.exit( "[%s] ERROR: error occurred while running read mapping (exit: %s, message: %s).\n" % (time_spend(begin_t), exitcode, msg) )
        else:
            print_message(f"Done mapping reads to {argvs.dbLevel} signature database.", argvs.silent, begin_t, logfile)
            print_message(f"Mapped SAM file saved to {samfile}.", argvs.silent, begin_t, logfile)
            sam_fp = open( samfile, "r" )
        gc.collect()

    # remove multiple hits
    if argvs.removeMultipleHits == 'yes':
        # remove multiple hits from the SAM file
        print_message("Removing multiple hits from SAM file...", argvs.silent, begin_t, logfile)
        samfile_output = f"{argvs.outdir}/{argvs.prefix}.gottcha_{argvs.dbLevel}.sam"
        samfile_temp = f"{argvs.outdir}/{argvs.prefix}.gottcha_{argvs.dbLevel}.sam.temp"
        flag = remove_multiple_hits(samfile, samfile_temp)
        if flag:
            os.rename(samfile_temp, samfile_output)
            samfile = samfile_output
            # Note:
            # When input of the gottcha2 is a SAM file and new outdir/prefix is provided, the output will be saved to that location.
            # If not, the output will overwrite the original SAM file.
        gc.collect()

    # preprocess SAM file for nanopore reads
    if argvs.nanopore:
        # remove inconsistent read chunks from the SAM file
        print_message("Removing inconsistent read chunks from SAM file...", argvs.silent, begin_t, logfile)
        samfile_output = f"{argvs.outdir}/{argvs.prefix}.gottcha_{argvs.dbLevel}.sam"
        samfile_temp = f"{argvs.outdir}/{argvs.prefix}.gottcha_{argvs.dbLevel}.sam.temp"
        flag, tol_chunks_count, tol_chunks_qualified = split_reads_samfile_postprocessing(samfile, samfile_temp)
        if flag:
            os.rename(samfile_temp, samfile_output)
            samfile = samfile_output
        print_message(f" - {tol_chunks_count} mapped read chunks processed", argvs.silent, begin_t, logfile)
        print_message(f" - {tol_chunks_count-tol_chunks_qualified} inconsistent hits removed", argvs.silent, begin_t, logfile)
        gc.collect()

    # processing SAM file and generate results
    if not argvs.extractOnly:
        print_message("Processing SAM file...", argvs.silent, begin_t, logfile)
        (res, mapped_r_cnt, tol_alignment_count, tol_invalid_match_count) = process_sam_file( os.path.abspath(samfile), argvs.threads, argvs.matchFraction, argvs.matchIdentity, argvs.matchLength, split_read_flag)
        gc.collect()

        print_message(f" - {tol_alignment_count} alignments processed", argvs.silent, begin_t, logfile)
        print_message(f" - {tol_invalid_match_count} alignments did not meet matching criteria", argvs.silent, begin_t, logfile)
        print_message(f" - {mapped_r_cnt} qualified mapped reads", argvs.silent, begin_t, logfile)

        if mapped_r_cnt:
            # Set SNI-SCORE default to 0.8, species 0.95, strain 0.99
            if ',' not in argvs.sniScore:
                argvs.sniScore = ','.join([argvs.sniScore]*3)
            elif argvs.sniScore.count(',') == 1:
                argvs.sniScore = argvs.sniScore + ',0.99'

            (sni_score_cutoff, sni_score_species, sni_score_strain) = [float(x) for x in argvs.sniScore.split(',')]

            # agg signature fragments to strains, and the read count of the accession#s of interest
            str_df, aoi_read_count = group_refs_to_strains(res, acc_list, argvs.accListAction)

            # aggregate the results
            _args = (str_df,
                     argvs.relAbu,
                     argvs.dbLevel,
                     argvs.minCov,
                     argvs.minReads,
                     argvs.minLen,
                     argvs.maxZscore,
                     sni_score_species,
                     sni_score_strain,
                     sni_score_cutoff,
                     argvs.errorRate
                     )
            res_df = aggregate_taxonomy(*_args)

            if acc_list:
                print_message(f" - {aoi_read_count} reads mapped to accession-of-interest", argvs.silent, begin_t, logfile)
                read_count_after_aoi = mapped_r_cnt
                if argvs.accListAction == 'filter_out':
                    read_count_after_aoi = mapped_r_cnt - aoi_read_count
                elif argvs.accListAction == 'filter_in':
                    read_count_after_aoi = aoi_read_count
                print_message(f" - {read_count_after_aoi} reads after applying accession-of-interest action ({argvs.accListAction})", argvs.silent, begin_t, logfile)

            print_message("Done taxonomy aggregation.", argvs.silent, begin_t, logfile)

            if not len(res_df):
                print_message("No qualified taxonomy profiled.", argvs.silent, begin_t, logfile)
            else:
                # generate output results
                if argvs.format == "biom":
                    generate_biom_file(res_df, out_fp, argvs.dbLevel, argvs.prefix)
                else:
                    generate_taxonomy_file(res_df, out_fp, outfile_full, argvs.format)
                # generate lineage file
                target_idx = (res_df['LEVEL']==argvs.dbLevel) & \
                                (res_df['NOTE'].str.contains('Filtered out', na=False) == False) & \
                                (res_df['NOTE'].str.contains('Not shown', na=False) == False)
                target_df = res_df.loc[target_idx, ['ABUNDANCE','TAXID']]
                tax_num = len(target_df)

                print_message(f"{tax_num} qualified {argvs.dbLevel} profiled.", argvs.silent, begin_t, logfile)

                if tax_num:
                    generate_lineage_file(target_df, outfile_lineage)

                    if argvs.mpa:
                        target_df = res_df.loc[target_idx, ['TAXID', 'REL_ABUNDANCE', 'REL_ABUNDANCE_GC','READ_COUNT', 'SIG_COV']]
                        generate_mpa_file(target_df, outfile_mpa)
                        print_message(f"MPA format file saved to {outfile_mpa}.", argvs.silent, begin_t, logfile)

                print_message(f"Results saved to {outfile}.", argvs.silent, begin_t, logfile)
        else:
            print_message("GOTTCHA2 stopped.", argvs.silent, begin_t, logfile)
            sys.exit(0)

    # extracting reads
    if argvs.extract:
        (taxa_arg, max_per_taxon, out_format) = (argvs.extract.split(':', maxsplit=2) + ['all', 'fasta'])[:3]

        print_message(f"Extracting {len(taxa_arg)} taxa, {max_per_taxon} sequences per taxa to {out_format}...", argvs.silent, begin_t, logfile)

        full_report_file = ""

        if max_per_taxon.isdigit() or max_per_taxon == 'all':
            max_per_taxon = int(max_per_taxon) if max_per_taxon != 'all' else 0

        if argvs.extractOnly:
            full_report_file = '.'.join(os.path.abspath(samfile).split('.')[:-2]) + ".full.tsv"

        taxa_dict, qualified_taxids = parse_taxids(taxa_arg, res_df, full_report_file)

        if not len(qualified_taxids):
            print_message("No qualified taxonomy profiled.", argvs.silent, begin_t, logfile)

        if not argvs.stdout:
            outfile = f"{argvs.outdir}/{argvs.prefix}.extract.{out_format.lower()}"
            out_fp = open(outfile, 'w')

            # Extract FASTA sequences based on the taxonomy entries
            print_message(f"Extracting up to {max_per_taxon} {out_format} sequences per taxon...", argvs.silent, begin_t, logfile)

            _args = (os.path.abspath(samfile),
                       taxa_dict,
                       qualified_taxids,
                       out_fp,
                       argvs.threads,
                       argvs.matchFraction,
                       argvs.matchIdentity,
                       argvs.matchLength,
                       max_per_taxon,
                       acc_list,
                       argvs.accListAction,
                       out_format)
            taxon_count, seq_count = extract_sequences_by_taxonomy(*_args)
            print_message(f"Done extracting {seq_count} sequences from {taxon_count} taxa to '{outfile}'.",
                            argvs.silent, begin_t, logfile)

if __name__ == '__main__':
    main(sys.argv[1:])
