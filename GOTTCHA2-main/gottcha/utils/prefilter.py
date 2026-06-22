#!/usr/bin/env python3
"""
quant.py - Wrapper for running sylph query command

Can be used as a standalone script or imported as a module.
"""

import subprocess
import sys
import argparse
from pathlib import Path
from typing import List, Optional, Union
import logging


def run_sylph_sketch(
    read_file: str,
    outdir: str,
    threads: int = 1,
    subsampling_rate: int = 100,
    additional_args: Optional[List[str]] = None
) -> subprocess.CompletedProcess:
    """
    Run sylph sketch command.
    
    Args:
        database: Path to the sylph database
        read: Path to read file
        output: Output file path
        threads: Number of threads to use (default: 1)
        subsampling_rate: Subsampling rate (default: 100)
        additional_args: Additional arguments to pass to sylph sketch
    
    Returns:
        CompletedProcess object containing result information
    
    Raises:
        FileNotFoundError: If database or read file don't exist
        subprocess.CalledProcessError: If sylph sketch fails
    """
    # Handle read as single string
    if not Path(read_file).exists():
        raise FileNotFoundError(f"Reads file not found: {read_file}")
    
    # Build command
    cmd = [
        "sylph", "sketch",
        "-t", str(threads),
        "-c", str(subsampling_rate),
        "-d", str(outdir),
        "-r", str(read_file),
    ]
    
    # Add additional arguments if provided
    if additional_args:
        cmd.extend(additional_args)
    
    logging.info(f"COMMAND: {' '.join(cmd)}")
    
    # Run command
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        return result
    except subprocess.CalledProcessError as e:
        logging.info(f"Error running sylph sketch:")
        logging.info(f"Command: {' '.join(cmd)}")
        logging.info(f"Return code: {e.returncode}")
        logging.info(f"stderr: {e.stderr}")
        raise

def run_sylph_query(
    database: str,
    reads: Union[str, List[str]],
    output: str,
    threads: int = 1,
    subsampling_rate: int = 100,
    minimum_ani: float = 80.0,
    minimum_kmer: int = 50,
    read_seq_id: float = 99.5,
    additional_args: Optional[List[str]] = None
) -> subprocess.CompletedProcess:
    """
    Run sylph query command.
    
    Args:
        database: Path to the sylph database
        reads: Path to reads file(s) (can be a string or list of strings)
        output: Output file path
        threads: Number of threads to use (default: 1)
        subsampling_rate: Subsampling rate (default: 100)
        minimum_ani: Minimum ANI threshold (default: 80.0)
        additional_args: Additional arguments to pass to sylph query
    
    Returns:
        CompletedProcess object containing result information
    
    Raises:
        FileNotFoundError: If database or reads files don't exist
        subprocess.CalledProcessError: If sylph query fails
    """
    # Validate inputs
    if not Path(database).exists():
        raise FileNotFoundError(f"Database not found: {database}")
    
    # Handle reads as list or single string
    if isinstance(reads, str):
        reads = [reads]
    
    for read_file in reads:
        if not Path(read_file).exists():
            raise FileNotFoundError(f"Reads file not found: {read_file}")
    
    # Build command
    cmd = [
        "sylph", "query",
        "-t", str(threads),
        "-c", str(subsampling_rate),
        "--read-seq-id", str(read_seq_id),
        "--min-number-kmers", str(minimum_kmer),
        "--minimum-ani", str(minimum_ani),
        "-o", output,
    ]
    
    # Add additional arguments if provided
    if additional_args:
        cmd.extend(additional_args)
    
    # Add database and reads
    cmd.append(database)
    cmd.extend(reads)
    
    logging.info(f"COMMAND: {' '.join(cmd)}")
    
    # Run command
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        return result
    except subprocess.CalledProcessError as e:
        logging.info(f"Error running sylph query:")
        logging.info(f"Command: {' '.join(cmd)}")
        logging.info(f"Return code: {e.returncode}")
        logging.info(f"stderr: {e.stderr}")
        raise


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Run sylph query command",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -d database.syldb -r reads.fq -o output.tsv
  %(prog)s -d database.syldb -r reads1.fq reads2.fq -o output.tsv -t 8
  %(prog)s -d database.syldb -r reads.fq -o output.tsv --minimum-ani 70 -c 50
        """
    )
    
    parser.add_argument(
        "-d", "--database",
        required=True,
        help="Path to sylph database"
    )
    
    parser.add_argument(
        "-r", "--reads",
        required=True,
        nargs="+",
        help="Path to reads file(s)"
    )
    
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file path"
    )
    
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        help="Number of threads (default: 1)"
    )
    
    parser.add_argument(
        "-c", "--coverage-threshold",
        type=int,
        default=100,
        help="Coverage threshold (default: 100)"
    )
    
    parser.add_argument(
        "--minimum-ani",
        type=float,
        default=60.0,
        help="Minimum ANI threshold (default: 60.0)"
    )
    
    parser.add_argument(
        "--extra-args",
        nargs=argparse.REMAINDER,
        help="Additional arguments to pass to sylph query"
    )
    
    args = parser.parse_args()
    
    try:
        result = run_sylph_query(
            database=args.database,
            reads=args.reads,
            output=args.output,
            threads=args.threads,
            coverage_threshold=args.coverage_threshold,
            minimum_ani=args.minimum_ani,
            additional_args=args.extra_args,
        )
        
        logging.info("Command completed successfully!")
        if result.stdout:
            logging.info("stdout:", result.stdout)
    
        return 0
        
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        logging.info(f"Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())