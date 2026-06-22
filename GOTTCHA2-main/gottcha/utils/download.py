#!/usr/bin/env python3

#This script allows the user to pull the existing gottcha database
import requests
import sys
import os
import tarfile
import hashlib
from tqdm import *
import argparse
from gottcha import GOTTCHA_DB_FULL_LATEST, GOTTCHA_DB_FAST_LATEST
import logging

def calculate_sha256(file_path, chunk_size=8192):
    sha256 = hashlib.sha256()
    with open(file_path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b''):
            sha256.update(chunk)
    return sha256.hexdigest()


def parse_params(args):
    parser = argparse.ArgumentParser(prog='download.py', description="""This script will pull the latest version of the Gottcha2 database.""")

    parser.add_argument('-u', '--url', required=False,
                    help='specify a URL to pull from (will override the default)')
    parser.add_argument('-d', '--database', default='regular',
                    help='specify the type of database to download (regular or fast)')
    parser.add_argument('-r', '--rank', default='species',
                    help='taxonomic rank of the database (superkingdom, phylum, class, order, famiily, genus, species)')
    return parser.parse_args(args)


def download_db(argvs):
    if os.path.isdir('database'):
        sys.exit('Please make sure a database directory does not exist.')

    os.mkdir('database')

    if not argvs.url:
        if argvs.database == 'regular':
            download_url = GOTTCHA_DB_FULL_LATEST
        elif argvs.database == 'fast':
            download_url = GOTTCHA_DB_FAST_LATEST
        else:
            sys.exit('Invalid database type specified. Please choose "regular" or "fast".')

    if argvs.url:
        download_url = argvs.url
        
    archive_name = os.path.basename(download_url)

    logging.info(f"Downloading GOTTCHA2 database from {download_url}...")

    with requests.get(download_url, stream=True) as r:
        r.raise_for_status()
        with open(archive_name, 'wb') as f:
            total_size = int(r.headers.get('Content-Length', 0))
            pbar = tqdm(total=total_size, unit='B', unit_scale=True, unit_divisor=1024)
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
            pbar.close()

    logging.info(f"Downloading GOTTCHA2 database SHA256 checksum from {download_url}.sha256...")

    with requests.get(f"{download_url}.sha256", stream=True) as r:
        r.raise_for_status()
        with open(f"{archive_name}.sha256", 'wb') as f:
            total_size = int(r.headers.get('Content-Length', 0))
            pbar = tqdm(total=total_size, unit='B', unit_scale=True, unit_divisor=1024)
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
            pbar.close()

    logging.info("Verifying SHA256 checksum...")

    with open(f"{archive_name}.sha256", 'r') as f:
        expected_sha256 = f.read().strip().split()[0]
    actual_sha256 = calculate_sha256(archive_name)
    if actual_sha256 != expected_sha256:
        os.remove(archive_name)
        os.remove(f"{archive_name}.sha256")
        sys.exit("SHA256 checksum does not match expected value. Download may be corrupted. Please try again.")

    logging.info("SHA256 checksum verified successfully.")

    logging.info(f"Extracting GOTTCHA2 database from {archive_name}...")
    with tarfile.open(archive_name) as tar:
        tar.extractall('database')
    logging.info("Database extraction completed.")
    os.remove(archive_name)
    os.remove(f"{archive_name}.sha256")
    logging.info("Temporary files removed. Database is ready to use.")

def main(args):
    argvs = parse_params(args)
    download_db(argvs)


if __name__ == '__main__':
    main(sys.argv[1:])
