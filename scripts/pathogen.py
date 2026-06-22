#!/usr/bin/env python3

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

try:
    # Try relative import first (for package usage)
    from . import taxonomy as gt
except ImportError:
    # Fall back to direct import (for script usage)
    import taxonomy as gt

DEFAULT_LEVELS = {"strain", "species", "genus", "family"}
DEFAULT_PATHOGEN_TSV = Path(__file__).resolve().parent.parent / "data" / "pathogen.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate GOTTCHA2 results with pathogen information.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input GOTTCHA2 TSV file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV file with pathogen annotations.",
    )
    parser.add_argument(
        "-p",
        "--pathogen-tsv",
        default=str(DEFAULT_PATHOGEN_TSV),
        help="Pathogen list TSV file.",
    )
    parser.add_argument(
        "-d",
        "--dbpath",
        default=None,
        help="Taxonomy database directory.",
    )
    parser.add_argument(
        "-c",
        "--custom-taxonomy",
        dest="cus_taxonomy_file",
        default=None,
        help="Custom taxonomy TSV file.",
    )
    parser.add_argument(
        "--levels",
        default=",".join(sorted(DEFAULT_LEVELS)),
        help="Comma-separated list of taxonomic levels to check.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )
    return parser.parse_args()


def setup_logging(verbose: bool) -> None:
    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def _validate_columns(df: pd.DataFrame, required: set, label: str) -> None:
    missing = required - set(df.columns)
    if missing:
        missing_str = ", ".join(sorted(missing))
        raise ValueError(f"{label} is missing required columns: {missing_str}")


def annotate_pathogens(
    pathogen_tsv: str,
    input_tsv: str,
    levels: set,
) -> pd.DataFrame:
    df_pg = pd.read_csv(pathogen_tsv, sep="\t")
    _validate_columns(df_pg, {"TAXID", "NAME", "HOST"}, "Pathogen TSV")
    df_pg["TAXID"] = df_pg["TAXID"].astype(str)

    df_taxa = pd.read_csv(input_tsv, sep="\t")
    # if empty dftaxa, return directly
    if df_taxa.empty:
        return df_taxa
    
    _validate_columns(df_taxa, {"TAXID", "NAME", "LEVEL"}, "Input TSV")
    df_taxa["TAXID"] = df_taxa["TAXID"].astype(str)
    df_taxa["PATHOGENIC_INFO"] = ""
    df_taxa["HUMAN_PATHOGEN"] = ""

    name_set = set(df_pg["NAME"].astype(str))
    taxid_set = set(df_pg["TAXID"])

    for idx, row in df_taxa.iterrows():
        if row["LEVEL"] not in levels:
            continue

        name = str(row["NAME"])
        taxid = str(row["TAXID"])
        is_pathogen = False
        is_human = False
        pathogenic_text = ""

        if name in name_set:
            pidx = df_pg["NAME"] == name
            pathogenic_text = df_pg[pidx].to_html(index=False).replace('\n', '')
            is_pathogen = True
            if "human" in str(df_pg[pidx]["HOST"].values).lower():
                is_human = True
        elif taxid in taxid_set:
            pidx = df_pg["TAXID"] == taxid
            pathogenic_text = df_pg[pidx].to_html(index=False).replace('\n', '')
            is_pathogen = True
            if "human" in str(df_pg[pidx]["HOST"].values).lower():
                is_human = True

        if is_pathogen==False and row["LEVEL"]=="strain":
            parent_name = gt.taxid2nameOnRank(taxid, "species")
            parent_taxid = gt.taxid2taxidOnRank(taxid, "species")
            if parent_name in name_set:
                pidx = df_pg["NAME"] == parent_name
                pathogenic_text = df_pg[pidx].to_html(index=False).replace('\n', '')
                is_pathogen = True
                if "human" in str(df_pg[pidx]["HOST"].values).lower():
                    is_human = True
            elif parent_taxid in taxid_set:
                pidx = df_pg["TAXID"] == parent_taxid
                pathogenic_text = df_pg[pidx].to_html(index=False).replace('\n', '')
                is_pathogen = True
                if "human" in str(df_pg[pidx]["HOST"].values).lower():
                    is_human = True

        df_taxa.at[idx, "PATHOGENIC_INFO"] = pathogenic_text if is_pathogen else ""
        
        if is_pathogen and is_human:
            df_taxa.at[idx, "HUMAN_PATHOGEN"] = "Yes"
        elif is_pathogen:
            df_taxa.at[idx, "HUMAN_PATHOGEN"] = "No"

    return df_taxa

def pathogen_summary(df_taxa: pd.DataFrame) -> str:
    summary = ""
    idx = (df_taxa["LEVEL"]=='species') & (df_taxa["PATHOGENIC_INFO"] != "")
    total_pathogens = idx.sum()
    idx = (df_taxa["LEVEL"]=='species') & (df_taxa["HUMAN_PATHOGEN"] == "Yes")
    human_pathogens = idx.sum()

    if human_pathogens > 0:
        summary += f"{human_pathogens} pathogenic species known to infect humans were detected: "

        df_h_patho = df_taxa.loc[idx].copy()
        df_h_patho['CATE'] = df_h_patho['TAXID'].apply(lambda x: gt.taxid2nameOnRank(x, 'superkingdom'))

        # eidx = (df_h_patho['CATE'] == "Eukaryota")
        # if eidx.sum() > 0:
        #     df_h_patho.loc[eidx, 'CATE'] = df_h_patho[eidx]['TAXID'].apply(lambda x: gt.taxid2nameOnRank(x, 'kingdom'))
        
        summary += "; ".join(
            f"{cate}: {', '.join(names)}"
            for cate, names in (
                df_h_patho.groupby("CATE")["NAME"].apply(list).items()
            )
        )
    elif total_pathogens > 0:
        summary += f"Total {total_pathogens} pathogens. None of them are human pathogens."
    else:
        summary += "No pathogens identified."

    return summary


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    levels = {lvl.strip() for lvl in args.levels.split(",") if lvl.strip()}
    if not levels:
        logging.error("No valid levels specified in --levels.")
        return 1

    try:
        gt.loadTaxonomy(
            dbpath=args.dbpath,
            cus_taxonomy_file=args.cus_taxonomy_file,
        )
    except Exception as exc:
        logging.error("Failed to load taxonomy: %s", exc)
        return 1

    try:
        df_taxa = annotate_pathogens(
            pathogen_tsv=args.pathogen_tsv,
            input_tsv=args.input,
            levels=levels,
        )

        if df_taxa.empty:
            logging.info("No taxa found in input file.")
            print("No taxa found in input file.")
            return 0

    except Exception as exc:
        logging.error("Failed to annotate pathogens: %s", exc)
        return 1

    summary = pathogen_summary(df_taxa)

    print(summary)

    df_taxa.to_csv(args.output, sep="\t", index=False)
    logging.info("Done writing output to %s", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
