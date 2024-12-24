"""
Summarize sequence counts grouped by date, submission_date, location, and clade.
"""
import argparse
import os
import sys
from datetime import datetime

import pandas as pd


def format_date(date_string, expected_format):
    """
    Format *date_string* to ISO 8601 date (YYYY-MM-DD).
    If *date_string* does not match *expected_format*, return None.

    >>> expected_format = '%Y-%m-%d'
    >>> format_date("2020", expected_format)
    None
    >>> format_date("2020-01", expected_format)
    None
    >>> format_date("XXXX-XX-XX", expected_format)
    None
    >>> format_date("2020-1-15", expected_format)
    '2020-01-15'
    >>> format_date("2020-01-15", expected_format)
    '2020-01-15'
    """
    try:
        return datetime.strptime(date_string, expected_format).strftime("%Y-%m-%d")
    except (TypeError, ValueError):
        return None


def count_sequences_with_submission_date(metadata):
    grouped = metadata.groupby(["date", "location", "variant", "date_submitted"])
    counts = grouped.size().reset_index(name="sequences")
    return counts


BIAS_BUFFER = 14

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        __doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--metadata",
        required=True,
        help="Path to metadata TSV that includes date, location, and clade information. "
        + "Local and remote (HTTP, S3, etc.) paths are valid.",
    )
    parser.add_argument(
        "--id-column",
        default="strain",
        help="Column in metadata TSV with the sequence ID.",
    )
    parser.add_argument(
        "--date-column",
        default="date",
        help="Column in metadata TSV with date information. "
        + "Dates are expected to be in ISO 8601 date format (i.e. YYYY-MM-DD). "
        + "Rows with incorrect date format or ambiguous dates will be dropped.",
    )
    parser.add_argument(
        "--location-column",
        default="country",
        help="Column in metadata TSV with location information",
    )
    parser.add_argument(
        "--clade-column",
        required=False,
        help="Column in metadata TSV with clade information",
    )
    parser.add_argument(
        "--filter-columns",
        nargs="+",
        help="Columns that will be used in the `--filter-query` option. "
        + "Must be provided if using the `--filter-query` option.",
    )
    parser.add_argument(
        "--filter-query",
        help="Filter sequences by attribute. "
        + "Uses Pandas Dataframe querying, "
        + "see https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query for syntax "
        + """(e.g., --filter-query "country == 'USA'") """,
    )
    parser.add_argument(
        "--metadata-chunk-size",
        type=int,
        default=500_000,
        help="Maximum metadata records to read into memory at once during initial pass."
        + "Increasing this value increases peak memory usage.",
    )
    parser.add_argument("--obs-date-min", help="First date of observation.")
    parser.add_argument("--obs-date-max", help="Last date of observation.")
    parser.add_argument(
        "--num-days-context",
        type=int,
        help="Optionally, the number of days that are included in context"
        + "By default, this includes all days after `obs_date_min`.",
    )
    parser.add_argument(
        "--output-path",
        help="Path to output TSV for sequence counts by observation date.",
    )

    args = parser.parse_args()

    if args.filter_query and not args.filter_columns:
        print(
            "ERROR: Filter columns must be provided if using a filter query.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Map of metadata columns to output columns
    metadata_column_map = {
        args.id_column: "sequences",
        args.date_column: "date",
        args.location_column: "location",
        args.clade_column: "variant",
        "date_submitted": "date_submitted",
    }

    # Only use required columns, adding filter columns if provided
    metadata_usecols = set(metadata_column_map.keys())
    if args.filter_columns:
        metadata_usecols.update(args.filter_columns)

    # Load metadata TSV.
    metadata_reader = pd.read_csv(
        args.metadata,
        sep="\t",
        usecols=list(metadata_usecols),
        dtype="object",
        chunksize=args.metadata_chunk_size,
    )

    # Iterate through metadata in chunks to control peak memory usage.
    print("Chunking")
    metadata_chunks = []
    for metadata in metadata_reader:
        # If provided filter query, apply query then subset to required columns
        if args.filter_query:
            try:
                metadata.query(args.filter_query, inplace=True)
            except Exception as e:
                print(
                    "ERROR: An error occurred when applying the filter query. "
                    "Most likely the filter query used columns that were not included in the filter columns. "
                    f"See detailed error: ({e})",
                    file=sys.stderr,
                )
                sys.exit(1)

            metadata = metadata[metadata_column_map.keys()]

        # Rename columns to output column names
        metadata.rename(columns=metadata_column_map, inplace=True)

        # Convert location and clade columns to category dtype
        metadata["location"] = metadata["location"].astype("category")
        metadata["variant"] = metadata["variant"].astype("category")

        # Convert date column to datetime, sets ambiguous dates to None
        metadata["date"] = pd.to_datetime(metadata["date"], errors="coerce")
        metadata.dropna(subset=["date"], inplace=True)

        # Filter to time period
        min_date = pd.Timestamp(args.obs_date_min) - pd.Timedelta(
            days=args.num_days_context + BIAS_BUFFER
        )
        max_date = pd.Timestamp(args.obs_date_max)
        metadata = metadata[(metadata.date >= min_date) & (metadata.date <= max_date)]
        metadata["date"] = pd.to_datetime(metadata["date"]).dt.strftime("%Y-%m-%d")
        metadata["date_submitted"] = pd.to_datetime(
            metadata["date_submitted"]
        ).dt.strftime("%Y-%m-%d")

        # Drop rows with null date, location, or variants
        metadata = metadata.dropna(subset=["date", "location", "variant"])
        metadata_chunks.append(metadata)

    # Merge all chunks that passed all filters.
    metadata = pd.concat(
        metadata_chunks,
        ignore_index=True,
    )

    print("Counting")

    sequence_count_by_submission = count_sequences_with_submission_date(metadata)

    # Make sure we have the folder
    path = args.output_path
    if not os.path.exists(path):
        os.makedirs(path)

    sequence_count_by_submission.to_csv(
        f"{path}/sequence_counts_by_submission.tsv", sep="\t", index=False
    )
