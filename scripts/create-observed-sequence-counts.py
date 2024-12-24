"""
Summarize sequence counts grouped by date, location, and clade.
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
    grouped = metadata.groupby(["date", "location", "variant"], group_keys=True)

    def compute_delay(x):
        delays = (
            x["date_submitted"]
            .value_counts(normalize=False)
            .rename_axis("date_submitted")
            .reset_index(name="sequences")
        )

        delays = delays.sort_values("date_submitted")
        return delays

    out = (
        grouped.apply(compute_delay)
        .reset_index(level=list(range(len(grouped.grouper.names))), drop=True)
        .reset_index()
    )
    return out


def observe_sequence_counts(delayed, obs_date=None):
    # Given an observation date as well as counts of sequences and their submission dates,
    # Reconstruct data available on observation date

    obs_seq = delayed.copy()

    # Filter to sequences submitted on or before date
    if obs_date:
        obs_seq = obs_seq[obs_seq["date_submitted"] < obs_date]

    # Sum across remaining sequences
    obs_seq = (
        obs_seq.groupby(["date", "location", "variant"])["sequences"].sum()
    ).reset_index()

    # Sort data
    obs_seq = obs_seq.sort_values(["location", "variant", "date"])

    # Remove entries with no observed sequences
    obs_seq = obs_seq[obs_seq["sequences"] > 0]

    return obs_seq


BIAS_BUFFER = 14


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        __doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--sequence-counts-by-submission",
        requred=True,
        help="Path to `sequence_counts_by_submission` TSV.",
    )
    parser.add_argument("--obs-date-min", help="First date of observation.")
    parser.add_argument("--obs-date-max", help="Last date of observation.")
    parser.add_argument(
        "--interval", help="The interval for the date range e.g. '2W' or '1D'."
    )
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

    sequence_count_by_submission = pd.read_csv(
        args.sequence_counts_by_submission, sep="\t"
    )

    # We need to create a directory to hold observed counts
    observation_dates = pd.date_range(
        start=args.obs_date_min, end=args.obs_date_max, freq=args.interval
    )

    for obs_date in observation_dates:
        # Observe sequences up to this date
        obs_date = obs_date.strftime("%Y-%m-%d")
        obs_seq = observe_sequence_counts(
            sequence_count_by_submission, obs_date=obs_date
        )

        # Filter to appropriate date, either to maintain consistent window or grow window in time
        if args.num_days_context is None:
            min_date = args.obs_date_min
        else:
            min_date = pd.to_datetime(obs_date) - pd.Timedelta(
                args.num_days_context + BIAS_BUFFER, "d"
            )
        obs_seq = obs_seq[obs_seq.date > min_date]

        # Remove most recent 14 days due to bias
        max_date = (
            pd.to_datetime(obs_date) - pd.Timedelta(BIAS_BUFFER, "d")
        ).dt.strftime("%Y-%m-%d")
        obs_seq = obs_seq[obs_seq.date <= max_date]

        # Export file
        path = args.output_path
        if not os.path.exists(path):
            os.makedirs(path)

        # Make sure we have the folder
        obs_seq.to_csv(
            f"{path}/prepared_seq_counts_{obs_date}.tsv", sep="\t", index=False
        )

    retrospective_seq_counts = observe_sequence_counts(
        sequence_count_by_submission, obs_date=None
    )

    # Make sure we have the folder
    path = args.output_path
    if not os.path.exists(path):
        os.makedirs(path)

    retrospective_seq_counts.to_csv(
        f"{path}/seq_counts_retrospective.tsv", sep="\t", index=False
    )
