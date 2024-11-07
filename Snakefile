import os
import pandas as pd

if not config:
    configfile: "config/config.yaml"

wildcard_constraints:
    data_provenance="[A-Za-z0-9_-]+",  # Allow letters, numbers, underscores, and dashes
    analysis_period="[A-Za-z0-9_-]+",  # Allow letters, numbers, underscores, and dashes
    geo_resolution="[A-Za-z0-9_-]+",  # Allow letters, numbers, underscores, and dashes
    variant_classifications="[A-Za-z0-9_-]+"  # Allow letters, numbers, underscores, and dashes


def _get_all_input(w):
    data_provenances = config["data_provenances"] if isinstance(config["data_provenances"], list) else [config["data_provenances"]]
    variant_classifications = config["variant_classifications"] if isinstance(config["variant_classifications"], list) else [config["variant_classifications"]]
    geo_resolutions = config["geo_resolutions"] if isinstance(config["geo_resolutions"], list) else [config["geo_resolutions"]]
    analysis_periods = config["analysis_period"]

    # date_ranges = {analysis: generate_dates(analysis_periods[analysis]) for analysis in analysis_periods}
    # obs_date = [date for analysis in analysis_periods.keys() for date in date_ranges[analysis]]

    all_input = [
        *expand(
            "results/{analysis_period}/growth_advantages.tsv",
            # data_provenance=data_provenances,
            # variant_classification=variant_classifications,
            # geo_resolution=geo_resolutions,
            analysis_period=analysis_periods.keys(),
            # obs_date=obs_date
        ),
    ]
    return all_input

rule all:
    input: _get_all_input

include: "workflow/snakemake_rules/prepare_data.smk"
include: "workflow/snakemake_rules/run_models.smk"
