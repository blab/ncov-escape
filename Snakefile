import os
import pandas as pd

if not config:
    configfile: "config/config.yaml"

def generate_dates(analysis_period):
    dates = pd.date_range(
        start=analysis_period['obs_date_min'],
        end=analysis_period['obs_date_max'],
        freq=analysis_period['obs_date_interval']
    )
    return dates.strftime('%Y-%m-%d').tolist()

def _get_all_input(w):
    data_provenances = config["data_provenances"] if isinstance(config["data_provenances"], list) else [config["data_provenances"]]
    variant_classifications = config["variant_classifications"] if isinstance(config["variant_classifications"], list) else [config["variant_classifications"]]
    geo_resolutions = config["geo_resolutions"] if isinstance(config["geo_resolutions"], list) else [config["geo_resolutions"]]
    analysis_periods = config["analysis_period"]

    date_ranges = {analysis: generate_dates(analysis_periods[analysis]) for analysis in analysis_periods}
    obs_date = [date for analysis in analysis_periods.keys() for date in date_ranges[analysis]]

    all_input = [
        *expand(
            "data/{data_provenance}/{variant_classification}/{geo_resolution}/{analysis_period}/prepared_cases.tsv",
            data_provenance=data_provenances,
            variant_classification=variant_classifications,
            geo_resolution=geo_resolutions,
            analysis_period=analysis_periods.keys()
        ),
        *expand(
            "data/{data_provenance}/{variant_classification}/{geo_resolution}/{analysis_period}/{obs_date}/collapsed_seq_counts_{obs_date}.tsv",
            data_provenance=data_provenances,
            variant_classification=variant_classifications,
            geo_resolution=geo_resolutions,
            analysis_period=analysis_periods.keys(),
            obs_date=obs_date
        ),
        *expand(
            "results/{data_provenance}/{variant_classification}/{geo_resolution}/{analysis_period}/{obs_date}/growth_advantages.tsv",
            data_provenance=data_provenances,
            variant_classification=variant_classifications,
            geo_resolution=geo_resolutions,
            analysis_period=analysis_periods.keys(),
            obs_date=obs_date
        ),
    ]
    return all_input


rule all:
    input: _get_all_input

include: "workflow/snakemake_rules/prepare_data.smk"
include: "workflow/snakemake_rules/run_models.smk"
