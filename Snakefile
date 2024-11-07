import os
import pandas as pd

if not config:
    configfile: "config/config.yaml"

wildcard_constraints:
    data_provenance="[A-Za-z0-9_-]+",  # Allow letters, numbers, underscores, and dashes
    analysis_period="[A-Za-z0-9_-]+",  # Allow letters, numbers, underscores, and dashes
    geo_resolution="[A-Za-z0-9_-]+",  # Allow letters, numbers, underscores, and dashes
    variant_classifications="[A-Za-z0-9_-]+"  # Allow letters, numbers, underscores, and dashes

def _get_date_range(ap):
    obs_date_min = config["analysis_period"].get(ap, {}).get('obs_date_min')
    obs_date_max = config["analysis_period"].get(ap, {}).get('obs_date_max')
    obs_date_interval = config["analysis_period"].get(ap, {}).get('interval')
    date_range = pd.date_range(
        start=obs_date_min,
        end=obs_date_max,
        freq=obs_date_interval,
    ).strftime("%Y-%m-%d").tolist()

    return date_range

def _get_all_input(w):
    data_provenances = config["data_provenances"] if isinstance(config["data_provenances"], list) else [config["data_provenances"]]
    variant_classifications = config["variant_classifications"] if isinstance(config["variant_classifications"], list) else [config["variant_classifications"]]
    geo_resolutions = config["geo_resolutions"] if isinstance(config["geo_resolutions"], list) else [config["geo_resolutions"]]
    analysis_periods = config["analysis_period"]

    # date_ranges = {analysis: generate_dates(analysis_periods[analysis]) for analysis in analysis_periods}
    all_input = [
        # Non-windowed analyses
        *expand(
            "results/{analysis_period}/growth_advantages.tsv",
            analysis_period=[ap for ap, cfg in analysis_periods.items() if not cfg.get("windowed", False)],
        )
    ]
    # Windowed analyses
    for ap, cfg in analysis_periods.items():
        if cfg.get("windowed", False):
            all_input += expand(
                "results/{analysis_period}/growth_advantages_{obs_date}.tsv",
                analysis_period=[ap],
                obs_date=_get_date_range(ap)
                )
    return all_input

rule all:
    input: _get_all_input

include: "workflow/snakemake_rules/prepare_data.smk"
include: "workflow/snakemake_rules/retrieve_phenotypes.smk"
include: "workflow/snakemake_rules/run_models.smk"
include: "workflow/snakemake_rules/run_models_over_period.smk"
