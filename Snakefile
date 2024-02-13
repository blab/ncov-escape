import os

if not config:
    configfile: "config/config.yaml"

def _get_all_input(w):
    data_provenances = config["data_provenances"] if isinstance(config["data_provenances"], list) else [config["data_provenances"]]
    variant_classifications = config["variant_classifications"] if isinstance(config["variant_classifications"], list) else [config["variant_classifications"]]
    geo_resolutions = config["geo_resolutions"] if isinstance(config["geo_resolutions"], list) else [config["geo_resolutions"]]
    analysis_periods = config["analysis_periods"] if isinstance(config["analysis_periods"], list) else [config["analysis_periods"]]

    all_input = [
        # Prepared data sets
        # TODO: Remove cases?
        *expand(
            "data/{data_provenance}/{variant_classification}/{geo_resolution}/{period}/prepared_cases.tsv",
            data_provenance=data_provenances,
            variant_classification=variant_classifications,
            geo_resolution=geo_resolutions,
            period=analysis_periods
        ),
        *expand(
            "data/{data_provenance}/{variant_classification}/{geo_resolution}/{period}/prepared_seq_counts.tsv",
            data_provenance=data_provenances,
            variant_classification=variant_classifications,
            geo_resolution=geo_resolutions,
            period=analysis_periods
        ),
        # Resulting growth advantages
        *expand(
            "results/{data_provenance}/{variant_classification}/{geo_resolution}/{period}/growth_advantages.tsv",
            data_provenance=data_provenances,
            variant_classification=variant_classifications,
            geo_resolution=geo_resolutions,
            period=analysis_periods
        ),

    ]

rule all:
    input: _get_all_input

include: "workflow/snakemake_rule/prepare_data.smk"
include: "workflow/snakemake_rule/run_models.smk"