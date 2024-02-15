rule download_case_counts:
    output:
        cases = "data/cases/{geo_resolution}.tsv.gz"
    params:
        cases_url = "https://data.nextstrain.org/files/workflows/forecasts-ncov/cases/{geo_resolution}.tsv.gz"
    shell:
        """
        curl -fsSL --compressed {params.cases_url:q} --output {output.cases}
        """

rule download_sequence_counts:
    output:
        clades = "data/{data_provenance}/{variant_classification}/{geo_resolution}.tsv.gz"
    params:
        clades_url = "https://data.nextstrain.org/files/workflows/forecasts-ncov/data/{data_provenance}/{variant_classification}/global.tsv.gz"
    shell:
        """
        curl -fsSL --compressed {params.clades_url:q} --output {output.clades}
        """

def _get_prepare_data_option(wildcards, option_name):
    """
    Return the option for prepare data from the config based on the
    wildcards.data_provenance, wildcards.variant_classification and the wildcards.geo_resolution values.

    If the *option* exists as a key within config['prepare_data'][wildcard.data_provenance][wildcard.geo_resolution]
    then return as "--{option-name} {option_value}". Or else return an empty string.
    """
    option_value = config.get('prepare_data', {}) \
                         .get(wildcards.data_provenance, {}) \
                         .get(wildcards.variant_classification, {}) \
                         .get(wildcards.geo_resolution, {}) \
                         .get(option_name)


    if option_value is not None:
        # Change underscores of YAML keys to dashes for proper CLI option names
        option_name = option_name.replace('_', '-')
        return f'--{option_name} {option_value}'
    return ''

def _get_analysis_period_option(wildcards, option_name):
    """
    Return the option for analysis period from the config based on the analysis period values.

    If the *option* exists as a key within config['analysis_periods'][wildcard.analysis_period]
    then return as "--{option-name} {option_value}". Or else return an empty string.
    """
    option_value = config.get('analysis_periods', {}) \
                         .get(wildcards.analysis_periods, {}) \
                         .get(option_name)

    if option_value is not None:
        # Change underscores of YAML keys to dashes for proper CLI option names
        option_name = option_name.replace('_', '-')
        return f'--{option_name} {option_value}'
    return ''

rule prepare_clade_data:
    "Preparing clade counts for analysis"
    input:
        cases = "data/cases/{geo_resolution}.tsv.gz",
        sequence_counts = "data/{data_provenance}/{variant_classification}/{geo_resolution}.tsv.gz"
    output:
        cases = "data/{data_provenance}/{variant_classification}/{geo_resolution}/{period}/prepared_cases.tsv",
        sequence_counts = "data/{data_provenance}/{variant_classification}/{geo_resolution}/{period}/prepared_seq_counts.tsv"
    params:
        max_date = lambda wildcards: _get_analysis_period_option(wildcards, 'max_date'),
        included_days = lambda wildcards: _get_analysis_period_option(wildcards, 'included_days'),
        location_min_seq = lambda wildcards: _get_prepare_data_option(wildcards, 'location_min_seq'),
        location_min_seq_days = lambda wildcards: _get_prepare_data_option(wildcards, 'location_min_seq_days'),
        excluded_locations = lambda wildcards: _get_prepare_data_option(wildcards, 'excluded_locations'),
        prune_seq_days = lambda wildcards: _get_prepare_data_option(wildcards, 'prune_seq_days'),
        clade_min_seq = lambda wildcards: _get_prepare_data_option(wildcards, 'clade_min_seq'),
        clade_min_seq_days = lambda wildcards: _get_prepare_data_option(wildcards, 'clade_min_seq_days'),
        force_include_clades = lambda wildcards: _get_prepare_data_option(wildcards, 'force_include_clades')
    shell:
        """
        python ./scripts/prepare-data.py \
            --seq-counts {input.sequence_counts} \
            --cases {input.cases} \
            {params.max_date} \
            {params.included_days} \
            {params.location_min_seq} \
            {params.location_min_seq_days} \
            {params.excluded_locations} \
            {params.prune_seq_days} \
            {params.clade_min_seq} \
            {params.clade_min_seq_days} \
            {params.force_include_clades} \
            --output-seq-counts {output.sequence_counts} \
            --output-cases {output.cases}
        """

rule collapse_sequence_counts:
    "Collapsing Pango lineages, based on sequence count threshold"
    input:
        sequence_counts = "data/{data_provenance}/{variant_classification}/{geo_resolution}/{period}/prepared_seq_counts.tsv"
    output:
        collapsed_counts = "data/{data_provenance}/{variant_classification}/{geo_resolution}/{period}/collapsed_seq_counts.tsv"
    params:
        collapse_threshold = lambda wildcards: _get_prepare_data_option(wildcards, 'collapse_threshold')
    shell:
        """
        python ./scripts/collapse-lineage-counts.py \
            --seq-counts {input.sequence_counts} \
            --output-seq-counts {output.collapsed_counts} \
            {params.collapse_threshold}
        """

rule get_pango_relationships:
    input:
        sequence_counts = "data/{data_provenance}/{variant_classification}/{geo_resolution}/{period}/collapsed_seq_counts.tsv"
    output:
        pango_relationships = "data/{data_provenance}/{variant_classification}/{geo_resolution}/{period}/pango_variant_relationships.tsv"
    shell:
        """
        python ./scripts/prepare-pango-relationships.py \
            --seq-counts {input.sequence_counts} \
            --output-relationships {output.pango_relationships}
        """

# sets of measurements to compare to natural evolution
#TODO: Need rule to obtain these from repo?
phenos_compare_natural = {
    "current_dms": {
        "input_data": "results/summaries/summary.csv",
        "rename_cols": {
            "human sera escape": "sera escape",
            "spike mediated entry": "cell entry",
        },
        "phenotype_cols": ["sera escape", "ACE2 binding", "cell entry"],
        "title": "XBB.1.5 full-spike DMS phenotypes",
        "missing_muts": "drop",  # drop clades with missing muts
    },
    "yeast_RBD_DMS": {
        "input_data": "data/compare_natural_datasets/yeast_RBD_DMS.csv",
        "rename_cols": {},
        "phenotype_cols": ["escape", "ACE2 affinity", "RBD expression"],
        "title": "yeast RBD DMS phenotypes",
        "missing_muts": "zero",  # set missing (non-RBD) mutations to zero
    },
    "muts_from_Wuhan-Hu-1": {
        "input_data": "data/compare_natural_datasets/incremental_Hamming_distance_from_Wuhan-Hu-1.csv",
        "rename_cols": {"incremental Hamming distance": "distance"},
        "phenotype_cols": ["distance"],
        "title": "relative distance from Wuhan-Hu-1",
        "missing_muts": "drop",  # drop clades with missing muts
    },
    "EVEscape": {
        "input_data": "data/compare_natural_datasets/EVEscape_XBB_single_mutation_predictions.csv",
        "rename_cols": {},
        "phenotype_cols": ["EVEscape"],
        "title": "EVEscape",
        "missing_muts": "drop",  # drop clades with missing muts
    },
    "EVEscape_components": {
        "input_data": "data/compare_natural_datasets/EVEscape_XBB_single_mutation_predictions.csv",
        "rename_cols": {"fitness_evol_indices": "EVE fitness", "dissimilarity_charge_hydrophobicity": "aa dissimilarity", "accessibility_wcn": "accessibility"},
        "phenotype_cols": ["EVE fitness", "aa dissimilarity", "accessibility"],
        "title": "EVEscape components",
        "missing_muts": "drop",  # drop clades with missing muts
    },
}

# TODO: Want to be able to repeat for BA.2 pivot? This will ideally just be extending phenotypes and making it only run for particular values

def pass_phenotype_config(wildcards):
    config = yaml.round_trip_dump({
                "starting_clades": ["XBB"],  # clades descended from this
                "exclude_muts": [],  # exclude clades w these mutations
                "split_by_rbd": False,  # whether to treat RBD and non-RBD mutations separately
                "dms_clade": "XBB.1.5",  # clade used for DMS
                # rename columns in input data
                "rename_cols": phenos_compare_natural[wildcards.pheno]["rename_cols"],
                # "basic" means not split by RBD, which is done later in code if `split_by_rbd`
                "phenotype_cols": phenos_compare_natural[wildcards.pheno]["phenotype_cols"],
                # "drop" clades with missing mutations, or set missing mutations to "zero"
                "missing_muts": phenos_compare_natural[wildcards.pheno]["missing_muts"],
                "exclude_clades": [],  # exclude these clades
        })
    return config

rule compute_phenotypes:
    input:
        input_data = lambda wildcards: phenos_compare_natural.get(wildcards.phenotype), # Q. Where is wc.phenotype defined?
        pango_consensus_jsons= "data/placeholder.txt"
    output:
        clade_pair_dms="predictors/{pheno}_clade_pair.csv",
        clade_dms="predictors/{pheno}_clade.csv"
    params:
        # Q: How do I properly process this to pass the entire config? 
        # Would it be better just to pass things one by one? How would this work for dictionaries?
       config = pass_phenotype_config
    shell:
        """
        python ./scripts/.compute-phenotypes \
            --input-data {input.input_data} \
            --pango-consensus-jsons {input.pango_consensus_jsons} \
            --clade-pair-dms {output.clade_pair_dms} \
            --clade-dms {output.clade_dms} \
            --config '{params.yaml}'
        """