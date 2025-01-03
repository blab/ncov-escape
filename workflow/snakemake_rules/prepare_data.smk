# TODO: Switch this to load metadata and use `create-observed-sequence-counts`
get_analysis_config = lambda wildcards: config.get("analysis_period", {}).get(wildcards.analysis_period, {})

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


def _get_prepare_data_option_analysis(wildcards, option_name):
    """
    Return the option for prepare data from the config based on the
    wildcards.data_provenance, wildcards.variant_classification and the wildcards.geo_resolution values.

    If the *option* exists as a key within config['prepare_data'][wildcard.data_provenance][wildcard.geo_resolution]
    then return as "--{option-name} {option_value}". Or else return an empty string.
    """
    data_provenance = get_analysis_config(wildcards).get("data_provenance", "gisaid")
    variant_classification = get_analysis_config(wildcards).get("variant_classification", "pango_lineages")
    geo_resolution = get_analysis_config(wildcards).get("geo_resolution", "global"),
    option_value = config.get('prepare_data', {}) \
                         .get(data_provenance, {}) \
                         .get(variant_classification, {}) \
                         .get(geo_resolution, {}) \
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
    option_value = config.get('analysis_period', {}) \
                         .get(wildcards.analysis_period, {}) \
                         .get(option_name)

    if option_value is not None:
        # Change underscores of YAML keys to dashes for proper CLI option names
        option_name = option_name.replace('_', '-')
        return f'--{option_name} {option_value}'
    return ''


rule provision_sequence_counts:
    output:
        "data/{data_provenance}/{variant_classification}/{geo_resolution}.tsv.gz"
    shell:
        """
        python ./scripts/provision-files.py \
            --data-provenance {wildcards.data_provenance} \
            --variant-classification {wildcards.variant_classification} \
            --geo-resolution {wildcards.geo_resolution} \
            --output-path {output}
        """

rule prepare_clade_data:
    "Preparing clade counts for analysis"
    input:
        sequence_counts = lambda wildcards: "data/{data_provenance}/{variant_classification}/{geo_resolution}.tsv.gz".format(
            data_provenance= get_analysis_config(wildcards).get("data_provenance", "gisaid"),
            variant_classification=get_analysis_config(wildcards).get("variant_classification", "pango_lineages"),
            geo_resolution=get_analysis_config(wildcards).get("geo_resolution", "global")
        )
    output:
        sequence_counts = "data/{analysis_period}/prepared_seq_counts.tsv"
    log:
        "logs/{analysis_period}/prepare_clade_data.txt"
    params:
        min_date = lambda wildcards: _get_analysis_period_option(wildcards, 'min_date'),
        max_date = lambda wildcards: _get_analysis_period_option(wildcards, 'max_date'),
        location_min_seq = lambda wildcards: _get_prepare_data_option_analysis(wildcards, 'location_min_seq'),
        excluded_locations = lambda wildcards: _get_prepare_data_option_analysis(wildcards, 'excluded_locations'),
        prune_seq_days = lambda wildcards: _get_prepare_data_option_analysis(wildcards, 'prune_seq_days'),
        clade_min_seq = lambda wildcards: _get_prepare_data_option_analysis(wildcards, 'clade_min_seq'),
        force_include_clades = lambda wildcards: _get_prepare_data_option_analysis(wildcards, 'force_include_clades'),
        force_exclude_clades = lambda wildcards: _get_prepare_data_option_analysis(wildcards, 'force_exclude_clades')
    shell:
        """
        python ./scripts/prepare-data.py \
            --seq-counts {input.sequence_counts} \
            {params.min_date} \
            {params.max_date} \
            {params.location_min_seq} \
            {params.excluded_locations} \
            {params.clade_min_seq} \
            {params.force_include_clades} \
            {params.force_exclude_clades} \
            --output-seq-counts {output.sequence_counts} 2>&1 | tee {log}
        """

rule collapse_sequence_counts:
    "Collapsing Pango lineages, based on sequence count threshold"
    input:
        sequence_counts = "data/{analysis_period}/prepared_seq_counts.tsv"
    output:
        collapsed_counts = "data/{analysis_period}/collapsed_seq_counts.tsv"
    params:
        collapse_threshold = lambda wildcards: _get_prepare_data_option_analysis(wildcards, 'collapse_threshold'),
        force_include = lambda wildcards: _get_analysis_period_option(wildcards, 'force_include')
    shell:
        """
        python ./scripts/collapse-lineage-counts.py \
            --seq-counts {input.sequence_counts} \
            --output-seq-counts {output.collapsed_counts} \
            {params.collapse_threshold} \
            {params.force_include}
        """

rule get_pango_relationships:
    input:
        sequence_counts = "data/{analysis_period}/collapsed_seq_counts.tsv"
    output:
        pango_relationships = "data/{analysis_period}/pango_variant_relationships.tsv"
    shell:
        """
        python ./scripts/prepare-pango-relationships.py \
            --seq-counts {input.sequence_counts} \
            --output-relationships {output.pango_relationships}
        """
