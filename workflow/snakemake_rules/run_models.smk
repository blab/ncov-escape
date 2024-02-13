
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

rule innovation_model:
    input:
        sequence_counts = "data/{data_provenance}/{variant_classification}/{geo_resolution}/{period}/collapsed_seq_counts.tsv"
        pango_relationships = "data/{data_provenance}/{variant_classification}/{geo_resolution}/{period}/pango_variant_relationships.tsv"
    params:
        min_date = lambda wildcards: _get_analysis_period_option(wildcards, 'min_date'),
        max_date = lambda wildcards: _get_analysis_period_option(wildcards, 'max_date'),
        pivot = lambda wildcards: _get_analysis_period_option(wildcards, 'pivot')
    output:
        results: "results/{data_provenance}/{variant_classification}/{geo_resolution}/{period}/growth_advantages.tsv"
    shell:
        """
        python ./scripts/run-innovation-model.py \
            --seq-counts {input.sequence_counts} \
            --pango-relationships {input.pango_relationships} \
            --growth-advantage-path {out.results} \
            {params.min_date} \
            {params.max_date} \
            {params.pivot}
        """
