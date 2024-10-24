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

rule innovation_model:
    input:
        sequence_counts = "data/{data_provenance}/{variant_classification}/{geo_resolution}/{analysis_period}/{obs_date}/collapsed_seq_counts_{obs_date}.tsv",
        pango_relationships = "data/{data_provenance}/{variant_classification}/{geo_resolution}/{analysis_period}/{obs_date}/pango_variant_relationships_{obs_date}.tsv",
    params:
        pivot = lambda wildcards: _get_analysis_period_option(wildcards, 'pivot')
    output:
    	posteriors = "results/{data_provenance}/{variant_classification}/{geo_resolution}/{analysis_period}/{obs_date}/posteriors",
        growth_advantages = "results/{data_provenance}/{variant_classification}/{geo_resolution}/{analysis_period}/{obs_date}/growth_advantages.tsv"
    shell:
        """
        python ./scripts/run-innovation-model.py \
            --seq-counts {input.sequence_counts} \
            --pango-relationships {input.pango_relationships} \
            --growth-advantage-path {output.growth_advantages} \
            --posterior-path {output.posteriors} \
            {params.pivot}
        """
