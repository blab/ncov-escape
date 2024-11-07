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
        sequence_counts = "data/{analysis_period}/collapsed_seq_counts.tsv",
        pango_relationships = "data/{analysis_period}/pango_variant_relationships.tsv",
    params:
        pivot = lambda wildcards: _get_analysis_period_option(wildcards, 'pivot'),
    	posteriors = "results/{analysis_period}/posteriors"
    output:
        growth_advantages = "results/{analysis_period}/growth_advantages.tsv",
        growth_advantages_delta = "results/{analysis_period}/growth_advantages_delta.tsv"
    shell:
        """
        python ./scripts/run-innovation-model.py \
            --seq-counts {input.sequence_counts} \
            --pango-relationships {input.pango_relationships} \
            --growth-advantage-path {output.growth_advantages} \
            --growth-advantage-delta-path {output.growth_advantages_delta} \
            --posterior-path {params.posteriors} \
            {params.pivot}
        """

rule innovation_model_informed:
    input:
        sequence_counts = "data/{analysis_period}/collapsed_seq_counts.tsv",
        pango_relationships = "data/{analysis_period}/pango_variant_relationships.tsv",
        predictor_path = "data/{analysis_period}/lineage_phenotypes.csv"
    params:
        predictor_names = lambda wildcards: _get_predictor_names(wildcards)
        pivot = lambda wildcards: _get_analysis_period_option(wildcards, 'pivot'),
    	posteriors = "results/{analysis_period}/posteriors/informed"
    output:
        growth_advantages = "results/{analysis_period}/informed/growth_advantages.tsv",
        growth_advantages_delta = "results/{analysis_period}/informed/growth_advantages_delta.tsv"
    shell:
        """
        python ./scripts/run-innovation-model.py \
            --seq-counts {input.sequence_counts} \
            --pango-relationships {input.pango_relationships} \
            --predictor-path {input.predictor_path} \
            --predictor-names {params.predictor_names} \
            --growth-advantage-path {output.growth_advantages} \
            --growth-advantage-delta-path {output.growth_advantages_delta} \
            --posterior-path {params.posteriors} \
            {params.pivot}
        """
