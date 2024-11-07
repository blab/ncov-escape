#TODO: Run collapsing logic on this

rule observed_over_period:
    input:
        # sequence_counts = "data/{data_provenance}/{variant_classification}/{geo_resolution}/{analysis_period}/collapsed_seq_counts.tsv"
        metadata = lambda wildcards: "data/{data_provenance}_metadata.tsv.gz".format(
            data_provenance=get_analysis_config(wildcards).get("data_provenance", "gisaid")
        )
    output:
        sequence_counts_dated = "data/{analysis_period}/prepared_seq_counts_{obs_date}.tsv", 
    params:
        obs_date_min = lambda wildcards: _get_analysis_period_option(wildcards, 'obs_date_min'),
        obs_date_max = lambda wildcards: _get_analysis_period_option(wildcards, 'obs_date_max'),
        obs_date_interval = lambda wildcards: _get_analysis_period_option(wildcards, 'interval'),
        num_days_context = lambda wildcards: _get_analysis_period_option(wildcards, 'num_days_context'),
        output_path = lambda wildcards: f"data/{wildcards.analysis_period}"
    shell:
        """
        python ./scripts/create-observed-sequence-counts.py \
            --metadata {input.metadata} \
            --clade-column "Nextclade_pango" \
            --output-path {params.output_path} \
            --filter-columns "QC_overall_status" \
            --filter-query "QC_overall_status != 'bad'" \
            {params.obs_date_min}\
            {params.obs_date_max} \
            {params.obs_date_interval} \
            {params.num_days_context}
        """

rule collapse_over_period:
    "Collapsing Pango lineages, based on sequence count threshold"
    input:
        sequence_counts = "data/{analysis_period}/prepared_seq_counts_{obs_date}.tsv"
    output:
        collapsed_counts = "data/{analysis_period}/collapsed_seq_counts_{obs_date}.tsv"
    params:
        collapse_threshold = lambda wildcards: _get_prepare_data_option_analysis(wildcards, 'collapse_threshold')
    shell:
        """
        python ./scripts/collapse-lineage-counts.py \
            --seq-counts {input.sequence_counts} \
            --output-seq-counts {output.collapsed_counts} \
            {params.collapse_threshold}
        """
    
rule run_innovation_model_over_period:
    input:
        sequence_counts = lambda wildcards: "data/{analysis_period}/collapsed_seq_counts_{obs_date}.tsv".format(
            analysis_period = wildcards.analysis_period,
            obs_date = wildcards.obs_date
            ),
        pango_relationships = "data/{analysis_period}/pango_variant_relationships.tsv",
    params:
        pivot = lambda wildcards: _get_analysis_period_option(wildcards, 'pivot'),
    	posteriors = lambda wildcards: "results/{analysis_period}/posteriors_{obs_date}".format(
            analysis_period = wildcards.analysis_period,
            obs_date = wildcards.obs_date
            ),
    output:
        growth_advantages = "results/{analysis_period}/growth_advantages_{obs_date}.tsv",
        growth_advantages_delta = "results/{analysis_period}/growth_advantages_delta_{obs_date}.tsv"
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

rule run_innovation_model_informed_over_period:
    input:
        sequence_counts = lambda wildcards: expand(
            "data/{analysis_period}/collapsed_seq_counts_{obs_date}.tsv",
            analysis_period = wildcards.analysis_period,
            obs_date = wildcards.obs_date
        ),
        pango_relationships = "data/{analysis_period}/pango_variant_relationships.tsv",
        predictor_path = "data/{analysis_period}/lineage_phenotypes.csv",
    params:
        predictor_names = lambda wildcards: _get_predictor_names(wildcards),
        pivot = lambda wildcards: _get_analysis_period_option(wildcards, 'pivot'),
    	posteriors = "results/{analysis_period}/posteriors_{obs_date}/informed"
    output:
        growth_advantages = "results/{analysis_period}/informed/growth_advantages_{obs_date}.tsv",
        growth_advantages_delta = "results/{analysis_period}/informed/growth_advantages_delta_{obs_date}.tsv"
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