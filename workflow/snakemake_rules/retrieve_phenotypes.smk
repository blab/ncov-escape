# TODO: Pull from main brach later once merged
rule download_phenotypes:
    output:
        "data/{analysis_period}/mutation_phenotypes.csv",
        "data/{analysis_period}/mutation_phenotypes_randomized.csv",
        "data/{analysis_period}/lineage_phenotypes.csv",
        "data/{analysis_period}/lineage_phenotypes_randomized.csv"
    shell:
        """
        python provision-phenotypes.py \
            --config-path XBB_config.yaml \
            --base-url https://raw.githubusercontent.com/jbloomlab/SARS2-spike-predictor-phenos/pull/1/files
        """
def _get_predictor_names(wildcards):
    predictor_names = config.get(wildcards).get('analysis_period', {}) \
                         .get(wildcards.analysis_period, {}) \
                         .get("predictor_names", None)
    return json.dumps(predictor_names)