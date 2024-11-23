# TODO: Pull from main brach later once merged
import json
rule download_phenotypes:
    output:
        "data/{analysis_period}/phenotypes/mutation_phenotypes.csv",
        "data/{analysis_period}/phenotypes/mutation_phenotypes_randomized.csv",
        "data/{analysis_period}/phenotypes/lineage_phenotypes.csv",
        "data/{analysis_period}/phenotypes/lineage_phenotypes_randomized.csv"
    params:
        output_path = "data/{analysis_period}/phenotypes"
    shell:
        """
        python ./scripts/provision-phenotypes.py \
            --config-path XBB_config.yaml \
            --base-url https://raw.githubusercontent.com/jbloomlab/SARS2-spike-predictor-phenos/refs/heads/Get-Lineage-Phenotypes \
            --output-path {params.output_path}
        """
def _get_predictor_names(wildcards):
    predictor_names = config.get('analysis_period', {}) \
                        .get(wildcards.analysis_period, {}) \
                        .get("predictor_names", None)
    return json.dumps(predictor_names)
