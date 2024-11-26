# TODO: Pull from main brach later once merged
import json

BASE_URL = "https://raw.githubusercontent.com/jbloomlab/SARS2-spike-predictor-phenos/refs/heads/Get-Lineage-Phenotypes"
BASE_URL = "https://raw.githubusercontent.com/jbloomlab/SARS2-spike-predictor-phenos/refs/heads/main"

rule download_phenotypes:
    output:
        "data/{analysis_period}/phenotypes/mutation_phenotypes.csv",
        "data/{analysis_period}/phenotypes/mutation_phenotypes_randomized.csv",
        "data/{analysis_period}/phenotypes/lineage_phenotypes.csv",
        "data/{analysis_period}/phenotypes/lineage_phenotypes_randomized.csv"
    params:
        output_path = "data/{analysis_period}/phenotypes",
        base_url = BASE_URL
    shell:
        """
        python ./scripts/provision-phenotypes.py \
            --config-path config.yaml \
            --base-url {params.base_url} \
            --output-path {params.output_path}
        """
def _get_predictor_names(wildcards):
    predictor_names = config.get('analysis_period', {}) \
                        .get(wildcards.analysis_period, {}) \
                        .get("predictor_names", None)
    return f"'{json.dumps(predictor_names)}'"
