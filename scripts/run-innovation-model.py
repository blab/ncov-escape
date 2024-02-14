import argparse
import pandas as pd
import evofr as ef

LOCATIONS = ["USA"]
CI_COVERAGE = [0.8]

ITERS = 50_000
LEARNING_RATE = 4e-3
NUM_SAMPLES = 1000

def main():
    parser = argparse.ArgumentParser(description = "Given input sequence counts, estimate parent-child relative fitness differences along branches.")
    parser.add_argument("--seq-counts", type=str, required=True, help="output TSV of collapsed sequence counts")
    parser.add_argument("--pango-relationships", type=str, required=True, help="output TSV of pango-variant-relationships")
    parser.add_argument("--growth-advantage-path", type=str, required=True, help="output path for the estimated growth advantages by location")
    parser.add_argument("--pivot", type=str, required=True, help="The variant with which to be consider all advantages relative to.")
    args = parser.parse_args()

    # Load data 
    raw_seq = pd.read_csv(args.seq_counts, sep="\t")
    raw_variant_parents = pd.read_csv(args.pango_relationships, sep="\t")

    def _get_growth_advantage(location, pivot):
        # Filtering to location of interest
        _raw_seq = raw_seq[raw_seq.location == location].copy()
        _raw_variant_parents = raw_variant_parents[raw_variant_parents.location == location].copy()

        data = ef.InnovationSequenceCounts(_raw_seq, _raw_variant_parents, pivot=pivot)

        # Defining model
        model = ef.InnovationMLR(tau=4.2)

        # Defining inference method
        inference_method = ef.InferFullRank(ITERS, LEARNING_RATE, NUM_SAMPLES)

        # Fitting model
        posterior = inference_method.fit(model, data)
        ga_df = pd.DataFrame(
            ef.posterior.get_growth_advantage(
                posterior.samples,
                posterior.data,
                name=location,
                ps=CI_COVERAGE,
                rel_to=pivot)
            )
        return ga_df

    ga_dfs = [_get_growth_advantage(location, pivot=args.pivot) for location in LOCATIONS]
    ga_df = pd.concat(ga_dfs)
    ga_df.to_csv(args.growth_advantage_path, sep="\t")
    return None

if __name__ == "__main__":
    main()