import argparse

import evofr as ef
import numpy as np
import pandas as pd

LOCATIONS = ["USA"]
CI_COVERAGE = [0.8]

ITERS = 50_000
LEARNING_RATE = 4e-3
NUM_SAMPLES = 1000


def get_growth_advantage(posterior, pivot):
    ga_df = pd.DataFrame(
        ef.posterior.get_growth_advantage(
            posterior.samples,
            posterior.data,
            name=posterior.name,
            ps=CI_COVERAGE,
            rel_to=pivot,
        )
    )
    return ga_df


def prep_predictors(predictors, variant_freqs, predictor_names):
    # Index by variant name
    predictors = predictors.rename(columns={"seqName": "variant"}).set_index("variant")
    predictors = predictors.replace("?").astype(
        {name: "float" for name in predictor_names}
    )

    # Find scores of interest and parents
    var_names = [v for v in variant_freqs.var_names if v in predictors.index]
    predictors = predictors.loc[var_names]  # Need all variants to be present...
    predictors["parent"] = predictors.index.map(variant_freqs.parent_map)

    # Get delta between parents and children
    def get_parent_delta(x, col="immune_escape"):
        variant = x.name
        parent = x.parent
        # If parent and child are present generate contrast
        if parent in predictors.index:
            return predictors.loc[variant][col] - predictors.loc[parent][col]
        # Gotta figure out how to deal with the nans
        return np.nan

    # Generate delta columns
    for name in predictor_names:
        predictors[f"delta_{name}"] = predictors.apply(
            lambda x: get_parent_delta(x, name), axis=1
        )
    return predictors


def make_features(predictors, variant_freqs, feature_names=None, intercept=True):

    if feature_names is None:
        feature_names = ["delta_immune_escape", "delta_ace2_binding"]

    n_features = len(feature_names)

    # Fill with features from data frame
    N_variants = len(variant_freqs.var_names)
    features = np.empty((N_variants, n_features))

    for v, var in enumerate(variant_freqs.var_names):
        if var in predictors.index:
            features[v, :] = predictors.loc[var][feature_names].values
        else:
            features[v, :] = np.nan

    # Add intercept if desired
    if intercept:
        features = np.column_stack((features, np.ones(N_variants)))
    return features


# TODO: Load in ITERS, LEARNING_RATE, NUM_SAMPLES from a config file as in forecasts-ncov
# Similarly, load in CI_COVERAGE

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Given input sequence counts, estimate parent-child relative fitness differences along branches."
    )
    parser.add_argument(
        "--seq-counts",
        type=str,
        required=True,
        help="input TSV of collapsed sequence counts",
    )
    parser.add_argument(
        "--pango-relationships",
        type=str,
        required=True,
        help="input TSV of pango-variant-relationships",
    )
    parser.add_argument(
        "--predictor_path",
        type=str,
        default=None,
        help="input TSV of predictors of variant fitness",
    )
    parser.add_argument(
        "--growth-advantage-path",
        type=str,
        required=True,
        help="output path for the estimated growth advantages by location",
    )
    parser.add_argument(
        "--pivot",
        type=str,
        required=True,
        help="The variant with which to be consider all advantages relative to.",
    )
    parser.add_argument(
        "--posterior-path",
        type=str,
        default=None,
        help="The path to save the posteriors by location.",
    )
    args = parser.parse_args()

    # Load data
    raw_seq = pd.read_csv(args.seq_counts, sep="\t")
    raw_variant_parents = pd.read_csv(args.pango_relationships, sep="\t")

    # Use all location present unless instructed otherwise
    locations = pd.unique(raw_seq["location"])

    def _get_posterior(location, pivot):
        # Filtering to location of interest
        _raw_seq = raw_seq[raw_seq.location == location].copy()
        data = ef.InnovationSequenceCounts(_raw_seq, raw_variant_parents, pivot=pivot)

        # Defining model
        if args.predictor_path is None:
            model = ef.InnovationMLR(tau=4.2)
        else:
            # Define predictors
            predictors = pd.read_csv(args.predictor_path, sep="\t")
            # TODO: All config to specify predictors to use
            predictor_names = ["immune_escape", "ace2_binding"]
            predictors = prep_predictors(
                predictors, data, predictor_names=predictor_names
            )
            features = make_features(predictors, data, feature_names=predictor_names)
            prior = ef.models.DeltaRegressionPrior(features)
            model = ef.InnovationMLR(tau=4.2, delta_prior=prior)

        # Defining inference method
        inference_method = ef.InferFullRank(ITERS, LEARNING_RATE, NUM_SAMPLES)

        # Fitting model
        posterior = inference_method.fit(model, data, name=location)
        if args.posterior_path is not None:
            posterior.save_posterior(args.posterior_path + f"/{location}.pkl")
        return posterior

    posteriors = [_get_posterior(location, pivot=args.pivot) for location in locations]

    ga_dfs = [
        get_growth_advantage(posterior, pivot=args.pivot) for posterior in posteriors
    ]
    ga_df = pd.concat(ga_dfs)
    ga_df.to_csv(args.growth_advantage_path, sep="\t")
