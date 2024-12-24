# Analysis of SARS-CoV-2 antigenic drift and evolutionary fitness

## Installation

All dependencies can be met by running through the Nextstrain runtime.
See [docs.nextstrain.org](https://docs.nextstrain.org/en/latest/install.html) for installation instructions.

## Provision metadata locally

Windowed analyses require access to SARS-CoV-2 metadata.
This can be acquired via

```bash
aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.gz data/gisaid_metadata.tsv.gz
zstd -c -d data/gisaid_metadata.tsv.gz \
   | tsv-select -H -f strain,date,date_submitted,country,clade_nextstrain,Nextclade_pango,QC_overall_status \
   | gzip -c > data/gisaid_metadata_filtered.tsv.gz
```

Access to this S3 bucket is restricted based on GISAID data sharing policies.

Non-windowed analyses will provision sequence counts from [forecasts-ncov](https://github.com/nextstrain/forecasts-ncov). These reduced sequence counts are publicly available.

Windowed analyses require working from detailed metadata as the rules [`process_metadata` and `observe_over_period`](/workflow/snakemake_rules/run_models_over_period.smk) restrict sequence counts to follow what was actually available based on submission dates, rather than what's available at the current moment.


## Workflow

Once metadata is provision locally, run the workflow with

```
nextstrain build . all
```

This will generate sequence count files for all analyses specificied in `./config/config.yaml`.
Under `analysis_period` in the config file, you can specify an analysis name as well as minimum and maximum dates for the analysis, the pivot variant of interest, lineages to forcibly include, and predictors to use for the regression-prior model e.g.

```
analysis_period:
  xbb15:
    min_date: "2023-01-01"
    max_date: "2023-12-01"
    pivot: "XBB.1.5"
    force_include: "defaults/xbb15/force_include_lineages.txt"
    predictor_names:
      - "spike pseudovirus DMS human sera escape relative to XBB.1.5"
      - "spike pseudovirus DMS ACE2 binding relative to XBB.1.5"
      - "RBD yeast-display DMS ACE2 affinity relative to XBB.1.5"
      - "RBD yeast-display DMS RBD expression relative to XBB.1.5"
      - "RBD yeast-display DMS escape relative to XBB.1.5"
```

### Sequence counts and variant relationships

It produces sequence count files for windowed analyses, non-windowed analyses, collapsing lineages into their parents based on count thresholds specified in `./config/config.yaml`, as well as generating variant relationship files, so that the innovation model can be fit.

Collapsed sequence counts follow the form `data/{analysis_period}/collapsed_seq_counts.tsv` and variant relationships follow the form `data/{analysis_period}/pango_variant_relationships.tsv`.

### MLR innovation model

The collapsed sequence count files and variant relationships are used to estimate relative fitness using the uninformed (normal-prior) innovation model.

For each analysis period, this produces `results/{analysis_period}/growth_advantage.tsv` and `results/{analysis_periods}/growth_advantage_delta.tsv`.

The model posteriors will also be saved under `results/{analysis_period}/posteriors/data_{location}.pkl` and `results/{analysis_period}/posteriors/samples_{location}.pkl`.

If predictor names are provided, they are used to fit the regression prior innovation model.
These predictor-informed results are stored under `results/{analysis_period}/informed/growth_advantages.tsv` and `results/{analysis_period}/informed/growth_advantages.tsv`.
