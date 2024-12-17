# Analysis of SARS-CoV-2 antigenic drift and evolutionary fitness

## Provision metadata locally

Windowed analyses require access to SARS-CoV-2 metadata.
This can be acquired via

```bash
aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.gz data/gisaid/gisaid_metadata.tsv.gz
zstd -c -d data/gisaid/giasaid_metadata.tsv.gz \
   | tsv-select -H -f strain,date,date_submitted,country,clade_nextstrain,Nextclade_pango,QC_overall_status \
   | zstd -c > data/gisaid/giasaid_metadata.tsv.gz
```

Non-windowed analyses will provision sequence counts from [forecasts-ncov](https://github.com/nextstrain/forecasts-ncov).


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

It will then produce sequence count files for windowed analyses, non-windowed analyses, collapsing lineages into their parents based on count thresholds specified in `./config/config.yaml`.

This will also generate variant relationship files, so that the innovation model can be fit.

The sequence count files and variant relationships will then be used to estimate relative fitness using the uninformed (normal-prior) innovation model.

If predictor names are provided, they will be also be used to fit the regression prior innovation model.

For each analysis, this will produce `results/{analysis_period}/growth_advantage.tsv` and `results/{analysis_periods}/growth_advantage_delta.tsv` for each analysis period.

The model posteriors will also be saved under `results/{analysis_period}/posteriors/data_{location}.pkl` and `results/{analysis_period}/posteriors/samples_{location}.pkl`.
