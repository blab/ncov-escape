# Analysis of SARS-CoV-2 antigenic drift and evolutionary fitness

## Nextstrain workflow

There is a Nextstrain workflow forked from v12 of the canonical Nextstrain [ncov](https://github.com/nextstrain/ncov) workflow in the [`ncov-workflow`](ncov-workflow/) directory.

## Sequence counts and MLR analyis

### Provisioning data locally

Provision data locally by cloning https://github.com/nextstrain/forecasts-ncov/ and then running
```
nextstrain build . -j 1 --configfile config/config.yaml --config data_provenances=gisaid variant_classification=pango_lineages geo_resolutions=global
```
This will create `forecasts-ncov/data/cases/global.tsv.gz` and `forecasts-ncov/data/gisaid/pango_lineages/global.tsv.gz` locally.

Go to the `ncov-escape` directory and run
```
mkdir -p data/cases/
mkdir -p data/gisaid/pango_lineages/
```
Copy each `global.tsv.gz` file from `forecasts-ncov` to `ncov-escape` in the respective directories.
