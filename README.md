# Analysis of SARS-CoV-2 antigenic drift and evolutionary fitness

## Sequence counts and MLR analyis

### Provisioning data locally

Provision data locally by cloning https://github.com/nextstrain/forecasts-ncov/ and then running
```
nextstrain build . -j 1 --configfile config/config.yaml --config data_provenances=gisaid variant_classification=pango_lineages geo_resolutions=global
```
This will create `forecasts-ncov/data/gisaid/pango_lineages/global.tsv.gz` locally.

Go to the `ncov-escape` directory and run
```
mkdir -p data/gisaid/pango_lineages/
```
Copy the `global.tsv.gz` file from `forecasts-ncov` to `ncov-escape` in the corresponding directory.

### Sequence counts

Construct collapsed Pango sequence counts for analysis period and Pango relationships with
```
nextstrain build . -j 1 data/gisaid/pango_lineages/global/xbb15/pango_variant_relationships.tsv
```
