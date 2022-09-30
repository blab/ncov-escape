## Pango lineages dataset

This dataset only includes Pango lineages with >100 sequences in US dataset.

Data preparation followed:

1. Nextstrain-curated metadata TSV of GISAID database was downloaded. Uncompressing and renaming this file resulted in `gisaid_metadata.tsv` via:
```
nextstrain remote download s3://nextstrain-ncov-private/metadata.tsv.gz
gzip -d metadata.tsv.gz -c > gisaid_metadata.tsv
```

2. The metadata file was pruned to only relevant columns via:
```
tsv-select -H -f strain,date,country,division,QC_overall_status,Nextclade_pango gisaid_metadata.tsv > gisaid_metadata_pruned.tsv
```

3. This `gisaid_metadata_pruned.tsv` is processed in Mathematica by running the notebook `pango_data-prep.nb`. This results in the file `pango_location-variant-sequence-counts.tsv` versioned here. These files represent heavily derived GISAID data and are equivalent to downloadable results from [outbreak.info](https://outbreak.info), [cov-spectrum.org](https://cov-spectrum.org) and [covariants.org](https://covariants.org). This use is allowable under the [GISAID Terms of Use](https://www.gisaid.org/registration/terms-of-use/).

There will be dates that are missing sequence counts. These should be assumed to be 0.