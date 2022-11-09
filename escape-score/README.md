# Escape scores

Escape scores are calculated by [Nextclade](https://github.com/nextstrain/nextclade) using approach from [Escape calculator for SARS-CoV-2 RBD](https://jbloomlab.github.io/SARS2_RBD_Ab_escape_maps/escape-calc/). These calculations have been implemented in Nextclade and are available across Pango lineages at https://nextstrain.org/nextclade/sars-cov-2/21L?c=immune_escape.

1. Download JSON from the Nextclade build:
```
curl -fsSL https://data.nextstrain.org/nextclade_sars-cov-2_21L.json -o nextclade_sars-cov-2_21L.json.gz
```

2. Decompress this file:
```
gzip -d nextclade_sars-cov-2_21L.json.gz -c > nextclade_sars-cov-2_21L.json
```

3. Flatten this JSON:
```
python flatten_auspice_json.py --json nextclade_sars-cov-2_21L.json --output nextclade_sars-cov-2_21L_flat.json
```

4. Extract relevant tip attributes to TSV:
```
python extract_tip_attributes.py --json nextclade_sars-cov-2_21L.json > escape_scores.tsv
```

This TSV looks like
```
seqName	clade	Nextclade_pango	partiallyAliased	immune_escape	ace2_binding
BQ.1	22E (Omicron)	BQ.1	BA.5.3.1.1.1.1.1	0.630694401494532	0.66401
BQ.1.1	22E (Omicron)	BQ.1.1	BA.5.3.1.1.1.1.1.1	0.8313780783921618	0.77543
```
