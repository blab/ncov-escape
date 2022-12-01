# Predictors

This script calculates escape score and ACE-2 binding for RBD using [Nextclade](https://github.com/nextstrain/nextclade) based on approach from [Escape calculator for SARS-CoV-2 RBD](https://jbloomlab.github.io/SARS2_RBD_Ab_escape_maps/escape-calc/). These calculations have been implemented in Nextclade and are available across Pango lineages at https://nextstrain.org/nextclade/sars-cov-2/21L?c=immune_escape.

Additionally, this script uses Nextclade tree to calculate RBD mutation count relative to BA.2, non-RBD S1 mutation count relative to BA.2 and non-S1 mutation count relative to BA.2. S1 rather than full spike was chosen based on the results from [Kistler et al](https://bedford.io/papers/kistler-sarscov2-adaptive-evolution/).

1. Download JSON from the Nextclade build:
```
curl -fsSL https://data.nextstrain.org/nextclade_sars-cov-2_21L.json -o nextclade_sars-cov-2_21L.json.gz
```

2. Decompress this file:
```
gzip -d nextclade_sars-cov-2_21L.json.gz -c > nextclade_sars-cov-2_21L.json
```

3. Extract relevant tip attributes to TSV:
```
python extract_tip_attributes.py --json nextclade_sars-cov-2_21L.json > predictors.tsv
```

This TSV looks like
```
seqName	clade	Nextclade_pango	partiallyAliased	immune_escape	ace2_binding	rbd_count	non_rbd_spike_count	non_s1_count
BQ.1	22E (Omicron)	BQ.1	BA.5.3.1.1.1.1.1	0.6306944	0.66401	5	6	15
BQ.1.1	22E (Omicron)	BQ.1.1	BA.5.3.1.1.1.1.1.1	0.83137808	0.77543	6	6	16
```
