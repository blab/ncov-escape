# Escape scores

Escape scores are calculated by [Nextclade](https://github.com/nextstrain/nextclade) using approach from [Escape calculator for SARS-CoV-2 RBD](https://jbloomlab.github.io/SARS2_RBD_Ab_escape_maps/escape-calc/).

1. Install `nextclade` CLI following instructions at [docs.nextstrain.org](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli.html)

2. Provision `sars-cov-2-21L` dataset with:
```
nextclade dataset get --name 'sars-cov-2-21L' --output-dir 'data/sars-cov-2-21L'
```

3. Download canonical Pango sequences provisioned by [@corneliusroemer](https://github.com/corneliusroemer) via:
```
curl -fsSL https://github.com/corneliusroemer/pango-sequences/blob/main/data/pango_consensus_sequences.fasta.zstd?raw=true -o data/pango_consensus_sequences.fasta.zstd
```
And decompress with:
```
zstdcat data/pango_consensus_sequences.fasta.zstd > data/pango_consensus_sequences.fasta
```

4. Run Nextclade on each Pango lineage with:
```
nextclade run \
   --input-dataset data/sars-cov-2-21L \
   --output-all=output/ \
   data/pango_consensus_sequences.fasta
```

This will generate the file `output/nextclade.tsv` containing columns `immune_escape` and `ace2_binding` and rows for each Pango lineage.

5. Prune this file to relevant columns:
```
tsv-select -H -f seqName,clade,Nextclade_pango,partiallyAliased,immune_escape,ace2_binding output/nextclade.tsv > escape_scores.tsv
```
