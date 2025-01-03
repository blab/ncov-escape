# About

Forked from v12 of the canonical ncov workflow [github.com/nextstrain/ncov](https://github.com/nextstrain/ncov) on Aug 12, 2022.

Visible live at [nextstrain.org/groups/blab/ncov/escape](https://nextstrain.org/groups/blab/ncov/escape).

# Running

Run with:
```
snakemake -c 1 -j 1 -p --profile escape_profiles/escape auspice/ncov_escape.json
```
if Nextstrain toolkit is installed locally. Or if using Docker and CLI, run with:
```
nextstrain build . --profile escape_profiles/escape auspice/ncov_escape.json
```

Visualize locally with:
```
auspice view --datasetDir auspice/
```

Push live with:
```
nextstrain remote upload nextstrain.org/groups/blab auspice/ncov_escape*.json
```
