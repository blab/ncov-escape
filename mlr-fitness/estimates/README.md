# MLR fitness estimates using evofr

1. Install latest version of [evofr](https://github.com/blab/evofr) via pip:
```
pip install evofr
```

2. Run script `run-mlr-model.py` via:
```
python -u run-mlr-model.py --config mlr-config.yaml --export-path . --data-name "pango"
```

This produces the file `pango_results.json`.

3. This file includes unnecessary information of full frequencies trajectories. The script `prune-mlr-results.py` walks through this JSON and extracts just the useful `site` of `ga`:
```
python prune-mlr-results.py --results pango_results.json > growth_advantages.tsv
```

And exports the TSV file called `growth_advantages.tsv` that contains the median growth advantage and 80% credible interval per Pango lineage.
