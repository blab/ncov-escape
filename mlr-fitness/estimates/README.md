# MLR fitness estimates using evofr

1. Clone [evofr](https://github.com/blab/evofr) repo elsewhere in the file system and install via
```
poetry build
pip install <path-to-wheel>
```

2. Run the notebook `evofr_mlr.ipynb` via `jupyter notebook`. This will export a file called `growth_advantages.tsv` that contains a series of posterior samples (1 to 100) with growth advantage per Pango lineage.
