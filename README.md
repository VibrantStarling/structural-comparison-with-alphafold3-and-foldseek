[![DOI](https://zenodo.org/badge/976134394.svg)](https://doi.org/10.5281/zenodo.17379945)

# structural-comparison-with-alphafold3-and-foldseek
various scripts and collections of code for processing alphafold, foldseek, metapredict, and interproscan data

This is largely a collection of code, NOT finished or polished standalone scripts.
You will need to edit file paths and uncomment functions to make the plots you desire.

Raw data can be found in `toxo-afum-fgram-merged-df-MODIFICATION-CHANGED-2-interproscan-FILTERED-3.csv`. This contains a whole lot of raw data and is not particularly user friendly. However, it acts as the source for scripts like `plotting-merged-data.py`,`plddt-counts-plots.py`,`pairplotting.py`, `foldseek-filtering.py`, and `scoring-model-scores-v2.py`

`plotting-merged-data.py`,`plddt-counts-plots.py`,`pairplotting.py`, `foldseek-filtering.py`, and `scoring-model-scores-v2.py` are the scripts that make all of the figures in the paper to be titled '_The promise of AlphaFold for gene structure annotation_'
