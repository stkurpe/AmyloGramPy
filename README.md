# AmyloGram-Py

Fast, testable Python reimplementation of the original R **AmyloGram** prediction
path.

AmyloGram-Py is designed to reproduce the predictions of the reference R
package while making high-throughput FASTA prediction substantially faster and
easier to run inside automated pipelines.

Original AmyloGram repository: [michbur/AmyloGram](https://github.com/michbur/AmyloGram)

## What AmyloGram-Py Reimplements

The Python implementation follows the same prediction logic as the R reference
path:

```text
protein sequence
  -> overlapping 6-mers
  -> AmyloGram amino-acid group encoding
  -> binary multigram features
  -> exported ranger random forest
  -> max 6-mer amyloid probability per protein
```

The R implementation remains the reference oracle. AmyloGram-Py is a compatible
Python prediction path built from exported reference fixtures and validated
against R outputs.

## What Was Changed

Compared with the original R execution path, AmyloGram-Py changes the execution
strategy, not the biological prediction target:

- Reimplemented AmyloGram amino-acid cleaning and group encoding in Python.
- Reimplemented 6-mer feature generation to match R/biogram fixtures.
- Exported the trained `ranger` random forest prediction path from R fixtures.
- Precomputed all degenerate 6-mer probabilities: `6^6 = 46656` possible encoded
  windows.
- Added a compact binary lookup-table format for production prediction.
- Replaced per-window forest traversal during FASTA prediction with direct array
  lookup.
- Added streaming FASTA prediction with rolling 6-mer encoding to avoid storing
  all windows for each protein in memory.
- Added Markdown and JSON run reports with input and lookup-table SHA256 hashes.
- Added skipped-record TSV output with machine-readable skip reasons.
- Added top-hit TSV output without requiring all predictions to remain in memory.

## Why The Outputs Are Considered Equivalent

AmyloGram-Py was compared against the original R AmyloGram path on three FASTA
samples. The validation checked both biological classification agreement and
numeric probability agreement.

### Count-Level Agreement

R AmyloGram and AmyloGram-Py produced the same number of valid predictions, the
same number of skipped records, and the same number of amyloid/non-amyloid
predictions.

| Sample | Input records | R predicted | Py predicted | R skipped | Py skipped | R amyloid | Py amyloid | R non-amyloid | Py non-amyloid | Discordant labels |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SRR32060215 | 27980 | 27958 | 27958 | 22 | 22 | 26772 | 26772 | 1186 | 1186 | 0 |
| SRR32060233 | 25196 | 25178 | 25178 | 18 | 18 | 24359 | 24359 | 819 | 819 | 0 |
| SRR32060234 | 28327 | 28312 | 28312 | 15 | 15 | 27248 | 27248 | 1064 | 1064 | 0 |

### Prediction-Level Agreement

The shared sequence IDs had zero discordant binary labels. Probability
agreement was also extremely close: mean absolute probability differences were
about `5.9e-7`, and the maximum observed absolute difference was `5.06e-6`.
These small differences are consistent with numeric/export precision differences
between the R path and the Python lookup-table path.

| Sample | Shared unique IDs | R-only unique IDs | Py-only unique IDs | Discordant labels | Mean abs prob delta | Median abs prob delta | Max abs prob delta | Near-threshold records |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SRR32060215 | 27956 | 0 | 0 | 0 | 0.00000059 | 0.00000048 | 0.00000506 | 189 |
| SRR32060233 | 25178 | 0 | 0 | 0 | 0.00000059 | 0.00000048 | 0.00000411 | 155 |
| SRR32060234 | 28312 | 0 | 0 | 0 | 0.00000059 | 0.00000047 | 0.00000411 | 166 |

### Skipped-Record Agreement

Both implementations failed or skipped the same biological/input edge cases:
records that were empty after amino-acid cleaning or shorter than six amino
acids. AmyloGram-Py reports these records explicitly in a skipped TSV.

| Sample | R invalid/short | Py skipped | R sequence errors | Py skipped reasons |
|---|---:|---:|---:|---|
| SRR32060215 | 22 | 22 | 0 | shorter_than_6: 22 |
| SRR32060233 | 18 | 18 | 0 | shorter_than_6: 18 |
| SRR32060234 | 15 | 15 | 0 | shorter_than_6: 15 |

Together, these checks support the interpretation that AmyloGram-Py reproduces
the R AmyloGram prediction path for the tested inputs while providing a faster
execution backend.

## Runtime Performance

AmyloGram-Py was substantially faster than the R execution path in the tested
samples.

| Sample | R runtime | Py runtime | Speedup | Py records/sec |
|---|---:|---:|---:|---:|
| SRR32060215 | 67.62 min | 3.79 sec | 1069.5x | 7370.1 |
| SRR32060233 | 51.11 min | 4.02 sec | 763.1x | 6265.7 |
| SRR32060234 | 51.77 min | 3.37 sec | 922.9x | 8411.6 |

The speedup comes from using a precomputed 6-mer probability table and a
streaming rolling-code implementation instead of repeatedly traversing the
random forest for every 6-mer window during prediction.

## Installation

Install the latest public code directly from GitHub:

```bash
python3 -m pip install "amylogram-py @ git+https://github.com/stkurpe/AmyloGramPy.git"
```

The package requires Python `>=3.10`.


## Predict A FASTA

```bash
PYTHONPATH=src python3 -m amylogram_py.cli \
  input.fasta \
  amylogram_predictions.csv \
  --sixmer-table tests/fixtures/amylogram_sixmer_probabilities.bin \
  --report-md amylogram_prediction_report.md \
  --report-json amylogram_prediction_report.json \
  --skipped-tsv amylogram_skipped.tsv \
  --top-k 100 \
  --top-tsv amylogram_top_hits.tsv
```

Output columns:

- `Sequence_ID`
- `AmyloGram_Prob`
- `AmyloGram_Pred`

The Markdown and JSON reports record run-level metrics such as total FASTA
records, predicted records, skipped records, elapsed seconds, throughput,
threshold, and SHA256 hashes for the input FASTA and lookup table. The skipped
TSV lists records that were too short or empty after amino-acid cleaning. The
top-hit TSV keeps the highest-probability proteins without storing all
predictions in memory.


