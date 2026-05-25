# The Bias Fingerprint

[![ci](https://github.com/mahmood726-cyber/bias-forensics/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/mahmood726-cyber/bias-forensics/actions/workflows/ci.yml) [![codeql](https://github.com/mahmood726-cyber/bias-forensics/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/mahmood726-cyber/bias-forensics/actions/workflows/codeql.yml) [![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) [![python: 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)

When eight publication-bias methods are applied to the same meta-analysis, how often do they reach the same conclusion? We applied four detection tests and four correction methods to 307 Cochrane reviews with at least five studies from the Pairwise70 dataset. Each review received a bias fingerprint summarising which tests flagged bias and how corrections shifted the pooled estimate, classified as Clean, Suspected, Confirmed, or Discordant. Only 54 reviews (17.6%) were Clean while 42 were Discordant, giving a discordance rate of 13.7% and a median risk-ratio shift below 0.05 (95% CI 0.02 to 0.08). Agreement ranged from 98% between PET-PEESE and limit meta-analysis to 73% between trim-and-fill and regression corrections, with mean relative effect shift of 1.01. In roughly one in seven Cochrane reviews, the apparent direction of the pooled conclusion depends on which bias method is chosen. These methods detect statistical asymmetry rather than bias itself and cannot separate publication bias from clinical or methodological heterogeneity.

**Live dashboard:** <https://mahmood726-cyber.github.io/biasforensics/>

## Run

Open `index.html` (or `index.html`) in any modern browser. No build step.

For local development:

```bash
python -m http.server 8000
# then open http://localhost:8000/
```

## Test

```bash
python -m pytest -q
```

The suite under `tests/` includes 2 test file(s).

## Repo layout

| Path | Purpose |
|---|---|
| `index.html` | the dashboard (main artifact) |
| `index.html` | landing page |
| `tests/` | pytest tests |
| `e156-submission/` | E156 micro-paper bundle |
| `E156-PROTOCOL.md` | project metadata (E156 entry #12) |

## License

See `LICENSE` (MIT).
