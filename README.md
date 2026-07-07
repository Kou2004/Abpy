
# abpy — Frequentist A/B Testing Library

A Python library for rigorous frequentist A/B and multivariate testing, with built-in support for CUPED variance reduction, Delta method for ratio metrics, and multiple comparison corrections.

---

## Installation

```bash
pip install abpy
```

---

## Quick Start

```python
from abpy.router import Experiment

# binary metric — direct input
exp = Experiment(type='ab', metric='binary', input='direct', alpha=0.05)
exp.fit({
    'control':   {'n': 1000, 'successes': 100},
    'treatment': {'n': 1000, 'successes': 130}
})
exp.results()
exp.plot()
```

---

## Features

- **Three metric types** — binary, continuous, count
- **Two input modes** — CSV file or direct summary statistics
- **Automatic test selection** — z-test, t-test, Mann-Whitney U, Poisson rate test
- **CUPED** — variance reduction using pre-experiment data
- **Delta method** — correct variance estimation for ratio metrics
- **Multiple comparison corrections** — Bonferroni, Holm, Benjamini-Hochberg
- **Multivariate testing** — all pairwise comparisons with corrections
- **Power analysis** — achieved power for every test
- **Dashboard visualisation** — 6-panel results plot

---

## API Reference

### `Experiment(type, metric, input, alpha, cuped, correction, is_ratio,power,alternative)`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `type` | str | `'ab'` | `'ab'` or `'multivariate'` |
| `metric` | str | `'continous'` | `'binary'`, `'continous'`, `'discrete'` |
| `input` | str | required | `'csv'` or `'direct'` |
| `alpha` | float | `0.05` | significance level |
| `cuped` | bool | `False` | variance reduction (CSV only) |
| `correction` | str | `None` | `'bonf'`, `'holm'`, `'bh'` |
| `is_ratio` | bool | `False` | ratio metric — applies Delta method (CSV only) |
| `power` | float | `0.8` | power |
| `alternative` | str | `two-tail` | type of alternative hypothesis |

---

## Usage Examples

### Binary Metric — CSV Input

```python
import pandas as pd
from abpy.router import Experiment

exp = Experiment(type='ab', metric='binary', input='csv', alpha=0.05)
exp.fit('experiment.csv')   # requires columns: variant, metric
exp.results()
exp.plot()
```

CSV structure:
```
variant,metric
control,0
control,1
treatment,1
treatment,0
...
```

### Continuous Metric — Direct Input

```python
exp = Experiment(type='ab', metric='continous', input='direct', alpha=0.05)
exp.fit({
    'control':   {'n': 1000, 'mean': 45.2, 'std': 12.3},
    'treatment': {'n': 1000, 'mean': 48.7, 'std': 11.9}
})
exp.results()
```

### Discrete / Count Metric

```python
exp = Experiment(type='ab', metric='discrete', input='direct', alpha=0.05)
exp.fit({
    'control':   {'n': 1000, 'events': 300},
    'treatment': {'n': 1000, 'events': 375}
})
exp.results()
```

### Multivariate with Holm Correction

```python
exp = Experiment(
    type='multivariate',
    metric='binary',
    input='direct',
    alpha=0.05,
    correction='holm'
)
exp.fit({
    'control':    {'n': 2000, 'successes': 200},
    'treatment1': {'n': 2000, 'successes': 240},
    'treatment2': {'n': 2000, 'successes': 260}
})
exp.results()
```

### CUPED — Variance Reduction

Reduces variance using pre-experiment data, requiring fewer users to reach significance.

```python
exp = Experiment(type='ab', metric='continous', input='csv', alpha=0.05, cuped=True)
exp.fit('experiment_with_pre.csv')
exp.results()
```

CSV structure for CUPED:
```
variant,metric,pre_metric
control,45.2,44.1
treatment,48.7,43.9
...
```

### Ratio Metric — Delta Method

For metrics like CTR = clicks / impressions where both vary per user.

```python
exp = Experiment(type='ab', metric='continous', input='csv', alpha=0.05, is_ratio=True)
exp.fit('ratio_experiment.csv')
exp.results()
```

CSV structure for ratio metrics:
```
variant,metric1,metric2
control,3,10
treatment,5,10
...
```
CSV structure for ratio metrics (if Cuped applied)
```
variant,metric1,metric2,pre_metric1,pre_metric2
control,3,10,2,6
treatment,5,10,3,5
...
```
---

## Test Selection Logic

```
Binary metric    → z-test (pooled proportion)
Continuous
  n >= 30        → t-test (Welch's approximation)
  n < 30         → Shapiro-Wilk normality check
    normal       → t-test
    not normal   → Mann-Whitney U
  is_ratio=True  → Delta method → t-test
Discrete/Count   → Poisson rate test
```

---

## Output

Every test returns:

```
Test used            which statistical test was selected
Alternative          type 
Test statistic       z, t, or U value
P-value              raw p-value
Adjusted p-value     if correction applied
Confidence interval  95% CI on the lift
Effect size          Cohen's h (binary) or d (continuous)
Power                achieved power given sample size
Required sample size The sample size for the given power to be achieved
Recommendation       "2nd variant wins" / "no significant difference" / "underpowered"
```

---

## Corrections Reference

| Method | Controls | Best For |
|---|---|---|
| `bonf` | FWER | Few comparisons, false positives very costly |
| `holm` | FWER | More powerful than Bonferroni, same control |
| `bh` | FDR | Many comparisons, some false positives acceptable |

---

## Visualisation

```python
exp.plot()   # renders 6-panel dashboard
```

Dashboard panels:
1. P-value vs significance threshold
2. Confidence interval on lift
3. Effect size (Cohen's benchmark)
4. Power curve vs sample size
5. Metric distribution by variant
6. Cuped applied metric distribution by variant

---

## Dependencies

```
numpy >= 1.21.0
scipy >= 1.7.0
pandas >= 1.3.0
matplotlib >= 3.4.0
```

---

## License

MIT License — see LICENSE file for details.
