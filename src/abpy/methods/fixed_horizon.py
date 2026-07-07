import numpy as np
from scipy import stats
from itertools import combinations
def _run_ttest(a, b, experiment):
    alpha = experiment.alpha
    alternative = getattr(experiment, 'alternative', 'two-tail')

    if 'corrected_var' in a and 'corrected_var' in b:
        var_a, var_b = a['corrected_var'], b['corrected_var']
    else:
        var_a, var_b = a['std'] ** 2, b['std'] ** 2

    if var_a == 0 and var_b == 0:
        raise ZeroDivisionError("your variance for each variant is zero")

    se = np.sqrt(var_a / a['n'] + var_b / b['n'])
    df = (var_a / a["n"] + var_b / b["n"]) ** 2 / (
        (var_a / a["n"]) ** 2 / (a["n"] - 1) +
        (var_b / b["n"]) ** 2 / (b["n"] - 1)
    )
    diff = b["mean"] - a["mean"]
    stat = diff / se

    if alternative == 'one-tail':
        # H1: variant b > variant a
        p_value = 1 - stats.t.cdf(stat, df=df)
        t_crit = stats.t.ppf(1 - alpha, df=df)
        ci = [diff - t_crit * se, None]  # one-sided lower bound, unbounded above
    else:
        p_value = 2 * (1 - stats.t.cdf(abs(stat), df=df))
        margin = stats.t.ppf(1 - alpha / 2, df=df) * se
        ci = [diff - margin, diff + margin]

    pooled_std = np.sqrt((var_a + var_b) / 2)
    effect_size = diff / pooled_std
    z = diff / se

    from abpy.power.power_analysis import compute_power, compute_required_sample_size
    power = compute_power(z, alpha, experiment)
    required_n = compute_required_sample_size(var_a + var_b, diff, alpha, experiment.power, alternative)

    return {
        "test_used": "ttest",
        "stat": stat,
        "p_value": p_value,
        "confidence_interval": ci,
        "effect_size": effect_size,
        "power": power,
        "required_sample_size": required_n,
        "alternative": alternative
    }


def _run_ztest(a, b, experiment):
    alpha = experiment.alpha
    alternative = getattr(experiment, 'alternative', 'two-tail')

    p_pool = (a["successes"] + b["successes"]) / (a["n"] + b["n"])
    se = np.sqrt(p_pool * (1 - p_pool) * (1 / a["n"] + 1 / b["n"]))
    diff = b["p"] - a["p"]
    stat = diff / se

    if alternative == 'one-tail':
        p_value = 1 - stats.norm.cdf(stat)
        z_crit = stats.norm.ppf(1 - alpha)
        ci = [diff - z_crit * se, None]
    else:
        p_value = 2 * (1 - stats.norm.cdf(abs(stat)))
        margin = stats.norm.ppf(1 - alpha / 2) * se
        ci = [diff - margin, diff + margin]

    effect_size = diff / np.sqrt(p_pool * (1 - p_pool))
    z = diff / se

    from abpy.power.power_analysis import compute_power, compute_required_sample_size
    power = compute_power(z, alpha, experiment)
    variance_sum = a["p"] * (1 - a["p"]) + b["p"] * (1 - b["p"])
    required_n = compute_required_sample_size(variance_sum, diff, alpha, experiment.power, alternative)

    return {
        "test_used": "ztest",
        "stat": stat,
        "p_value": p_value,
        "confidence_interval": ci,
        "effect_size": effect_size,
        "power": power,
        "required_sample_size": required_n,
        "alternative": alternative
    }


def _run_mannwhitney(a, b, experiment):
    alpha = experiment.alpha
    alternative = getattr(experiment, 'alternative', 'two-tail')
    # scipy convention: alternative='less' tests that x (a) is stochastically
    # smaller than y (b), i.e. b > a -- matching this library's one-tail convention.
    scipy_alt = 'less' if alternative == 'one-tail' else 'two-sided'

    stat, p_value = stats.mannwhitneyu(a['raw_values'], b['raw_values'], alternative=scipy_alt)
    effect_size = stat / (a["n"] * b["n"])
    z = (stat - (a['n'] * b['n']) / 2) / np.sqrt(a['n'] * b['n'] * (a['n'] + b['n'] - 1) / 12)

    from abpy.power.power_analysis import compute_power, compute_required_sample_size
    power = compute_power(z, alpha, experiment)

    # normal-theory approximation of the required sample size for Mann-Whitney;
    # exact non-parametric sample size formulas require distributional
    # assumptions that aren't available here, so this is reported as an
    # approximation based on the observed variance of the raw values.
    var_a = a['raw_values'].var()
    var_b = b['raw_values'].var()
    diff = b['raw_values'].mean() - a['raw_values'].mean()
    required_n = compute_required_sample_size(var_a + var_b, diff, alpha, experiment.power, alternative)

    return {
        "test_used": "mannwhitney",
        "stat": -z,
        "p_value": p_value,
        "confidence_interval": None,
        "effect_size": effect_size,
        "power": power,
        "required_sample_size": required_n,
        "alternative": alternative
    }


def _run_poisson(a, b, experiment):
    alpha = experiment.alpha
    alternative = getattr(experiment, 'alternative', 'two-tail')

    se = np.sqrt(a["rate"] / a["n"] + b["rate"] / b["n"])
    diff = b["rate"] - a["rate"]
    stat = diff / se

    if alternative == 'one-tail':
        p_value = 1 - stats.norm.cdf(stat)
        z_crit = stats.norm.ppf(1 - alpha)
        ci = [diff - z_crit * se, None]
    else:
        p_value = 2 * (1 - stats.norm.cdf(abs(stat)))
        margin = stats.norm.ppf(1 - alpha / 2) * se
        ci = [diff - margin, diff + margin]

    effect_size = diff / np.sqrt((a["rate"] + b["rate"]) / 2)
    z = diff / se

    from abpy.power.power_analysis import compute_power, compute_required_sample_size
    power = compute_power(z, alpha, experiment)
    required_n = compute_required_sample_size(a["rate"] + b["rate"], diff, alpha, experiment.power, alternative)

    return {
        "test_used": "poisson",
        "stat": stat,
        "p_value": p_value,
        "confidence_interval": ci,
        "effect_size": effect_size,
        "power": power,
        "required_sample_size": required_n,
        "alternative": alternative
    }


def _run_pair(a, b, experiment):
    test_map = {
        "z_test":       _run_ztest,
        "t_test":       _run_ttest,
        "mann_whitney": _run_mannwhitney,
        "poisson":      _run_poisson
    }
    test_fn = test_map.get(experiment.test_type)
    if test_fn is None:
        raise ValueError(
            f"Unknown test type '{experiment.test_type}'. "
            f"Choose from: {list(test_map.keys())}"
        )
    return test_fn(a, b, experiment)


def apply_fixed_horizon(extracted, experiment):
    variants = list(extracted.keys())
    results = {}
    for var_a, var_b in combinations(variants, 2):
        key = f"{var_a} vs {var_b}"
        results[key] = _run_pair(
            extracted[var_a],
            extracted[var_b],
            experiment
        )
    return results
