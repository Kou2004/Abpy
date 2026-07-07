import numpy as np
from scipy import stats


def compute_power(z, alpha, experiment):
    '''
    Achieved power for the observed effect.
    Respects experiment.alternative:
      - two-tail: uses alpha/2 critical value and |z| (direction agnostic)
      - one-tail: uses alpha critical value and signed z (effect must be
        in the hypothesized direction, otherwise power correctly drops
        toward 0 instead of being inflated by taking abs())
    '''
    alternative = getattr(experiment, 'alternative', 'two-tail')
    if alternative == 'one-tail':
        z_alpha = stats.norm.ppf(1 - alpha)
        power = stats.norm.cdf(z - z_alpha)
    else:
        z_alpha = stats.norm.ppf(1 - alpha / 2)
        power = stats.norm.cdf(abs(z) - z_alpha)
    return power


def compute_required_sample_size(variance_sum, diff, alpha, power, alternative='two-tail'):
    '''
    Required sample size PER VARIANT to detect `diff` given the summed
    variance of the two groups, at the requested alpha/power.

    n = (z_alpha + z_beta)^2 * variance_sum / diff^2
    '''
    if diff == 0:
        return float('inf')

    if alternative == 'one-tail':
        z_alpha = stats.norm.ppf(1 - alpha)
    else:
        z_alpha = stats.norm.ppf(1 - alpha / 2)

    z_beta = stats.norm.ppf(power)
    n = ((z_alpha + z_beta) ** 2) * variance_sum / (diff ** 2)
    return int(np.ceil(n))