import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats


def _alt_label(experiment):
    return "one-tailed" if getattr(experiment, "alternative", "two-tail") == "one-tail" else "two-tailed"


def plot_pvalue(result, experiment, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    p_value = result["p_value_adjusted"] if result.get("p_value_adjusted") is not None else result["p_value"]
    alpha   = experiment.alpha

    colors = ["#2ecc71" if p_value < alpha else "#e74c3c"]

    ax.barh(["p-value"], [p_value], color=colors)
    ax.axvline(x=alpha, color="black", linestyle="--", label=f"α = {alpha}")
    ax.set_xlim(0, 1)
    ax.set_xlabel("p-value")
    ax.set_title(f"P-value vs Significance Threshold ({_alt_label(experiment)})")
    ax.legend()

    return ax


def plot_ci(result, experiment, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    ci = result["confidence_interval"]

    if ci is None:
        ax.text(0.5, 0.5, "CI not available\n(Mann-Whitney test)",
                ha="center", va="center", transform=ax.transAxes)
        return ax

    lower, upper = ci

    if upper is None:
        # one-tail test: the upper bound is unbounded (+inf). Draw the
        # known lower bound as a point with an arrow extending right,
        # rather than crashing on a None upper bound.
        span = max(abs(lower) * 0.75, 1.0)

        ax.plot([lower], [0], "o", color="#3498db", markersize=8)
        ax.annotate(
            "", xy=(lower + span, 0), xytext=(lower, 0),
            arrowprops=dict(arrowstyle="-|>", color="#3498db", lw=2)
        )
        ax.axvline(x=0, color="black", linestyle="--", label="no effect")
        ax.set_xlim(min(lower, 0) - span * 0.3, lower + span * 1.3)
        ax.set_xlabel("Lift")
        ax.set_title(f"Effect Size — One-Sided Lower Bound [{lower:.4f}, ∞)")
        ax.set_yticks([])
        ax.legend()
        return ax

    center = (lower + upper) / 2

    ax.errorbar(
        x          = [center],
        y          = [0],
        xerr       = [[center - lower], [upper - center]],
        fmt        = "o",
        color      = "#3498db",
        capsize    = 10,
        markersize = 8
    )

    ax.axvline(x=0, color="black", linestyle="--", label="no effect")
    ax.set_xlabel("Lift")
    ax.set_title(f"Effect Size with Confidence Interval ({_alt_label(experiment)})")
    ax.set_yticks([])
    ax.legend()

    return ax


def plot_power_curve(result, experiment, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    # generate sample sizes range
    n_current = experiment.data[list(experiment.data.keys())[0]]["n"] \
                if experiment.input == "direct" \
                else len(experiment.data) // 2

    n_range   = np.linspace(10, n_current * 2, 100)
    effect    = abs(result["effect_size"])
    alpha     = experiment.alpha

    from abpy.power.power_analysis import compute_power
    power_values = []
    for n in n_range:
        z = effect * np.sqrt(n / 2)
        power_values.append(compute_power(z, alpha, experiment))

    ax.plot(n_range, power_values, color="#3498db", linewidth=2)
    ax.axhline(y=experiment.power,         color="green",  linestyle="--", label=f"target power ({experiment.power})")
    ax.axhline(y=result["power"], color="red", linestyle="--", label=f"achieved power ({result['power']:.2f})")
    ax.axvline(x=n_current,   color="orange", linestyle="--", label=f"current n ({n_current})")
    ax.set_xlabel("Sample Size per Variant")
    ax.set_ylabel("Power")
    ax.set_title(f"Power Curve ({_alt_label(experiment)})")
    ax.set_ylim(0, 1)
    ax.legend()

    return ax


def plot_effect_size(result, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    effect = abs(result["effect_size"])

    # cohen's benchmark thresholds
    thresholds = {"small": 0.2, "medium": 0.5, "large": 0.8}

    if effect < thresholds["small"]:
        label = "small"
        color = "#e74c3c"
    elif effect < thresholds["medium"]:
        label = "small-medium"
        color = "#e67e22"
    elif effect < thresholds["large"]:
        label = "medium"
        color = "#f1c40f"
    else:
        label = "large"
        color = "#2ecc71"

    ax.barh(["effect size"], [effect], color=color)
    ax.axvline(x=0.2, color="gray", linestyle="--", alpha=0.5, label="small (0.2)")
    ax.axvline(x=0.5, color="gray", linestyle="-.",  alpha=0.5, label="medium (0.5)")
    ax.axvline(x=0.8, color="gray", linestyle=":",   alpha=0.5, label="large (0.8)")
    ax.set_xlabel("Effect Size")
    ax.set_title(f"Effect Size — {label}")
    ax.legend()

    return ax


def _get_distribution_columns(experiment):
    '''
    Work out which column(s) hold per-row metric values suitable for a
    histogram/bar distribution plot, for both the original metric and,
    if CUPED was applied, the CUPED-adjusted metric.

    Returns (original_col, cuped_col_or_None).
    Adds a synthetic 'metric_ratio'/'metric_ratio_cuped' column when the
    experiment is a ratio metric, since ratio experiments don't have a
    single 'metric' column (they use metric1/metric2 instead).
    '''
    data = experiment.data

    if experiment.metric == "continous" and experiment.is_ratio:
        if "metric_ratio" not in data.columns:
            data["metric_ratio"] = data["metric1"] / data["metric2"]
        original_col = "metric_ratio"

        cuped_col = None
        if experiment.cuped and "metric1_cuped" in data.columns and "metric2_cuped" in data.columns:
            if "metric_ratio_cuped" not in data.columns:
                data["metric_ratio_cuped"] = data["metric1_cuped"] / data["metric2_cuped"]
            cuped_col = "metric_ratio_cuped"

        return original_col, cuped_col

    original_col = "metric"
    cuped_col = "metric_cuped" if (experiment.cuped and "metric_cuped" in data.columns) else None
    return original_col, cuped_col


def plot_distributions(experiment, ax=None, title_suffix=""):
    if experiment.input == "direct":
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "Distribution plot requires CSV input",
                ha="center", va="center", transform=ax.transAxes)
        return ax

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    data = experiment.data

    if experiment.metric == "continous":
        col_name, _ = _get_distribution_columns(experiment)

    for variant, group in data.groupby("variant"):
        if experiment.metric == "binary":
            rate = group["metric"].mean()
            ax.bar(variant, rate, alpha=0.7, label=variant)
            ax.set_ylabel("Conversion Rate")

        elif experiment.metric == "continous":
            ax.hist(group[col_name], bins=30, alpha=0.6, label=variant)
            ax.set_ylabel("Frequency")

        elif experiment.metric == "discrete":
            rate = group["metric"].mean()
            ax.bar(variant, rate, alpha=0.7, label=variant)
            ax.set_ylabel("Average Count")

    ax.set_title(f"Metric Distribution by Variant{title_suffix}")
    ax.legend()

    return ax


def plot_cuped_comparison(experiment, ax=None):
   
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    if experiment.input == "direct":
        ax.text(0.5, 0.5, "CUPED comparison requires CSV input",
                ha="center", va="center", transform=ax.transAxes)
        return ax

    if experiment.metric != "continous":
        ax.text(0.5, 0.5, "CUPED is only supported for\ncontinuous metrics",
                ha="center", va="center", transform=ax.transAxes)
        return ax

    original_col, cuped_col = _get_distribution_columns(experiment)

    if not experiment.cuped or cuped_col is None:
        ax.text(0.5, 0.5, "CUPED not applied\n(no adjusted column found)",
                ha="center", va="center", transform=ax.transAxes)
        return ax

    data = experiment.data
    variants = list(data.groupby("variant").groups.keys())
    palette = plt.cm.tab10(np.linspace(0, 1, max(len(variants), 3)))

    for color, (variant, group) in zip(palette, data.groupby("variant")):
        ax.hist(group[cuped_col], bins=30, alpha=0.9, color=color,
                label=f"{variant} — CUPED", histtype="step", linewidth=2)

    ax.set_ylabel("Frequency")
    ax.set_xlabel("Metric value")
    ax.set_title(" CUPED-Adjusted Distribution")
    ax.legend(fontsize=8)

    return ax


def plot_dashboard(result, experiment):
    show_cuped = (
        experiment.cuped
        and experiment.metric == "continous"
        and experiment.input == "csv"
    )

    if show_cuped:
        fig = plt.figure(figsize=(16, 10))
        gs  = gridspec.GridSpec(2, 3, figure=fig)

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[0, 2])
        ax4 = fig.add_subplot(gs[1, 0])
        ax5 = fig.add_subplot(gs[1, 1])
        ax6 = fig.add_subplot(gs[1, 2])

        plot_pvalue(result, experiment, ax=ax1)
        plot_ci(result, experiment,     ax=ax2)
        plot_effect_size(result,        ax=ax3)
        plot_power_curve(result, experiment, ax=ax4)
        plot_distributions(experiment,  ax=ax5, title_suffix=" (Original)")
        plot_cuped_comparison(experiment, ax=ax6)
    else:
        fig = plt.figure(figsize=(14, 10))
        gs  = gridspec.GridSpec(2, 3, figure=fig)

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[0, 2])
        ax4 = fig.add_subplot(gs[1, 0])
        ax5 = fig.add_subplot(gs[1, 1:])

        plot_pvalue(result, experiment, ax=ax1)
        plot_ci(result, experiment,     ax=ax2)
        plot_effect_size(result,        ax=ax3)
        plot_power_curve(result, experiment, ax=ax4)
        plot_distributions(experiment,  ax=ax5)

    fig.suptitle("A/B Test Results Dashboard", fontsize=16, fontweight="bold")
    plt.tight_layout()
    plt.show()
    plt.close(fig)
    return fig


def plot_results(experiment):
    if not hasattr(experiment, "_raw_results"):
        raise RuntimeError("No results available. Call results() before plot()")
    result = experiment._raw_results
    if experiment.type == "ab":
        key    = list(result.keys())[0]
        return plot_dashboard(result[key], experiment)

    else:
        # multivariate — one dashboard per comparison
        figs = {}
        for key, comparison_result in result.items():
            figs[key] = plot_dashboard(comparison_result, experiment)
        return figs
