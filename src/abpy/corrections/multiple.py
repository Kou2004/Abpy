import numpy as np
def _bonferroni(results):
    n_tests = len(results)
    adjusted = {}
    for key, result in results.items():
        adjusted_p = min(result["p_value"] * n_tests, 1.0)
        adjusted[key] = {
            **result,
            "p_value_adjusted": adjusted_p,
            "correction":       "bonferroni"
        }
    return adjusted
def _holm(results):
    n_tests = len(results)
    sorted_items = sorted(results.items(), key=lambda x: x[1]["p_value"])
    adjusted = {}
    for i, (key, result) in enumerate(sorted_items):
        adjusted_p = min(result["p_value"] * (n_tests - i), 1.0)
        adjusted[key] = {
            **result,
            "p_value_adjusted": adjusted_p,
            "correction":       "holm"
        } 
    return adjusted
def _benjamini_hochberg(results):
    n_tests = len(results)
    sorted_items = sorted(results.items(), key=lambda x: x[1]["p_value"])
    adjusted = {}
    for i, (key, result) in enumerate(sorted_items):
        rank      = i + 1
        adjusted_p = min(result["p_value"] * n_tests / rank, 1.0)
        adjusted[key] = {
            **result,
            "p_value_adjusted": adjusted_p,
            "correction":       "bh"
        } 
    return adjusted
def apply_corrections(results, correction):
    if len(results) < 2:
        raise ValueError(
            "Multiple comparison correction requires at least 2 comparisons. "
            "Use correction=None for a single A/B test."
        )
    correction_map = {
        "bonf": _bonferroni,
        "holm":       _holm,
        "bh":         _benjamini_hochberg
    }
    correction_fn = correction_map.get(correction)
    if correction_fn is None:
        raise ValueError(
            f"Unknown correction '{correction}'. "
            f"Choose from: {list(correction_map.keys())}"
        )
    return correction_fn(results)