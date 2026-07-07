import numpy as np
def apply_delta(data,col1,col2):
    '''
    only for ratio metrics
    metric= x/y
    csv (minimal strcture):
    variant || x || y || metric  
    '''
    required_columns = ["variant", col1, col2]
    for col in required_columns:
        if col not in data.columns:
            raise ValueError(f"Missing required column for Delta method: '{col}'")

    if (data[col2] == 0).any():
        raise ValueError(
            f"Denominator column '{col2}' contains zeros. "
            "Cannot compute ratio with zero denominator."
        )

    results = {}

    for variant, group in data.groupby("variant"):

        Y = group[col1]
        X = group[col2]
        n = len(group)

        mu_Y = Y.mean()
        mu_X = X.mean()

        if mu_X == 0:
            raise ValueError(
                f"Mean of denominator is zero for variant '{variant}'. "
                "Cannot apply Delta method."
            )

        var_Y  = Y.var()
        var_X  = X.var()
        cov_YX = Y.cov(X)
        mean   = mu_Y / mu_X

        corrected_var = (1 / n) * (
            var_Y / mu_X**2
            - 2 * (mu_Y / mu_X**3) * cov_YX
            + (mu_Y**2 / mu_X**4) * var_X
        )

        std = np.sqrt(corrected_var * n)

        results[variant] = {
            "n":             n,
            "mean":          mean,
            "std":           std,
            "corrected_var": corrected_var
        }

    return results