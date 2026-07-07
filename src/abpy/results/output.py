import pandas as pd
'''
input :
experiment object
raw results => nested dictionary
output : depends on the test type (ab or multivariate)
'''
class ExperimentResult:

    def __init__(self, raw_results, experiment):
        self.raw_results = raw_results
        self.experiment  = experiment
        self._parse()

    def _parse(self):
        # for a/b — single comparison
        if self.experiment.type == "ab":
            key          = list(self.raw_results.keys())[0]
            result       = self.raw_results[key]

            self.comparison           = key
            self.test_used            = result["test_used"]
            self.stat                 = result["stat"]
            self.p_value              = result["p_value"]
            self.p_value_adjusted     = result.get("p_value_adjusted", None)
            self.confidence_interval  = result["confidence_interval"]
            self.effect_size          = result["effect_size"]
            self.power                = result["power"]
            self.required_sample_size = result.get("required_sample_size", None)
            self.alternative          = result.get("alternative", self.experiment.alternative)
            self.correction           = result.get("correction", None)
            self.recommendation       = self._recommend(result)

        # multivariate — multiple comparisons
        else:
            self.comparisons    = {}
            self.recommendations = {}

            for key, result in self.raw_results.items():
                self.comparisons[key]     = result
                self.recommendations[key] = self._recommend(result)

            self.winner = self._find_winner()

    def _recommend(self, result):
        alpha   = self.experiment.alpha
        power   = result["power"]
        p_value = result["p_value_adjusted"] if result.get("p_value_adjusted") is not None else  result["p_value"]
        stat    = result["stat"]

        if power < self.experiment.power:
            return "underpowered — inconclusive"
        elif p_value < alpha:
            if stat > 0:
                return "2nd variant in that key wins"
            else:
                return "1st variant in that key wins"
        else:
            return "no significant difference"

    def _find_winner(self):
        if self.experiment.type != "multivariate":
            return None

        # find comparison with smallest adjusted p-value;
        # ties broken by largest absolute test statistic
        best_key    = None
        best_pvalue = 1.0
        best_stat_abs = -1.0

        for key, result in self.comparisons.items():
            p = result["p_value_adjusted"] if result.get("p_value_adjusted") is not None else  result["p_value"]
            stat_abs = abs(result["stat"])
            if p < best_pvalue or (p == best_pvalue and stat_abs > best_stat_abs):
                best_pvalue   = p
                best_stat_abs = stat_abs
                best_key      = key

        if best_pvalue < self.experiment.alpha:
            if self.comparisons[best_key]['stat']>0:
                return f"later one among the '{best_key}' wins"
            else:
                return f"first one among the '{best_key}' wins"
        return "no winner — no significant differences found"

    @staticmethod
    def _format_ci(ci):
        if ci is None:
            return "n/a (Mann-Whitney test)"
        lower, upper = ci
        if upper is None:
            return f"[{lower:.4f}, \u221e)  (one-sided)"
        return f"[{lower:.4f}, {upper:.4f}]"

    def summary(self):
        if self.experiment.type == "ab":
            print(f"\n{'='*50}")
            print(f"  A/B Test Results — {self.comparison}")
            print(f"{'='*50}")
            print(f"  Test used:            {self.test_used}")
            print(f"  Alternative:          {self.alternative}")
            print(f"  Test statistic:       {self.stat:.4f}")
            print(f"  P-value:              {self.p_value:.4f}")
            if self.p_value_adjusted is not None:
                print(f"  Adjusted p-value:     {self.p_value_adjusted:.4f} ({self.correction})")
            print(f"  Confidence interval:  {self._format_ci(self.confidence_interval)}")
            print(f"  Effect size:          {self.effect_size:.4f}")
            print(f"  Power (achieved):     {self.power:.4f}")
            if self.required_sample_size is not None:
                print(f"  Required sample size: {self.required_sample_size} per variant "
                      f"(to detect the observed effect at power={self.experiment.power}, alpha={self.experiment.alpha})")
            print(f"  Recommendation:       {self.recommendation}")
            print(f"{'='*50}\n")

        else:
            print(f"\n{'='*50}")
            print(f"  Multivariate Test Results")
            print(f"{'='*50}")
            for key, result in self.comparisons.items():
                p = result["p_value_adjusted"] if result.get("p_value_adjusted") is not None else  result["p_value"]
                print(f"  {key}:")
                print(f"    alternative:    {result.get('alternative', self.experiment.alternative)}")
                print(f"    p-value:        {result['p_value']:.4f}")
                if self.experiment.correction:
                    print(f"    adjusted p-value: {result['p_value_adjusted']:.4f} ({result['correction']})")
                print(f"    confidence interval: {self._format_ci(result.get('confidence_interval'))}")
                print(f"    effect size:    {result['effect_size']:.4f}")
                print(f"    power (achieved): {result['power']:.4f}")
                if result.get("required_sample_size") is not None:
                    print(f"    required sample size: {result['required_sample_size']} per variant")
                print(f"    recommendation: {self.recommendations[key]}")
            print(f"\n  Winner: {self.winner}")
            print(f"{'='*50}\n")

    def to_dict(self):
        if self.experiment.type == "ab":
            return {
                "comparison":          self.comparison,
                "test_used":           self.test_used,
                "alternative":         self.alternative,
                "stat":                self.stat,
                "p_value":             self.p_value,
                "p_value_adjusted":    self.p_value_adjusted,
                "confidence_interval": self.confidence_interval,
                "effect_size":         self.effect_size,
                "power":               self.power,
                "required_sample_size": self.required_sample_size,
                "correction":          self.correction,
                "recommendation":      self.recommendation
            }
        else:
            return {
                "comparisons":    self.comparisons,
                "recommendations": self.recommendations,
                "winner":         self.winner
            }

    def to_dataframe(self):
        if self.experiment.type == "ab":
            return pd.DataFrame([self.to_dict()])
        else:
            rows = []
            for key, result in self.comparisons.items():
                rows.append({
                    "comparison":     key,
                    "alternative":    result.get("alternative", self.experiment.alternative),
                    "p_value":       result["p_value_adjusted"] if result.get("p_value_adjusted") is not None else  result["p_value"],
                    "effect_size":    result["effect_size"],
                    "power":          result["power"],
                    "required_sample_size": result.get("required_sample_size"),
                    "recommendation": self.recommendations[key]
                })
            return pd.DataFrame(rows)

