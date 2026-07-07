import os
import pandas as pd
import numpy as np
'''
Input : 9 parameters
Methods:
fit(data)
results()
plot()
'''
class Experiment:
    valid_types=['ab','multivariate']
    valid_metrics=['continous','discrete','binary']
    valid_inputs=['csv','direct']
    valid_corrections=['bonf','bh','holm',None]
    valid_alternatives=['one-tail','two-tail']

    def __init__(self,type='ab',metric='continous',input=None,alpha=0.05,power=0.8,cuped=False,correction=None,is_ratio=False,alternative='two-tail'):
        self.input=input
        self.metric=metric
        self.alpha=alpha
        self.cuped=cuped
        self.type=type
        self.correction=correction
        self.is_ratio=is_ratio
        self.power=power
        self.alternative=alternative
        self.data=None # its a guard check
        self._validate_inputs()
    def _validate_inputs(self):
        if self.type.lower() not in self.valid_types:
            raise ValueError(f"Invalid Input Type:'{self.type}' Choose from:'{self.valid_types}'")
        if self.metric.lower() not in self.valid_metrics:
            raise ValueError(f"Invalid Input Type:'{self.metric}' Choose from:'{self.valid_metrics}'")
        if self.input is None or not isinstance(self.input,str):
            raise ValueError(f"'input' is required. Choose from:'{self.valid_inputs}'")
        if self.input.lower() not in self.valid_inputs:
            raise ValueError(f"Invalid Input Type:'{self.input}' Choose from:'{self.valid_inputs}'")
        if self.correction is not None and self.correction.lower() not in self.valid_corrections:
            raise ValueError(f"Invalid Input Type:'{self.correction}' Choose from:'{self.valid_corrections}'")
        if self.cuped not in [True,False]:
            raise ValueError(f"Invalid Input Type:'{self.cuped}' Choose boolean type")
        if self.is_ratio not in [True,False]:
            raise ValueError(f"Invalid Input Type:'{self.is_ratio}' Choose boolean type")
        if self.cuped and self.metric.lower()!='continous':
            raise ValueError(f"'cuped' is only supported for metric='continous' (it silently had no effect for '{self.metric}' before this check)")
        if self.is_ratio and self.metric.lower()!='continous':
            raise ValueError(f"'is_ratio' is only supported for metric='continous' (it silently had no effect for '{self.metric}' before this check)")
        if not isinstance(self.alternative,str) or self.alternative.lower() not in self.valid_alternatives:
            raise ValueError(f"Invalid Input Type:'{self.alternative}' Choose from:'{self.valid_alternatives}'")
        numeric_inputs=[self.power,self.alpha]
        for i in numeric_inputs:
            if i<=0 or i>=1:
                raise ValueError(f"Invalid input type :'{i}'. alpha and power must be strictly between 0 and 1")
        self.type=self.type.lower()
        self.metric=self.metric.lower()
        self.input=self.input.lower()
        self.alternative=self.alternative.lower()
    def fit(self,data):
        if self.input=='csv':
            self.load_csv(data)
        else:
            self.load_direct(data)
    def load_direct(self,data):
        if not isinstance(data,dict):
            raise ValueError(f"direct input must be a dictionary")
        if len(data)<2:
            raise ValueError(f"Minimum 2 variants are required")
        keys={
            'continous':['n','mean','std'],
            'binary':['n','successes'],
            'discrete':['n','events']
        }
        current_item=keys[self.metric]
        for var,stats in data.items():
            for key in current_item:
                if key not in stats:
                    raise KeyError(f"'{key}' Data is missing in variant'{var}'")
        for var,stats in data.items():
            for key,val in stats.items():
                if not isinstance(val,(int,float)):
                    raise ValueError(f"data type is wrong in '{val}' from variant'{var}'")
                if val<0 and key!='mean':
                    raise ValueError(f"data  is negative in '{val}' from variant'{var}'")
        self.data=data

    def load_csv(self,data):
        if not os.path.exists(data):
            raise FileNotFoundError(f"file not found")
        if not data.endswith('.csv'):
            raise ValueError(f"file must be csv")
        self.data=pd.read_csv(data)
        if self.data.isnull().any().any():
            raise ValueError(f"it contains missing value")
        requried_cols=['variant']# minimum check for columns
        if self.is_ratio:
            requried_cols+=['metric1','metric2']
        else:
            requried_cols+=['metric']
        for col in requried_cols:
            if col not in self.data.columns:
                raise KeyError(f"'{col}'column is missing")
        measurement_cols=[c for c in ['metric1','metric2'] if c in self.data.columns] if self.is_ratio else [c for c in ['metric'] if c in self.data.columns]
        for col in measurement_cols:
            if (self.data[col]<0).any() and self.metric!='continous':
                raise ValueError(f"data is negative in column '{col}'")
    def results(self):
        if self.data is None:
            raise RuntimeError(f"No data loaded")
        if self.metric=='binary':
            from abpy.metrics.binary import run_binary_test
            result= run_binary_test(self)# dict of dict
        if self.metric=='continous':
            from abpy.metrics.continous import run_continous_test
            result= run_continous_test(self)
        if self.metric=='discrete':
            from abpy.metrics.discrete import run_discrete_test
            result= run_discrete_test(self)
        self._raw_results=result
        from abpy.results.output import ExperimentResult
        result_final=ExperimentResult(result,self)
        result_final.summary()
        return result_final
    def plot(self):
        if self.data is None:
            raise RuntimeError(f"No data loaded")
        from abpy.visualizations.plots import plot_results
        return plot_results(self)
        

       
