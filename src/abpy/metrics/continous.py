from scipy import stats
'''
order would be like this:
variant || metric1 || metric2 || pre_metric1 || pre_metric2 
'''
def _extract_from_direct(data):
    result={}
    for var,col in data.items():
        result[var]=col
    return result
def _csv_to_dict(data,col_name):
    result={}
    for var,col in data.groupby('variant'):
        result[var]={}
        result[var]['n']=len(col)
        result[var]['mean']=col[col_name].mean()
        result[var]['std']=col[col_name].std()
        result[var]['raw_values']=data[data['variant']==var][col_name]
    return result
def _check_normality(data,col_name,experiment):
    sum=0
    total=0
    for var,col in data.groupby('variant'):
        total+=1
        n=len(col)
        if n>30:
            sum+=1
        else:
            _,p_value=stats.shapiro(col[col_name])
            if p_value>experiment.alpha:
                sum+=1
    return sum==total
def _extract_from_csv(data,experiment):
    results={}
    is_normal=False
    if experiment.cuped:
        if experiment.is_ratio:
            from abpy.preprocessing.delta import apply_delta
            results=apply_delta(data,'metric1_cuped','metric2_cuped')
            is_normal=True
        else:
            is_normal=_check_normality(data,'metric_cuped',experiment)
            results=_csv_to_dict(data,'metric_cuped')
    else:
        if experiment.is_ratio:
            from abpy.preprocessing.delta import apply_delta
            results=apply_delta(data,'metric1','metric2')
            is_normal=True
        else:
            is_normal=_check_normality(data,'metric',experiment)
            results=_csv_to_dict(data,'metric')
    return (results,is_normal)

def run_continous_test(experiment):
    is_normal=False
    if experiment.input=='direct':
        results1=_extract_from_direct(experiment.data)
        is_normal=True
    else:
        results=(experiment.data)

    if experiment.cuped:
        if experiment.input=='csv':
         from abpy.preprocessing.cuped import apply_cuped
         if experiment.is_ratio:
             results,experiment.theta1=apply_cuped(results,'metric1','pre_metric1')
             results,experiment.theta2=apply_cuped(results,'metric2','pre_metric2')
         else:
             results,experiment.theta=apply_cuped(results,'metric','pre_metric')
        else:
         raise ValueError(f"It will not work with the predefined sample estimates")
    if experiment.input=='csv':
        results1,is_normal=_extract_from_csv(results,experiment)
    if is_normal:
        experiment.test_type='t_test'
    else:
        experiment.test_type='mann_whitney'
    from abpy.methods.fixed_horizon import apply_fixed_horizon
    results_final=apply_fixed_horizon(results1,experiment)
    if experiment.correction:
        from abpy.corrections.multiple import apply_corrections
        results_final=apply_corrections(results_final,experiment.correction)
    return results_final






    