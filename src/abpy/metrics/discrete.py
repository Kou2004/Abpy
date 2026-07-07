def _extract_from_direct(data):
    result={}
    for var,col in data.items():
        result[var]=col
        result[var]['rate']=col['events']/col['n']
    return result
def _extract_from_csv(data):
    result={}
    for var,col in data.groupby('variant'):
        result[var]={}
        result[var]['n']=len(col)
        result[var]['rate']=col['metric'].sum()/len(col)
    return result
def run_discrete_test(experiment):
    experiment.test_type='poisson'
    if experiment.input=='direct':
        results=_extract_from_direct(experiment.data)
    else:
        results=_extract_from_csv(experiment.data)
   
    from abpy.methods.fixed_horizon import apply_fixed_horizon
    results_final=apply_fixed_horizon(results,experiment)
    if experiment.correction:
        from abpy.corrections.multiple import apply_corrections
        results_final=apply_corrections(results_final,experiment.correction)
    return results_final





