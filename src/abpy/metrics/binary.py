def _extract_from_direct(data):
    result={}
    for var,col in data.items():
        proportions=col['successes']/col['n']
        result[var]=col
        result[var]['p']=proportions
    return result
def _extract_from_csv(data):
    result={}
    for var,col in data.groupby('variant'):
        result[var]={}
        result[var]['n']=len(col)
        result[var]['successes']=col['metric'].sum()
        result[var]['p']=col['metric'].sum()/len(col)
    return result
def run_binary_test(experiment):
    experiment.test_type='z_test'
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





