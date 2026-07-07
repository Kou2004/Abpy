def apply_cuped(data,target_col,pre_col):
    '''
    only valid for csv files
    structure of the csv:
    user_id || variant || metric || pre_metric
    '''
    expected_x=data[pre_col].mean()
    var_y=data[pre_col].var()
    if var_y==0:
        raise ValueError(f"your exp data must have variability")
    theta=data[target_col].cov(data[pre_col])/var_y
    data[f'{target_col}_cuped']=data[target_col]-theta*(data[pre_col]-expected_x)
    return (data,theta)

    

    