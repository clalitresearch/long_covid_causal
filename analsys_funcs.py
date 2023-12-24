import pandas as pd # v1.5.3
import numpy as np # v1.23.5


def single_laboratory_analysis_SMDiD(lab_df: pd.DataFrame, lab_name: str):
    """Calculate SMDiD between infected and control for months 0-12 for a single laboratory
    See "Methods" section for further information about formulas.

    Args:
        lab_df (pd.DataFrame): Single laboratory results for matched cohort, both infected and control, after censoring
            columns: rowid(int, identifier), adjusted_lab_month/time(int, adjusted month/day, infection day/month = 0),
            group(str, infected/control group), result(int)
        lab_name (str): _description_

    Returns:
            monthly_res(list): Ordered monthly SMDiD, vector of ints, 0=infection month
            monthly_CI(list): Ordered vector of tuples that represents CIs
    """
    monthly_res = [lab_name]
    monthly_CI = [lab_name]
    # Get mean result of each patient in each month 
    mean_all_df = lab_df.groupby(['group','adjusted_lab_month','rowid'])['result'].mean().reset_index()
    # Keep laboratories prior to infection to calculate SMD diff
    lab_temp_df_pre  = lab_df[(-365 < lab_df.adjusted_lab_time) & (lab_df.adjusted_lab_time < 0)]
    std_inf = lab_temp_df_pre[lab_temp_df_pre.group=='Infected']["result"].std()
    std_con = lab_temp_df_pre[lab_temp_df_pre.group=='Control']["result"].std()
    # Get mean results pre-infection
    mean_pre_df = mean_pre_df[mean_pre_df.adjusted_lab_month < 0]
    # Calculate standardized mean diff pre infection 
    smd_diff = (mean_pre_df[mean_pre_df.group=='Infected'].result.mean() - 
                mean_pre_df[mean_pre_df.group=='Control'].result.mean())/np.sqrt((std_inf**2+std_con**2)/2)
    for m in range(0,13): # Calculate SMDiD and CI each month separately
        month_0 = mean_all_df[mean_all_df.adjusted_lab_month==m]
        lab_month_inf = month_0[month_0.group =='Infected']
        lab_month_con = month_0[month_0.group =='Control']
        month_std_inf = lab_month_inf.result.std()
        month_std_con = lab_month_con.result.std()
        month_mean_inf = lab_month_inf.result.mean()
        month_mean_con = lab_month_con.result.mean()
        # Calculate standardized mean difference in differences
        smdid = (month_mean_inf - month_mean_con) / np.sqrt((month_std_inf**2 + month_std_con**2) / 2) - smd_diff
        # Number of bootstrap samples
        n_bootstrap = 500
        np.random.seed(1) # set random seed to ensure reproductability
        bootstrap_smds = []
        # Generate bootstrap samples
        for ii in range(n_bootstrap):
            smp_inf = lab_month_inf.result.sample(frac=1, replace=True)
            smp_con = lab_month_con.result.sample(frac=1, replace=True)
            # Calculate the standardized mean difference in difference for each bootstrap sample
            smd_boot = (smp_inf.mean() - smp_con.mean()) / (np.sqrt((smp_inf.std()**2+smp_con.std()**2)/2)) - smd_diff
            bootstrap_smds.append(smd_boot)
        # Calculate the confidence interval
        lower_ci = np.percentile(bootstrap_smds, 2.5)
        upper_ci = np.percentile(bootstrap_smds, 97.5)
        monthly_res.append(smdid) 
        monthly_CI.append((lower_ci, upper_ci))
    return monthly_res, monthly_CI


def single_laboratory_analysis_SMD(lab_df: pd.DataFrame, lab_name: str):
    """
    Calculate the SMD between infected and control for 12 months pre-inf -> 12 months post for a single laboratory
    See "Methods" section for further information about formulas.
    
    Args:
        lab_df (pd.DataFrame): Single laboratory results for matched cohort, both infected and control, after censoring
            columns: rowid(int, identifier), adjusted_lab_month/time(int, adjusted month/day, infection day/month = 0),
            group(str, infected/control group), result(int)
        lab_name (str): _description_

    Returns:
            monthly_res(list): Ordered monthly SMDiD, vector of ints, 0=infection month
            monthly_CI(list): Ordered vector of tuples that represents CIs
    """
    monthly_res = [lab_name]
    monthly_CI = [lab_name]
    # Get mean result of each patient in each month (see "Methods")
    mean_all_df = lab_df.groupby(['group','adjusted_lab_month','rowid'])['result'].mean().reset_index()
    for m in range(-12,13): # Calculate SMD and CI each month separately
        month_0 = mean_all_df[mean_all_df.adjusted_lab_month==m]
        lab_month_inf = month_0[month_0.group =='Infected']
        lab_month_con = month_0[month_0.group =='Control']
        month_std_inf = lab_month_inf.result.std()
        month_std_con = lab_month_con.result.std()
        month_mean_inf = lab_month_inf.result.mean()
        month_mean_con = lab_month_con.result.mean()
        # Calculate SMD
        smd = (month_mean_inf - month_mean_con) / np.sqrt((month_std_inf**2 + month_std_con**2) / 2) 
        # Number of bootstrap samples
        n_bootstrap = 500
        np.random.seed(1) # set random seed to ensure reproductability
        bootstrap_smds = []
        # Generate bootstrap samples
        for ii in range(n_bootstrap):
            smp_inf = lab_month_inf.result.sample(frac=1, replace=True)
            smp_con = lab_month_con.result.sample(frac=1, replace=True)
            # Calculate the standardized mean difference for each bootstrap sample
            smd_boot = (smp_inf.mean() - smp_con.mean()) / (np.sqrt((smp_inf.std()**2+smp_con.std()**2)/2))
            bootstrap_smds.append(smd_boot)
        # Calculate the confidence interval
        lower_ci = np.percentile(bootstrap_smds, 2.5)
        upper_ci = np.percentile(bootstrap_smds, 97.5)
        monthly_res.append(smd, 2)
        monthly_CI.append((lower_ci, upper_ci))
    return monthly_res, monthly_CI


def calculate_diff_in_diff_for_lineplot(mean_monthly_res: pd.DataFrame):
    """Calculate difference in difference for mean lab values (per patient) in each month

    Args:
        mean_monthly_res (pd.DataFrame): Dataframe with laboratory results of a single laboratory biomarker

    Returns:
        mean_values_per_group(pd.DataFrame): mean lab values (per patient) in each month with correction
    """
    # calculate the differences in average monthly results between infected and control, return
    mean_values_per_group = mean_monthly_res.groupby(['group','adjusted_lab_month'])['result'].mean().reset_index()
    treated_avg = mean_values_per_group[mean_values_per_group.group=='Infected']['result']
    untreated_avg = mean_values_per_group[mean_values_per_group.group=='Control']['result']
    mean_values_per_group['difference_in_differences'] = treated_avg.reset_index(drop=True)-untreated_avg.reset_index(drop=True)
    return mean_values_per_group


def create_SMDiD_df_for_heatmap(df_smdid: pd.DataFrame, df_CI: pd.DataFrame):
    """Create a dataframe of annotation for SMDiD visualization, round each diff value where the CI intersect with 0

    Args:
        df_smdid (pd.DataFrame): Columns: lab_name, monthly SMDiDs
        df_CI (pd.DataFrame): Columns: lab_name, monthly CIs of SMDiDs

    Returns:
        pd.DataFrame: For annotation, will present the effect each month unless CI intersects with 0
    """
    all_rows = []
    for row in df_smdid.iterrows(): # row should be in order as they will appear on the heatmap later
        lab_name = row[0]
        row_CI = df_CI.loc[lab_name] # get the same row for the CIs df
        row_smds = row[1] # get smds
        row_fin = [lab_name] # add lab name to row
        for ii in range(0,13):
            low_lim,high_lim = row_CI[ii][0], row_CI[ii][1]
            if pd.Interval(low_lim,high_lim).overlaps(pd.Interval(0,0)):
                row_fin.append(0)
            else:
                row_fin.append(row_smds[ii])
        all_rows.append(row_fin)
    cols = ['lab_name',*range(0,13)]
    return pd.DataFrame(all_rows, columns=cols).set_index(['lab_name'])


def create_SMD_df_for_heatmap(df_smd: pd.DataFrame, df_CI: pd.DataFrame):
    """Create a dataframe of annotation for SMD visualization, round each diff value where the CI intersect with 0
    
    Args:
        df_smdid (pd.DataFrame): Columns: lab_name, monthly SMDiDs
        df_CI (pd.DataFrame): Columns: lab_name, monthly CIs of SMDiDs

    Returns:
        pd.DataFrame: For annotation, will present the effect each month unless CI intersects with 0
    """
    all_rows = []
    for row in df_smd.iterrows(): # row should be in order as they will appear on the heatmap later
        lab_name = row[0]
        row_CI = df_CI.loc[lab_name] # get the same row for the CIs df
        row_smds = row[1] # get smds
        row_fin = [lab_name] # add lab name to row
        for ii in range(-12,13):
            low_lim,high_lim = row_CI[ii][0], row_CI[ii][1]
            if pd.Interval(low_lim,high_lim).overlaps(pd.Interval(0,0)):
                row_fin.append(0)
            else:
                row_fin.append(row_smds[ii])
        all_rows.append(row_fin)
    cols = ['lab_name',*range(-12,13)]
    return pd.DataFrame(all_rows, columns=cols).set_index(['lab_name'])
