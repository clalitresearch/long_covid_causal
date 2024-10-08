# Trends of Common Laboratory Biomarkers after SARS-CoV-2 Infection

## SMD, SMDiD Analysis

The following code is an example of calculating each analysis mentioned in the manuscript. Each analysis is described based on the "Methods" section of the paper and further explained in the paper.

#### *Due to data privacy regulations, data is unavailable.*

### Import libraries


```python
# import open-source libraries, with relevant versions
import pandas as pd # v1.5.3
import numpy as np # v1.23.5
```

### Standardized mean difference (SMD) analysis

The following equation determines the calculation of SMD for month i, where μ stands for the mean laboratory biomarker across each study group during month i, and σ stand for the standard deviation of that biomarker prior to the infection:

$SMD_{month_i}=\frac{{μ_{infected-month_i}-μ_{control-month_i}}} {\sqrt{\frac{{{σ_{infected-month_i}^2+σ_{control-month_i}^2}}}{2}}}$

Results of SMD analysis are shown in Supplementary figures 1,2.

Example of calculation in Python:
```python
monthly_res = [] # contains the monthly SMD for each month
monthly_CI = [] # contains the monthly SMD CI for each month 
# Get mean result of each patient in each month (see "Methods")
mean_all_df = ... # table of the mean result value per patient per month of each group
for month in range(-12,13): # Calculate SMD and CI each month separately
    month_0 = mean_all_df[mean_all_df.month==month]
    lab_month_inf = month_0[month_0.group =='Infected']
    lab_month_con = month_0[month_0.group =='Control']
    month_std_inf = lab_month_inf["result"].std()
    month_std_con = lab_month_con["result"].std()
    month_mean_inf = lab_month_inf["result"].mean()
    month_mean_con = lab_month_con["result"].mean()
    # Calculate SMD
    smd = (month_mean_inf - month_mean_con) / np.sqrt((month_std_inf**2 + month_std_con**2) / 2) 
    # Number of bootstrap samples
    n_bootstrap = 500
    np.random.seed(1) # set random seed to ensure reproductability
    bootstrap_smds = []
    # Generate bootstrap samples
    for ii in range(n_bootstrap):
        smp_inf = lab_month_inf["result"].sample(frac=1, replace=True)
        smp_con = lab_month_con["result"].sample(frac=1, replace=True)
        # Calculate the standardized mean difference for each bootstrap sample
        smd_boot = (smp_inf.mean() - smp_con.mean()) / (np.sqrt((smp_inf.std()**2+smp_con.std()**2)/2))
        bootstrap_smds.append(smd_boot)
    # Calculate the confidence interval
    lower_ci = np.percentile(bootstrap_smds, 2.5)
    upper_ci = np.percentile(bootstrap_smds, 97.5)
    monthly_res.append(smd)
    monthly_CI.append((lower_ci, upper_ci))
```

### Standardized mean difference in differences (SMDiD) analyis
For this calculation, we first calculated the SMD between the study groups for the 12-month period prior to the index date as follows:  $SMD_{pre-infection}=\frac{{μ_{infected-pre}-μ_{control-pre}}} {\sqrt{\frac{{{σ_{infected-pre}^2+σ_{control-pre}^2}}}{2}}}$

The standardized mean difference-in-differences (SMDiD) for month i during the follow-up period was then calculated for each cohort as follows: 
$SMDiD_{month_i}=SMD_{month_i}-SMD_{pre-infection}$

For both SMD and SMDiD analysis, 95% confidence intervals for each estimate were calculated using 500 bootstraps. To account for months with small sample sizes, if the interval intersects with zero, no effect is assumed. Otherwise, we consider a positive effect if the $SMDiD_{month_i}$ is greater than 0, and a negative effect if the $SMDiD_{month_i}$ is smaller than 0.

Results of SMDiD analysis are shown in Figure 2, and Supplementary Figures 3,4 (for stratified analysis)

Example of calculation in Python:
```python
monthly_res = []
monthly_CI = []
# Get mean result of each patient in each month 
mean_all_df = ... # table of the mean result value per patient per month of each group
# Keep laboratories prior to infection to calculate SMD diff
lab_temp_df_pre  = lab_df[(-365 < lab_df.date) & (lab_df.date < 0)]
std_inf = lab_temp_df_pre[lab_temp_df_pre.group=='Infected']["result"].std()
std_con = lab_temp_df_pre[lab_temp_df_pre.group=='Control']["result"].std()
# Get mean results pre-infection
mean_pre_df = mean_pre_df[mean_pre_df.month < 0]
# Calculate standardized mean diff pre infection 
smd_diff = (mean_pre_df[mean_pre_df.group=='Infected']["result"].mean() - 
            mean_pre_df[mean_pre_df.group=='Control']["result"].mean())/np.sqrt((std_inf**2+std_con**2)/2)
for month in range(0,13): # Calculate SMDiD and CI each month separately
    month_0 = mean_all_df[mean_all_df.month==month]
    lab_month_inf = month_0[month_0.group =='Infected']
    lab_month_con = month_0[month_0.group =='Control']
    month_std_inf = lab_month_inf["result"].std()
    month_std_con = lab_month_con["result"].std()
    month_mean_inf = lab_month_inf["result"].mean()
    month_mean_con = lab_month_con["result"].mean()
    # Calculate standardized mean difference in differences
    smdid = (month_mean_inf - month_mean_con) / np.sqrt((month_std_inf**2 + month_std_con**2) / 2) - smd_diff
    # Number of bootstrap samples
    n_bootstrap = 500
    np.random.seed(1) # set random seed to ensure reproductability
    bootstrap_smds = []
    # Generate bootstrap samples
    for ii in range(n_bootstrap):
        smp_inf = lab_month_inf["result"].sample(frac=1, replace=True)
        smp_con = lab_month_con["result"].sample(frac=1, replace=True)
        # Calculate the standardized mean difference in difference for each bootstrap sample
        smd_boot = (smp_inf.mean() - smp_con.mean()) / (np.sqrt((smp_inf.std()**2+smp_con.std()**2)/2)) - smd_diff
        bootstrap_smds.append(smd_boot)
    # Calculate the confidence interval
    lower_ci = np.percentile(bootstrap_smds, 2.5)
    upper_ci = np.percentile(bootstrap_smds, 97.5)
    monthly_res.append(smdid) 
    monthly_CI.append((lower_ci, upper_ci))
```

### Monthly average results with DiD correction

The infection group is plotted both with the original average results, as well as a calculated DiD line, shifted to account for differences in baseline between the study groups prior to infection. The DiD in each month i is calculated as follows: $Infection_{DiD-correction_i}=μ_{i_{infected}}-diff_{pre}$ 

where $diff_{pre}=μ_{inf-pre}-μ_{control-pre}$ is the mean difference between the infected and control groups in the 12 months prior to infection

The results of the chosen laboratories with this analysis are shown in Figure 3.

Example of calculation in Python:
```python
mean_values_per_group = ... # table that contains monthly mean result per group 
treated_avg = mean_values_per_group[mean_values_per_group.group=='Infected']['result']
untreated_avg = mean_values_per_group[mean_values_per_group.group=='Control']['result']
mean_values_per_group['difference_in_differences'] = treated_avg.reset_index(drop=True)-untreated_avg.reset_index(drop=True)
```
