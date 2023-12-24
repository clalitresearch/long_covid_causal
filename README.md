# Longitudinal Effects of SARS-CoV-2 Infection on Laboratory Biomarkers

## SMD, SMDiD Analysis

The following code uses the analysis and visualization functions defined above to recreate the analysis and visualizations. Each analysis is described based on the "Methods" section of the paper and further explained on the paper.

To properly run the code, it is required to load a laboratory dataframe (*LAB_ALL*) with all 63 laboratory results throughout the study period. Two laboratory dataframes are required for complete analysis, one for for the general cohort and one for lab-specific cohorts.

The dataframe should also include the following columns:  
rowid(int, identifier), adjusted_lab_month/time(int, adjusted month/day, infection day/month = 0), group(str, infected/control group), result(float)

#### *Due to data privacy regulations, data is unavailable.*

### Import libraries


```python
# import open-source libraries, with relevant versions
import pandas as pd # v1.5.3
import numpy as np # v1.23.5
import seaborn as sns # v0.12.2
import matplotlib.pyplot as plt # v3.7.0

# Functions created for analysis and visualization:
import analsys_funcs as anl
import visualization_funcs as vis
```


```python
# LAB_ALL = pd.DataFrame which should include all laboratry results for the cohort throughout the study period. 
LAB_ALL = pd.read_parquet("path/to/labs/dataframe")
FIGS_PATH = "path/to/project/figures"
laboratories_list = LAB_ALL.lab_name.unique() # LAB_ALL assume laboratory biomarker names ordered for visualization
```

### Standardized mean difference (SMD) analysis

The following equation determines the calculation of SMD for month i, where μ stands for the mean laboratory biomarker across each study group during month i, and σ stand for the standard deviation of that biomarker prior to the infection:

$SMD_{month_i}=\frac{{μ_{infected-month_i}-μ_{control-month_i}}} {\sqrt{\frac{{{σ_{infected-month_i}^2+σ_{control-month_i}^2}}}{2}}}$

Results of SMD analysis are shown in Supplementary figures 1,2.


```python
# SMD Analysis
heatmap_data_SMD, heatmap_data_SMD_CI = [], [] # list for SMD analysis
for lab_name in laboratories_list: # calculate each lab separately
    lab_df = LAB_ALL[LAB_ALL.lab_name == lab_name]
    monthly_SMD, monthly_SMD_CI = anl.single_laboratory_analysis_SMD(lab_df, lab_name)
    heatmap_data_SMD.append(monthly_SMD); heatmap_data_SMD_CI.append(monthly_SMD_CI)
# Visualization preparation
smd_df_cols = ['Lab name',*range(-12,13)]
df_SMD = pd.DataFrame(heatmap_data_SMD, smd_df_cols)
df_SMD_rounded = anl.create_SMD_df_for_heatmap(df_SMD, pd.DataFrame(heatmap_data_SMD_CI, smd_df_cols))
# Visualization
vis.plot_SMD_heatmap(df_SMD, df_SMD_rounded, title: "SMD Plot", FIGS_PATH, save_name="SMDiD Plot.png")
```

### Standardized mean difference in differences (SMDiD) analyis
For this calculation, we first calculated the SMD between the study groups for the 12-month period prior to the index date as follows:  $SMD_{pre-infection}=\frac{{μ_{infected-pre}-μ_{control-pre}}} {\sqrt{\frac{{{σ_{infected-pre}^2+σ_{control-pre}^2}}}{2}}}$

The standardized mean difference-in-differences (SMDiD) for month i during the follow-up period was then calculated for each cohort as follows: 
$SMDiD_{month_i}=SMD_{month_i}-SMD_{pre-infection}$

For both SMD and SMDiD analysis, 95% confidence intervals for each estimate were calculated using 500 bootstraps. To account for months with small sample sizes, if the interval intersects with zero, no effect is assumed. Otherwise, we consider a positive effect if the $SMDiD_{month_i}$ is greater than 0, and a negative effect if the $SMDiD_{month_i}$ is smaller than 0.

Results of SMDiD analysis are shown in Figure 2, and Supplementary Figures 3,4 (for stratified analysis)


```python
# SMDiD Analysis
heatmap_data_SMDiD, heatmap_data_SMDiD_CI = [], [] # list for SMDiD
for lab_name in laboratories_list: # calculate each lab separately
    lab_df = LAB_ALL[LAB_ALL.lab_name == lab_name]
    # SMDiD Analysis
    monthly_SMDiD, monthly_SMDiD_CI = anl.single_laboratory_analysis_SMDiD(lab_df, lab_name)
    heatmap_data_SMDiD.append(monthly_SMDiD); heatmap_data_SMDiD_CI.append(monthly_SMDiD_CI)
# Visualization preparation
smdid_df_cols = ['Lab name',*range(0,13)]
df_SMDiD = pd.DataFrame(heatmap_data_SMDiD, smdid_df_cols)
df_SMDID_rounded = anl.create_SMDiD_df_for_heatmap(df_SMDiD, pd.DataFrame(heatmap_data_SMDiD_CI, smdid_df_cols))
# Visualization
vis.plot_SMDiD_heatmap(df_SMDiD, df_SMDID_rounded, title: "SMDiD Plot", FIGS_PATH, save_name="SMDiD Plot.png")
```

### Monthly average results with DiD correction

The infection group is plotted both with the original average results, as well as a calculated DiD line, shifted to account for differences in baseline between the study groups prior to infection. The DiD in each month i is calculated as follows: $Infection_{DiD-correction_i}=μ_{i_{infected}}-diff_{pre}$ 

where $ diff_{pre}=μ_{inf-pre}-μ_{control-pre}$ is the mean difference between the infected and control groups in the 12 months prior to infection

Results of chosen laboratories with this analysis are shown in Figure 3.


```python
for lab_name in laboratories_list: # calculate each lab separately
    lab_df = LAB_ALL[LAB_ALL.lab_name == lab_name]
    # Line plot visualization - mean values with DiD correction
    vis.visualize_DiD_line_plot(lab_df, lab_dir, lab_name)
```
