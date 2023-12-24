import seaborn as sns # v0.12.2
import matplotlib.pyplot as plt # v3.7.0
import pandas as pd # v1.5.3
import analsys_funcs as anl

# Assumed orderd of laboratories:
HEMATO_RANGE = range(0,21)
COAG_RANGE = range(21,25)
METAB_RANGE = range(25,50)
INFL_RANGE = range(50,52)
CHOLES_RANGE = range(52,57)
VITM_RANGE = range(57,60)
ENDO_RANGE = range(60,63)

def plot_SMDiD_heatmap(df: pd.DataFrame, df_rounded: pd.DataFrame, title: str, 
                       figs_folder: str, save_name=False,v_min=-0.5, v_max=0.5, figsize=(13,18)):
    """Plots the main visualization (heatmap) - SMDiD post infection
    main parameters include the dataframes that includes SMDiD and SMDiD after rounding based on CI
    Labs are considered to be in order for the visualization

    Args:
        df (pd.DataFrame): pd.Dataframe of effects for each month for each laboratory biomarker
        df_rounded (pd.DataFrame): Rounded pd.Dataframe of effects for annotation
        title (str):  Main title of the plot
        figs_folder (str): Path to the figures folder of the relevant experiment
        save_name (bool, optional): Relevant save name, will save in "figs_folder" path. Defaults when no saving.
        v_min (float, optional): Set the color bar range. Defaults to -0.5.
        v_max (float, optional): Set the color bar range. Defaults to 0.5.
        figsize (tuple, optional): Plot size. Defaults to (13,18).
    """
    # split into cbc/coagulation/metabolic/inflammatory/lipid/vitamins/thyroid
    hemato_df = df.iloc[HEMATO_RANGE]; hemato_rounded = df_rounded.iloc[HEMATO_RANGE]
    coagulation_df = df.iloc[COAG_RANGE]; coagulation_rounded = df_rounded.iloc[COAG_RANGE]
    metab_df = df.iloc[METAB_RANGE]; metab_rounded = df_rounded.iloc[METAB_RANGE]
    infl_df = df.iloc[INFL_RANGE]; infl_rounded = df_rounded.iloc[INFL_RANGE]
    choles_df = df.iloc[CHOLES_RANGE]; choles_rounded = df_rounded.iloc[CHOLES_RANGE]
    vit_df = df.iloc[VITM_RANGE]; vit_rounded = df_rounded.iloc[VITM_RANGE]
    endo_seru_df = df.iloc[ENDO_RANGE]; endo_rounded = df_rounded.iloc[ENDO_RANGE]
    relative_heights = [] # Ensure all sub-plots will contain rows of the same height
    for df_t in [hemato_df, coagulation_df, metab_df, infl_df, choles_df, vit_df, endo_seru_df]:
        relative_heights.append(df_t.shape[0])
    # create subplots
    fig, axes = plt.subplots(nrows=7,  figsize=figsize,
                             gridspec_kw={'height_ratios':relative_heights}, constrained_layout=True)
    # plot each heatmap separately, annotate based on df_rounded
    g1 = sns.heatmap(hemato_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[0], cmap="vlag", annot=hemato_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                     xticklabels=True, yticklabels=True, cbar=False)
    g2 = sns.heatmap(coagulation_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[1], cmap="vlag", annot=coagulation_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                     xticklabels=True, yticklabels=True, cbar=False)
    g3 = sns.heatmap(metab_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[2], cmap="vlag", annot=metab_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                        xticklabels=True, yticklabels=True, cbar_kws={"shrink": 1})  
    g4 = sns.heatmap(infl_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[3], cmap="vlag", annot=infl_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                     xticklabels=True, yticklabels=True, cbar=False)  
    g5 = sns.heatmap(choles_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[4], cmap="vlag", annot=choles_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                     xticklabels=True, yticklabels=True, cbar=False)    
    g6 = sns.heatmap(vit_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[5], cmap="vlag", annot=vit_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                     xticklabels=True, yticklabels=True, cbar=False)    
    g7 = sns.heatmap(endo_seru_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[6], cmap="vlag", annot=endo_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                     xticklabels=True, yticklabels=True, cbar=False)         
    # set title, remove y label, and ensure only text above the visualization threshold is presented in each heatmap
    for g, t in [(g1, "Complete Blood Count"), (g2, "Coagulation Panel") , (g3, "Metabolic Panel"), (g4, "Inflammatory Markers"),
                 (g5, "Lipid Panel") ,(g6, "Vitamins"), (g7, "Thyroid Panel")]:
        g.set_title(t)
        g.set_ylabel("")
        g.set_xticklabels(["0 → 1","1 → 2","2 → 3","3 → 4","4 → 5","5 → 6","6 → 7","7 → 8","8 → 9",
                           "9 → 10","10 → 11","11 → 12","12 → 13"], rotation=45, size=7)
        for t in g.texts:
            if abs(float(t.get_text()))>0:
                t.set_text(t.get_text()) #if the abs value is greater than 0 then we show the actual number on the heatmap
            else:
                t.set_text("") # if not it sets an empty text
    fig.suptitle(title, fontsize="x-large")     
    if save_name: # save figure to folder
        fig.get_figure().savefig(figs_folder+save_name+ '.png', format='png',dpi=400)


def plot_SMD_heatmap(df: pd.DataFrame, df_rounded: pd.DataFrame, title: str, 
                       figs_folder: str, save_name=False,v_min=-0.5, v_max=0.5, figsize=(13,18)):
    """Plots the main visualization (heatmap) - SMD of all months (pre and post infection)
    main parameters include the dataframes that includes SMDiD and SMDiD after rounding based on CI
    Labs are considered to be in order for the visualization

    Args:
        df (pd.DataFrame): pd.Dataframe of effects for each month for each laboratory biomarker
        df_rounded (pd.DataFrame): Rounded pd.Dataframe of effects for annotation
        title (str):  Main title of the plot
        figs_folder (str): Path to the figures folder of the relevant experiment
        save_name (bool, optional): Relevant save name, will save in "figs_folder" path. Defaults when no saving.
        v_min (float, optional): Set the color bar range. Defaults to -0.5.
        v_max (float, optional): Set the color bar range. Defaults to 0.5.
        figsize (tuple, optional): Plot size. Defaults to (13,18).
    """
    
    # split into cbc/coagulation/metabolic/inflammatory/lipid/vitamins/thyroid
    hemato_df = df.iloc[HEMATO_RANGE]; hemato_rounded = df_rounded.iloc[HEMATO_RANGE]
    coagulation_df = df.iloc[COAG_RANGE]; coagulation_rounded = df_rounded.iloc[COAG_RANGE]
    metab_df = df.iloc[METAB_RANGE]; metab_rounded = df_rounded.iloc[METAB_RANGE]
    infl_df = df.iloc[INFL_RANGE]; infl_rounded = df_rounded.iloc[INFL_RANGE]
    choles_df = df.iloc[CHOLES_RANGE]; choles_rounded = df_rounded.iloc[CHOLES_RANGE]
    vit_df = df.iloc[VITM_RANGE]; vit_rounded = df_rounded.iloc[VITM_RANGE]
    endo_seru_df = df.iloc[ENDO_RANGE]; endo_rounded = df_rounded.iloc[ENDO_RANGE]
    relative_heights = [] # Ensure all sub-plots will contain rows of the same height
    for df_t in [hemato_df, coagulation_df, metab_df, infl_df, choles_df, vit_df, endo_seru_df]:
        relative_heights.append(df_t.shape[0])
    # create subplots
    fig, axes = plt.subplots(nrows=7,  figsize=figsize,
                             gridspec_kw={'height_ratios':relative_heights}, constrained_layout=True)
    # plot each heatmap separately, annotate based on df_rounded
    g1 = sns.heatmap(hemato_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[0], cmap="vlag", annot=hemato_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                     xticklabels=True, yticklabels=True, cbar=False)
    g2 = sns.heatmap(coagulation_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[1], cmap="vlag", annot=coagulation_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                     xticklabels=True, yticklabels=True, cbar=False)
    g3 = sns.heatmap(metab_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[2], cmap="vlag", annot=metab_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                        xticklabels=True, yticklabels=True, cbar_kws={"shrink": 1})  
    g4 = sns.heatmap(infl_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[3], cmap="vlag", annot=infl_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                     xticklabels=True, yticklabels=True, cbar=False)  
    g5 = sns.heatmap(choles_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[4], cmap="vlag", annot=choles_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                     xticklabels=True, yticklabels=True, cbar=False)    
    g6 = sns.heatmap(vit_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[5], cmap="vlag", annot=vit_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                     xticklabels=True, yticklabels=True, cbar=False)    
    g7 = sns.heatmap(endo_seru_df.reset_index(), vmin=v_min, vmax=v_max,
                     ax=axes[6], cmap="vlag", annot=endo_rounded.round(2), fmt=".3", annot_kws={"fontsize":7}, 
                     xticklabels=True, yticklabels=True, cbar=False)         
    # set title, remove y label, and ensure only text above the visualization threshold is presented in each heatmap
    for g, t in [(g1, "Complete Blood Count"), (g2, "Coagulation Panel") , (g3, "Metabolic Panel"), (g4, "Inflammatory Markers"),
                 (g5, "Lipid Panel") ,(g6, "Vitamins"), (g7, "Thyroid Panel")]:
        g.set_title(t)
        g.set_ylabel("")
        g.axvline(12, c='red', ls='--') # mark the infection time (0) with dashed line
        g.set_xticklabels(["-12 → -11","-11 → -10","-10 → -9","-9 → -8","-8 → -7","-7 → -6","-6 → -5",
                           "-5 → -4","-4 → -3","-3 → -2","-2 → -1","-1 → 0","0 → 1","1 → 2","2 → 3","3 → 4","4 → 5",
                           "5 → 6","6 → 7","7 → 8","8 → 9","9 → 10","10 → 11","11 → 12","12 → 13"], size=7, rotation=45)
        for t in g.texts:
            if abs(float(t.get_text()))>0:
                t.set_text(t.get_text()) #if the abs value is greater than 0 then we show the actual number on the heatmap
            else:
                t.set_text("") # if not it sets an empty text
    fig.suptitle(title, fontsize="x-large")     
    if save_name: # save figure to folder
        fig.get_figure().savefig(figs_folder+save_name+'.png', format='png',dpi=400)        


def visualize_DiD_line_plot(lab_temp_df: pd.DataFrame, lab_name: str, lab_dir: str, figsize=(12,6)):
    """Plots lineplot of mean average values plus DiD correction for a single laboratory

    Args:
        lab_temp_df (pd.DataFrame): laboratory dataframe of a single lab, contains rowid identifier, month, and result.
        lab_name (str): laboratory biomarker name
        lab_dir (str): save path
        figsize (tuple, optional): plot size, defaults to (12,6).
    """
    # get the mean average result of each rowid monthly
    mean_monthly_res = lab_temp_df.groupby(['group','adjusted_lab_month','rowid'])['result'].mean().reset_index()
    # get the mean result of all cohort per group, monthly
    average_values = anl.calculate_diff_in_diff(mean_monthly_res)
    # set up the figure
    plt.figure(figsize=figsize)
    # Plot the mean values of infected and control separately first
    g = sns.lineplot(data=mean_monthly_res,x='adjusted_lab_month', y='result', 
                     hue='group', palette={"Control": 'green', "Infected": 'red'})
    # Plot the DiD correction of the infected
    treated_diff = mean_monthly_res[mean_monthly_res.group=='Infected']
    treated_diff['difference_in_differences'] = treated_diff.result - \
        average_values[average_values.adjusted_lab_month<0].difference_in_differences.mean()
    sns.lineplot(data=treated_diff, x='adjusted_lab_month', y='difference_in_differences', color='red', 
                 linestyle='dashed', label='Infected after DiD correction')
    plt.xticks(range(-12,13))
    g.set_xlim(-12,12)
    plt.margins(x=0)
    plt.title("Monthly Average Results, Lab: {}".format(lab_name))
    plt.xlabel("Month ('0' is infection date)")
    plt.ylabel("Average Result")
    plt.axvline(x=0, color='grey', linestyle='dashed')
    plt.tight_layout()
    # save figure
    g.get_figure().savefig(lab_dir+'lab_{}_DiD_plot'.format(lab_name)+'.png', 
                           format='png',dpi=300)
