
import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import yaml
import numpy as np
import os

def inputs():
    experimentID = 'exp-741748223893'
    date_time_stamp = '2023-05-16 14-29-45'

    # select the kind of plots you want. You need to answer 'Yes' to at least one
    plot_final = 'Yes'     #Answer 'Yes' or 'No'
    plot_final_minus_initial = 'Yes'   #Answer 'Yes' or 'No'

    # parameter extracted from Gaussian Process analysis to plot against HPLC
    # e.g. 'Max Specific Growth Rate (1/Hour)' or 'ln(Max Amplitude)'
    plot_x_label = 'Max Specific Growth Rate (1/Hour)'

    # This is for sources of organic acids other than those detailed in the definitions yaml and the fillstock yamls
    # Units must be Molar (M) !!!! NOT g/L !!!!!
    other_sources = [['1,2-propanediol', 0], ['1,3-propanediol', 0], ['2,3-butanediol', 0], ['acetate', 0], ['butyrate', 0], ['ethanol', 0],
                     ['fructose', 0], ['fumarate', 0], ['galactose', 0], ['glucose', 0], ['glycerol', 0], ['isobutanol', 0], ['isopropanol', 0],
                     ['isovalerate', 0], ['lactate', 0], ['lactose', 0], ['malate', 0], ['n-butanol', 0], ['n-propanol', 0], ['oxaloacetate', 0],
                     ['propionate', 0], ['succinate', 0], ['trehalose', 0], ['valerate', 0], ['xylose', 0]
                     ]

    # This should not need to be changed, but take a look at it and make sure nothing is missing
    # This is to match chemicals that are in the definitions and media yamls to the appropriate ions detected by HPLC
    # The number provides the number of moles of the ion of interest. e.g. for 1 mole of calciumLactatePentahydrate there are 2 moles of lactate
    match_chemical_w_ion = [['sodiumAcetate', 'acetate', 1], ['ammoniumAcetate', 'acetate', 1], ['potassiumAcetate', 'acetate', 1], ['aceticAcid', 'acetate', 1],
                            ['sodiumAcetateTrihydrate', 'acetate', 1], ['sodiumButyrate', 'butyrate', 1], ['butyricAcid', 'butyrate', 1],
                            ['sodiumFumarate', 'fumarate', 1], ['sodiumDLLactate', 'lactate', 1], ['sodiumLlactate', 'lactate', 1],
                            ['L-lacticAcid', 'lactate', 1], ['potassiumLLactate', 'lactate', 1], ['calciumLactatePentahydrate', 'lactate', 2],
                            ['oxaloaceticAcid', 'oxaloacetate', 1], ['sodiumPropionate', 'propionate', 1], ['propionicAcid', 'propionate', 1],
                            ['succinicAcid', 'succinate', 1], ['sodiumValerate', 'valerate', 1]
                            ]

    ####################################################################################
    # Stuff that should not need to be changed
    ####################################################################################

    # name of file containing hplc data
    filename = 'hplc_clean_data.csv'

    last_characters = 'g/L'  # this helps identify the columns that have the amount of the molecule in g/L
    num_last_characters = 3  # number of characters in last_characters

    # 'root' location of HPLC data
    hplc_directory_path = pathlib.Path('/Volumes/data/prod/Experiments/RobotData/HPLC/Results/organic_acid/')
    hplc_data_path = hplc_directory_path / experimentID / date_time_stamp / filename

    # location of the yamls
    yaml_layouts_path = pathlib.Path('/Volumes/data/prod/Experiments/Layouts/')
    yaml_definitions_path = pathlib.Path('/Volumes/data/prod/Experiments/Definitions/')
    yaml_made_path = pathlib.Path('/Volumes/data/prod/Experiments/Made/')

    # biolector analysis location
    bl_location = pathlib.Path('/Volumes/data/prod/Experiments/Analytics/')

    # Dictionary with molecular weights
    conversion_dict = dict(
        [('1,2-propanediol', 76.095), ('1,3-propanediol', 76.09), ('2,3-butanediol', 90.121), ('acetate', 60.052), ('acetoin', 88.11),
         ('butyrate', 88.11), ('ethanol', 46.069), ('fumarate', 116.07), ('fructose', 180.16), ('galactose', 180.156),
         ('glucose', 180.156), ('glycerol', 92.094), ('isobutanol', 74.122), ('isopropanol', 60.1), ('isovalerate', 102.13),
         ('lactate', 90.078), ('lactose', 342.3), ('malate', 134.09), ('n-butanol', 74.123), ('n-propanol', 60.096),
         ('oxaloacetate', 132.07), ('propionate', 74.079), ('succinate', 118.088), ('trehalose', 342.296), ('valerate', 102.133),
         ('xylose', 150.13)
         ])
    #################################################################################

    return plot_final, plot_final_minus_initial, experimentID, yaml_definitions_path, yaml_layouts_path, yaml_made_path, hplc_data_path, last_characters, num_last_characters, conversion_dict, other_sources, match_chemical_w_ion, bl_location, plot_x_label

def get_single_analysis(experiment, analysis_loc):
    files = os.listdir(analysis_loc / experiment)

    for i in range(len(files)):
        if files[i][:8] == 'analysis':
            analysis_filename = files[i]
            break

    analysis = pd.read_csv(analysis_loc / experiment / analysis_filename, keep_default_na=False)

    return analysis

def open_hplc_dataframe(hplc_data_path):
    hplc_df = pd.read_csv(hplc_data_path, keep_default_na=False)
    column_names = list(hplc_df.columns)

    return hplc_df, column_names

def remove_calibration_samples(hplc_df):
    indices_to_drop = []
    for i in range(len(hplc_df)):
        plate = hplc_df['plate_name'].iloc[i]
        if plate[:5] != 'plate':
            indices_to_drop.append(i)
    calibrations_removed_df = hplc_df.drop(indices_to_drop)

    return calibrations_removed_df

def extract_from_yaml(experimentID, yaml_layouts_path):
    layout = open_yaml_expcond(yaml_layouts_path, experimentID, 'layout')

    plates = []
    well_id = []
    inoculant = []
    cond = []
    keys1 = list(layout['Plates'].keys())
    for i in range(len(keys1)):
        keys2 = list(layout['Plates'][keys1[i]]['Wells'].keys())
        for j in range(len(keys2)):
            plates.append(keys1[i].split('.')[1])
            well_id.append(keys2[j])
            inoculant.append(layout['Plates'][keys1[i]]['Wells'][keys2[j]]['inoculant'])
            cond.append(layout['Plates'][keys1[i]]['Wells'][keys2[j]]['variant'])

    df = pd.DataFrame({'plate_name': plates, 'well_id': well_id, 'Inoculant': inoculant, 'Condition': cond})

    return df

def add_condition(hplc_df, experimentID, yaml_layouts_path):
    layout_df = extract_from_yaml(experimentID, yaml_layouts_path)
    hplc_df = pd.merge(hplc_df, layout_df, on=['well_id', 'plate_name'])

    return hplc_df

def stack_columns(headers, hplc_df, last_characters, num_last_characters):
    stacked_hplc_df = pd.DataFrame({'plate_name': [], 'well_id': [], 'Inoculant': [], 'Condition': [], 'Organic Acid': [], 'Detection Method': [] ,'HPLC Measured Concentration (g/L)': []})

    for i in range(len(headers)):
        if headers[i][-num_last_characters:] == last_characters:
            oa = []
            detection_method = []
            for j in range(len(hplc_df)):
                oa.append(headers[i].split('_')[0])
                detection_method.append((headers[i].split('_')[1]))
            df = pd.DataFrame({'plate_name': hplc_df['plate_name'], 'well_id': hplc_df['well_id'], 'Inoculant': hplc_df['Inoculant'], 'Condition': hplc_df['Condition'], 'Organic Acid': oa, 'Detection Method': detection_method, 'HPLC Measured Concentration (g/L)': hplc_df[headers[i]]})
            stacked_hplc_df = pd.concat([stacked_hplc_df, df])

    return stacked_hplc_df

def open_yaml_expcond(file_loc, file, type):
    # open the yaml files, if we are opening the layout yaml file '_layout.yaml' needs to be appended to the experiment ID
    if type == 'layout':
        filename = file + '_layout.yaml'
        file_to_open = file_loc / filename
    elif type == 'made':
        filename = file + '_made.yaml'
        file_to_open = file_loc / filename
    else:
        filename = file + '.yaml'
        file_to_open = file_loc / filename

    with open(file_to_open) as file:
        try:
            data = yaml.safe_load(file)
        except yaml.YAMLError as exception:
            print(exception)

    return data

def expcond(data):
    dict = {}

    variants = data['Variants']
    keys = list(variants.keys())
    for i in range(len(keys)):
        keys2 = list(variants[keys[i]].keys())
        for j in range(len(keys2)):
            keys3 = list(variants[keys[i]][keys2[j]].keys())
            values = list(variants[keys[i]][keys2[j]].values())
            for t in range(len(values)):
                values[t] = str(values[t])
            for k in range(len(keys3)):
                new_value = ''
                for m in range(len(values[k])):
                    new_value += values[k][m]
                    if values[k][m] == ' ':
                        temp = float(new_value)
                        temp = np.round(temp, 3)
                        new_value = str(temp)
                        new_value += values[k][m]

                dict[keys3[k]] = keys[i] + " " + new_value

    return dict

def add_xaxis_info(yaml_definitions_path, experimentID, hplc_df):
    data = open_yaml_expcond(yaml_definitions_path, experimentID, 'definition')
    experiment_type = data['Design']

    if experiment_type == 'Dilution':
        expcond_dict = expcond(data)
        conditions = []
        for i in range(len(hplc_df)):
            conditions.append(hplc_df['Inoculant'].iloc[i] + ' ' + expcond_dict[hplc_df['Condition'].iloc[i]])

    elif experiment_type == 'Desolation':
        expcond_dict = expcond(data)
        conditions = []
        for i in range(len(hplc_df)):
            if hplc_df['Condition'].iloc[i] != 'cond000':
                conditions.append(hplc_df['Inoculant'].iloc[i] + ' ' + expcond_dict[hplc_df['Condition'].iloc[i]])
            else:
                conditions.append(hplc_df['Inoculant'].iloc[i] + ' ' + 'Control')

    else:
        conditions = hplc_df['Condition']

    hplc_df.insert(4, 'Experimental Conditions', conditions)

    return hplc_df

def blank_subtraction_propionate(hplc_df):
    detection_methods = pd.unique(hplc_df['Detection Method'])

    none_df = hplc_df[(hplc_df['Inoculant'] == 'None') & (hplc_df['Organic Acid'] == 'propionate')]
    means = none_df.groupby('Detection Method')['HPLC Measured Concentration (g/L)'].mean()

    blanked_concentration = np.zeros(len(hplc_df))
    for i in range(len(hplc_df)):
        if hplc_df['Organic Acid'].iloc[i] == 'propionate':
            for j in range(len(detection_methods)):
                if hplc_df['Detection Method'].iloc[i] == detection_methods[j]:
                    blanked_concentration[i] = means[detection_methods[j]]
    hplc_df.insert(8, 'Blank (g/L)', blanked_concentration)

    hplc_df['HPLC Measured Concentration (g/L) - Blanked'] = hplc_df['HPLC Measured Concentration (g/L)'] - hplc_df['Blank (g/L)']

    return hplc_df

def get_concentration(value, oa, conversion_dict, multiplication):
    value = str(value)
    if value[-2:] == ' M':  # need to leave space before M otherwise it catches mM, uM and nM
        c = float(value[:-2])
        c *= conversion_dict[oa]
    elif value[-2:] == 'mM':
        c = float(value[:-3]) / 1e3
        c *= conversion_dict[oa]
    elif value[-2:] == 'uM':
        c = float(value[:-3]) / 1e6
        c *= conversion_dict[oa]
    elif value[-2:] == 'nM':
        c = float(value[:-3]) / 1e9
        c *= conversion_dict[oa]
    elif value[-3:] == 'g/L':
        c = float(value[:-4])
    elif value[-3:] == 'g/l':
        c = float(value[:-4])
    else:
        c = 0
        print('units on organic acid in made yaml had either (1) no units or (2) unusual units, value was set to zero')

    c *= multiplication

    return c

def make_dataframe_using_made(experimentID, yaml_made_path):
    made = open_yaml_expcond(yaml_made_path, experimentID, 'made')

    conditions = list(made['Actual'].keys())
    df = pd.DataFrame({'Condition': [], 'Chemical': [], 'Amount': [], 'Ion': []})
    for i in range(len(conditions)):
        chemicals = list(made['Actual'][conditions[i]].keys())
        ions = list(made['Actual'][conditions[i]].keys())
        amount = []
        for j in range(len(chemicals)):
            amount.append(made['Actual'][conditions[i]][chemicals[j]])
        df1 = pd.DataFrame({'Chemical': chemicals, 'Amount': amount, 'Ion': ions})
        df1.insert(0, 'Condition', conditions[i])
        df = pd.concat([df, df1])

    return df

def convert_to_ion(df, match_chemical_w_ion, organic_acids):
    df.insert(4, 'Multiplication', np.ones(len(df)))

    for i in range(len(df)):
        for j in range(len(match_chemical_w_ion)):
            if df['Chemical'].iloc[i] == match_chemical_w_ion[j][0]:
                df.iat[i, 3] = match_chemical_w_ion[j][1]
                df.iat[i, 4] = match_chemical_w_ion[j][2]
                break
    df_filtered = df[df['Ion'].isin(organic_acids)]

    return df_filtered

def convert_to_gperL(df, conversion_dict):
    df.insert(5, 'Concentration (g/L)', np.zeros(len(df)))
    for i in range(len(df)):
        df.iat[i, 5] = get_concentration(df['Amount'].iloc[i], df['Ion'].iloc[i], conversion_dict, df['Multiplication'].iloc[i])

    return df

def add_made_to_hplc_df(hplc_df, made_df):
    hplc_df.insert(10, 'Concentration (g/L) from Made', np.zeros(len(hplc_df)))
    for i in range(len(hplc_df)):
        for j in range(len(made_df)):
            if hplc_df['Condition'].iloc[i] == made_df['Condition'].iloc[j] and hplc_df['Organic Acid'].iloc[i] == \
                    made_df['Ion'].iloc[j]:
                hplc_df.iat[i, 10] = made_df['Concentration (g/L)'].iloc[j]
                break

    return hplc_df

def get_initial_from_made(hplc_df, experimentID, yaml_made_path, conversion_dict, match_chemical_w_ion):
    made_df = make_dataframe_using_made(experimentID, yaml_made_path)
    made_df = convert_to_ion(made_df, match_chemical_w_ion, hplc_df['Organic Acid'].unique())
    made_df = convert_to_gperL(made_df, conversion_dict)
    hplc_df = add_made_to_hplc_df(hplc_df, made_df)

    return hplc_df

def get_cond_from_other_sources(hplc_df, other_sources, conversion_dict):
    conc_other_sources = np.zeros(len(hplc_df))
    for i in range(len(hplc_df)):
        for j in range(len(other_sources)):
            if hplc_df['Organic Acid'].iloc[i] == other_sources[j][0]:
                conc_other_sources[i] = other_sources[j][1]
                conc_other_sources[i] *= conversion_dict[other_sources[j][0]]
                break

    hplc_df.insert(11, 'Concentration (g/L) from Other Sources', conc_other_sources)

    return hplc_df

def get_initial_concentration(hplc_df):
    hplc_df.insert(12, 'Initial Concentration (g/L)', hplc_df['Concentration (g/L) from Made'] + hplc_df['Concentration (g/L) from Other Sources'])

    return hplc_df

def update_HPLC_with_BL(hplc_df, bl_df):

    hplc_df.rename(columns={'plate_name': 'PlateNumber', 'well_id': 'Well Number'}, inplace=True)
    hplc_bl_df = pd.merge(hplc_df, bl_df, on=['PlateNumber', 'Well Number'])

    return hplc_bl_df

def colors_symbols():
    colors = ['blue', 'orange', 'skyblue', 'grey', 'black', 'red', 'green', 'purple', 'yellow', 'peru',
              'blue', 'orange', 'skyblue', 'grey', 'black', 'red', 'green', 'purple', 'yellow', 'peru',
              'blue', 'orange', 'skyblue', 'grey', 'black', 'red', 'green', 'purple', 'yellow', 'peru',
              'blue', 'orange', 'skyblue', 'grey', 'black', 'red', 'green', 'purple', 'yellow', 'peru',
              'blue', 'orange', 'skyblue', 'grey', 'black', 'red', 'green', 'purple', 'yellow', 'peru']
    symbols = ['x', 'v', '^', 'o', '>', '<', 'P', 'h', 'D', 'p', 'd', '*', '8', 'h', '1', '2', '3', '4',
               'x', 'v', '^', 'o', '>', '<', 'P', 'h', 'D', 'p', 'd', '*', '8', 'h', '1', '2', '3', '4',
               'x', 'v', '^', 'o', '>', '<', 'P', 'h', 'D', 'p', 'd', '*', '8', 'h', '1', '2', '3', '4']

    return colors, symbols

def get_final_minus_initial(hplc_df):
    hplc_df.insert(13, 'Final - Initial Concentration (g/L)', hplc_df['HPLC Measured Concentration (g/L) - Blanked'] - hplc_df['Initial Concentration (g/L)'])

    return hplc_df

def add_condition_for_sorting(hplc_df, df):
    condition = []
    for i in range(len(df)):
        for j in range(len(hplc_df)):
            if df['Experimental Conditions'].iloc[i] == hplc_df['Experimental Conditions'].iloc[j]:
                condition.append(hplc_df['Condition'].iloc[j])
                break

    df.insert(0, 'Condition', condition)

    return df

def individual_plots_w_BL(hplc_bl_df, plot_x_label):
    colors, symbols = colors_symbols()

    # get list of organic acids
    oa_tested = pd.unique(hplc_bl_df['Organic Acid'])
    oa_tested.sort()
    oa_number = len(oa_tested)

    for i in range(oa_number):
        subset = hplc_bl_df.loc[(hplc_bl_df['Organic Acid'] == oa_tested[i])]
        detection_method = pd.unique(subset['Detection Method'])
        for j in range(len(detection_method)):
            # set up the figure
            fig, axs = plt.subplots()
            fig.set_figwidth(8.8)
            fig.set_figheight(5.7)

            subset = hplc_bl_df.loc[(hplc_bl_df['Organic Acid'] == oa_tested[i]) & (hplc_bl_df['Detection Method'] == detection_method[j])]

            subset_mean = subset.groupby(['Experimental Conditions'], as_index=False).mean(numeric_only=True)
            subset_std = subset.groupby(['Experimental Conditions'], as_index=False).std(numeric_only=True)

            if subset_mean.empty == False:
                for k in range(len(subset_mean)):
                    leg = subset_mean['Experimental Conditions'][k]
                    if pd.isna(subset_std['HPLC Measured Concentration (g/L) - Blanked'].iloc[k]) == False:
                        axs.scatter(subset_mean[plot_x_label][k], subset_mean['HPLC Measured Concentration (g/L) - Blanked'][k],
                                    linestyle='None', marker=symbols[k], color=colors[k], label=leg)
                        axs.errorbar(subset_mean[plot_x_label][k], subset_mean['HPLC Measured Concentration (g/L) - Blanked'][k],
                                     xerr=subset_std[plot_x_label][k],
                                     yerr=subset_std['HPLC Measured Concentration (g/L) - Blanked'][k], linestyle='None',
                                     color=colors[k])
                    else:
                        axs.scatter(subset_mean[plot_x_label][k], subset_mean['HPLC Measured Concentration (g/L) - Blanked'][k],
                                        linestyle='None', marker=symbols[k], color=colors[k], label=leg)

            axs.set_ylabel('HPLC: Concentration (g/L) - Blanked')
            axs.set_xlabel(plot_x_label)
            axs.legend(bbox_to_anchor=(1, 1), loc='upper left')

            plt.subplots_adjust(top=0.95, bottom=0.25, left=0.09, right=0.6)
            plt.title(oa_tested[i] + ' ' + detection_method[j])
            plt.show()

def individual_plots_final_minus_initial_w_BL(hplc_bl_df, plot_x_label):
    colors, symbols = colors_symbols()

    # get list of organic acids
    oa_tested = pd.unique(hplc_bl_df['Organic Acid'])
    oa_tested.sort()
    oa_number = len(oa_tested)

    for i in range(oa_number):
        subset = hplc_bl_df.loc[(hplc_bl_df['Organic Acid'] == oa_tested[i])]
        detection_method = pd.unique(subset['Detection Method'])
        for j in range(len(detection_method)):
            # set up the figure
            fig, axs = plt.subplots()
            fig.set_figwidth(8.8)
            fig.set_figheight(5.7)

            subset = hplc_bl_df.loc[(hplc_bl_df['Organic Acid'] == oa_tested[i]) & (hplc_bl_df['Detection Method'] == detection_method[j])]

            subset_mean = subset.groupby(['Experimental Conditions'], as_index=False).mean(numeric_only=True)
            subset_std = subset.groupby(['Experimental Conditions'], as_index=False).std(numeric_only=True)

            if subset_mean.empty == False:
                for k in range(len(subset_mean)):
                    leg = subset_mean['Experimental Conditions'][k]
                    if pd.isna(subset_std['Final - Initial Concentration (g/L)'].iloc[k]) == False:
                        axs.scatter(subset_mean[plot_x_label][k], subset_mean['Final - Initial Concentration (g/L)'][k],
                                    linestyle='None', marker=symbols[k], color=colors[k], label=leg)
                        axs.errorbar(subset_mean[plot_x_label][k], subset_mean['Final - Initial Concentration (g/L)'][k],
                                     xerr=subset_std[plot_x_label][k],
                                     yerr=subset_std['Final - Initial Concentration (g/L)'][k], linestyle='None',
                                     color=colors[k])
                    else:
                        axs.scatter(subset_mean[plot_x_label][k], subset_mean['Final - Initial Concentration (g/L)'][k],
                                        linestyle='None', marker=symbols[k], color=colors[k], label=leg)

            axs.set_ylabel('Concentration (g/L): Final - Initial')
            axs.set_xlabel(plot_x_label)
            axs.legend(bbox_to_anchor=(1, 1), loc='upper left')
            axs.axhline(0, color='black', linestyle='dotted')

            plt.subplots_adjust(top=0.95, bottom=0.25, left=0.09, right=0.6)
            plt.title(oa_tested[i] + ' ' + detection_method[j])
            plt.show()

def save(hplc_df, experimentID):
    hplc_df.to_csv(experimentID + '.csv')

def main():
    # read in inputs
    plot_final, plot_final_minus_initial, experimentID, yaml_definitions_path, yaml_layouts_path, yaml_made_path, hplc_data_path, last_characters, num_last_characters, conversion_dict, other_sources, match_chemical_w_ion, bl_location, plot_x_label = inputs()
    # get HPLC data
    hplc_df, column_names = open_hplc_dataframe(hplc_data_path)
    # remove calibration samples
    hplc_df = remove_calibration_samples(hplc_df)
    # add condition to data frame
    hplc_df = add_condition(hplc_df, experimentID, yaml_layouts_path)
    # get organic acids & stack columns
    hplc_df = stack_columns(column_names, hplc_df, last_characters, num_last_characters)
    # add experimental conditions - details, e.g. glycerol 20 mM not cond001
    hplc_df = add_xaxis_info(yaml_definitions_path, experimentID, hplc_df)
    # Blank subtraction for propionate
    hplc_df = blank_subtraction_propionate(hplc_df)
    # get concentrations from made
    hplc_df = get_initial_from_made(hplc_df, experimentID, yaml_made_path, conversion_dict, match_chemical_w_ion)
    # get concentrations from other sources
    hplc_df = get_cond_from_other_sources(hplc_df, other_sources, conversion_dict)
    # get initial concentrations of organic acids
    hplc_df = get_initial_concentration(hplc_df)
    # read BioLector analysis
    bl_df = get_single_analysis(experimentID, bl_location)
    # update HPLC data with BioLector data
    hplc_bl_df = update_HPLC_with_BL(hplc_df, bl_df)
    # Plot stuff
    if plot_final == 'Yes':
        # individual plots final and initial
        # individual plots
        individual_plots_w_BL(hplc_bl_df, plot_x_label)
    if plot_final_minus_initial == 'Yes':
        # get final minus initial concentration
        hplc_bl_df = get_final_minus_initial(hplc_bl_df)
        # individual plots final minus initial
        individual_plots_final_minus_initial_w_BL(hplc_bl_df, plot_x_label)
    # Save data frame
    save(hplc_bl_df, experimentID)

main()