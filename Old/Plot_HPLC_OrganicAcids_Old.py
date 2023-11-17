
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import yaml
import numpy as np
import math

def inputs():
    experimentID = 'exp-741748223893'
    date_time_stamp = '2023-05-16 14-29-45'

    # select the kind of plots you want. You need to answer 'Yes' to at least one
    plot_final_and_initial = 'Yes'     #Answer 'Yes' or 'No'
    plot_final_minus_initial = 'Yes'   #Answer 'Yes' or 'No'

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
    # name of file containing hplc data
    filename = 'hplc_clean_data.csv'

    last_characters = 'g/L'  # this helps identify the columns that have the amount of the molecule in g/L
    num_last_characters = 3  # number of characters in last_characters

    # 'root' location of HPLC data
    hplc_directory_path = pathlib.Path('/Volumes/data/prod/Experiments/RobotData/HPLC/Results/organic_acid/')
    #hplc_directory_path = pathlib.Path('/Users/jennifer/Desktop/HPLC_testing/organic_acids/')
    hplc_data_path = hplc_directory_path / experimentID / date_time_stamp / filename

    # location of the yamls
    yaml_layouts_path = pathlib.Path('/Volumes/data/prod/Experiments/Layouts/')
    #yaml_layouts_path = pathlib.Path('/Users/jennifer/Desktop/HPLC_testing/yamls/Layouts/')
    yaml_definitions_path = pathlib.Path('/Volumes/data/prod/Experiments/Definitions/')
    #yaml_definitions_path = pathlib.Path('/Users/jennifer/Desktop/HPLC_testing/yamls/Definitions/')
    yaml_fillstock_path = pathlib.Path('/Volumes/data/prod/Chemistry/Combinations/')
    #yaml_fillstock_path = pathlib.Path('/Users/jennifer/Desktop/HPLC_testing/yamls/')

    # Colors
    colors = ['blue', 'grey', 'orange', 'skyblue',  'purple']

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

    return plot_final_and_initial, plot_final_minus_initial, colors, experimentID, yaml_definitions_path, yaml_layouts_path, yaml_fillstock_path, hplc_data_path, last_characters, num_last_characters, conversion_dict, other_sources, match_chemical_w_ion

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
    else:
        filename = file + '.yaml'
        file_to_open = file_loc / filename

    with open(file_to_open) as file:
        try:
            data = yaml.safe_load(file)
        except yaml.YAMLError as exception:
            print(exception)

    return data

def expcond(file_loc, experiment):
    data = open_yaml_expcond(file_loc, experiment, 'definition')
    experiment_type = data['Design']
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

    return dict, experiment_type

def add_xaxis_info(yaml_definitions_path, experimentID, hplc_df):
    expcond_dict, experiment_type = expcond(yaml_definitions_path, experimentID)
    conditions = []

    if experiment_type == 'Dilution':
        for i in range(len(hplc_df)):
            conditions.append(hplc_df['Inoculant'].iloc[i] + ' ' + expcond_dict[hplc_df['Condition'].iloc[i]])
    else:
        for i in range(len(hplc_df)):
            if hplc_df['Condition'].iloc[i] != 'cond000':
                conditions.append(hplc_df['Inoculant'].iloc[i] + ' ' + expcond_dict[hplc_df['Condition'].iloc[i]])
            else:
                conditions.append(hplc_df['Inoculant'].iloc[i] + ' ' + 'Control')

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

def get_concentration(value, oa, conversion_dict):
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
        print('units on organic acid in fillstock or definition yaml did not have any units attached, value was set to zero')

    return c

def convert_to_ion(chemicals_from_yaml, match_chemical_w_ion):
    chemicals_converted = chemicals_from_yaml.copy()
    multiplication = np.ones(len(chemicals_from_yaml))
    for i in range(len(chemicals_from_yaml)):
        for j in range(len(match_chemical_w_ion)):
            if chemicals_from_yaml[i] == match_chemical_w_ion[j][0]:
                chemicals_converted[i] = match_chemical_w_ion[j][1]
                multiplication[i] = match_chemical_w_ion[j][2]
                break

    return chemicals_converted, multiplication

def get_initial_from_fillstock(hplc_df, experimentID, yaml_definitions_path, yaml_fillstock_path, conversion_dict, match_chemical_w_ion):
    definitions_yaml = open_yaml_expcond(yaml_definitions_path, experimentID, 'definitions')

    if definitions_yaml['Fillstock'][:3] == 'Med':
        fillstock_yaml = open_yaml_expcond(yaml_fillstock_path / 'Media', definitions_yaml['Fillstock'], 'fillstock')
    elif definitions_yaml['Fillstock'][:3] == 'Fee':
        fillstock_yaml = open_yaml_expcond(yaml_fillstock_path / 'Feedstock', definitions_yaml['Fillstock'], 'fillstock')
    elif definitions_yaml['Fillstock'][:3] == 'Oth':
        fillstock_yaml = open_yaml_expcond(yaml_fillstock_path / 'Other', definitions_yaml['Fillstock'], 'fillstock')
    else:
        print('fillstock file cannot be found, not in Media, Feedstock or Other folder')

    dilution = definitions_yaml['Fillstock Dilution']
    initial_conc = np.zeros(len(hplc_df))

    chemicals_in_fillstock = list(fillstock_yaml['components'].keys())
    chemicals_in_fillstock_converted_to_ion, multiplication = convert_to_ion(chemicals_in_fillstock, match_chemical_w_ion)
    for i in range(len(hplc_df)):
        for j in range(len(chemicals_in_fillstock_converted_to_ion)):
            if hplc_df['Organic Acid'].iloc[i] == chemicals_in_fillstock_converted_to_ion[j]:
                c = get_concentration(fillstock_yaml['components'][chemicals_in_fillstock[j]]['solvation'], hplc_df['Organic Acid'].iloc[i], conversion_dict)
                initial_conc[i] += c * dilution * multiplication[j]
                # there is no break here, one could have two sources of an ion, e.g. sodiumAcetate and potassiumAcetate

    hplc_df.insert(10, 'Concentration (g/L) from Fillstock', initial_conc)

    return hplc_df

def get_mediafeedstockother_yaml(mediafeedstockother, yaml_path):
    yaml = []

    if mediafeedstockother[:3] == 'Med':
        yaml = open_yaml_expcond(yaml_path / 'Media', mediafeedstockother, 'fillstock')
    elif mediafeedstockother[:3] == 'Fee':
        yaml = open_yaml_expcond(yaml_path / 'Feedstock', mediafeedstockother, 'fillstock')
    elif mediafeedstockother[:3] == 'Oth':
        yaml = open_yaml_expcond(yaml_path / 'Other', mediafeedstockother, 'fillstock')
    else:
        print('File did not start with Med, Fed or Oth')

    return yaml

def replace_control_concentration(df, chemical, cond, variant_concentration):
    for i in range(len(df)):
        if df['Condition'].iloc[i] == cond and df['Chemicals'].iloc[i] == chemical:
            df['Concentration'].iloc[i] = variant_concentration
            break

    return df

def get_conc_from_definitions(hplc_df, experimentID, yaml_definitions_path, match_chemical_w_ion, conversion_dict, yaml_fillstock_path):
    definitions_yaml = open_yaml_expcond(yaml_definitions_path, experimentID, 'definitions')
    if definitions_yaml['Design'] == 'Desolation':
        get_conc_from_definitions_desolation(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict)
    elif definitions_yaml['Design'] == 'Dilution':
        get_conc_from_definitions_dilution(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict, yaml_fillstock_path)
    else:
        print('Definition yaml is not a Desolation or Dilution - no initial values have been entered')
        hplc_df.insert(11, 'Concentration (g/L) from Definitions', np.zeros(len(hplc_df)))

    return hplc_df

def get_conc_from_definitions_desolation(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict):
    # get list of conditions
    conditions = hplc_df['Condition'].unique()
    conditions.sort()

    # get list of chemicals
    chemicals = list(definitions_yaml['Variants'].keys())
    chemicals.sort()

    # create data frame
    conditions_array = []
    chemicals_array = []
    for i in range(len(chemicals)):
        for j in range(len(conditions)):
            conditions_array.append(conditions[j])
            chemicals_array.append(chemicals[i])
    df = pd.DataFrame({'Condition': conditions_array, 'Chemicals': chemicals_array})

    # get control concentrations - going to put these in every space, then will overwrite specific ones with 'variant' concentrations below
    control_concentrations = []
    for i in range(len(df)):
        control_concentrations.append(definitions_yaml['Variants'][df['Chemicals'].iloc[i]]['Control']['cond000'])
    df.insert(2, 'Concentration', control_concentrations)

    # get variant concentrations and overwrite control concentrations where appropriate
    for i in range(len(chemicals)):
        cond = list(definitions_yaml['Variants'][chemicals[i]]['Variants'].keys())
        for j in range(len(cond)):
            variant_concentration = definitions_yaml['Variants'][chemicals[i]]['Variants'][cond[j]]
            df = replace_control_concentration(df, chemicals[i], cond[j], variant_concentration)

    # convert chemicals to ions
    chemicals_converted_to_ion, multiplication = convert_to_ion(chemicals_array, match_chemical_w_ion)
    df.insert(3, 'Ions', chemicals_converted_to_ion)
    df.insert(4, 'Multiplication', multiplication)

    # convert concentrations to g/L
    concentrations_g_per_L = np.zeros(len(df))
    for i in range(len(df)):
        concentrations_g_per_L[i] = get_concentration(df['Concentration'].iloc[i], df['Ions'].iloc[i], conversion_dict)
        concentrations_g_per_L[i] *= multiplication[i]
    df.insert(5, 'Concentration (g/L) From Definitions', concentrations_g_per_L)

    # Add the concentrations from the definitions yaml to the data frame
    conc_from_def = np.zeros(len(hplc_df))
    for i in range(len(hplc_df)):
        for j in range(len(df)):
            if hplc_df['Condition'].iloc[i] == df['Condition'].iloc[j] and hplc_df['Organic Acid'].iloc[i] == df['Ions'].iloc[j]:
                conc_from_def[i] += df['Concentration (g/L) From Definitions'].iloc[j]
                break

    hplc_df.insert(11, 'Concentration (g/L) from Definitions', conc_from_def)

    return hplc_df

def get_conc_from_definitions_dilution(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict, yaml_path):
    variants = list(definitions_yaml['Variants'].keys())
    num_variants = len(variants)
    count = 0
    for i in range(num_variants):
        if variants[i][:3] == 'Med' or  variants[i][:3] == 'Fee' or variants[i][:3] == 'Oth':
            count += 1

    if count == num_variants:
        hplc_df = get_conc_from_definitions_dilution_mediafeedstocksother(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict, yaml_path)
    elif count == 0:
        hplc_df = get_conc_from_definitions_dilution_chemicals(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict)
    else:
        hplc_df.insert(11, 'Concentration (g/L) from Definitions', np.zeros(len(hplc_df)))
        print('!!! the definitions yaml is of a mixed type and could not be handled, values have all been set to zero')

    return hplc_df

def get_conc_from_definitions_dilution_mediafeedstocksother(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict, yaml_path):
    organic_acids = hplc_df['Organic Acid'].unique()
    mediafeedstockother = list(definitions_yaml['Variants'].keys())
    condition_array_bigloop = []
    chemcial_array_bigloop = []
    concentration_array_bigloop = []
    for i in range(len(mediafeedstockother)):
        mediafeedstockother_yaml = get_mediafeedstockother_yaml(mediafeedstockother[i], yaml_path)
        chemicals_in_mediafeedstockother = list(mediafeedstockother_yaml['components'].keys())
        chemicals_in_mediafeedstockother_converted_to_ion, multiplication_converted_to_ion = convert_to_ion(chemicals_in_mediafeedstockother, match_chemical_w_ion)
        solvation_array_small_loop = []
        for j in range(len(chemicals_in_mediafeedstockother)):
            solvation_array_small_loop.append(mediafeedstockother_yaml['components'][chemicals_in_mediafeedstockother[j]]['solvation'])

        for j in range(len(organic_acids)):
            # get the controls
            concentration = 0
            cond = list(definitions_yaml['Variants'][mediafeedstockother[i]]['Control'].keys())[0]
            for k in range(len(chemicals_in_mediafeedstockother_converted_to_ion)):
                if organic_acids[j] == chemicals_in_mediafeedstockother_converted_to_ion[k]:
                    temp = get_concentration(solvation_array_small_loop[k], organic_acids[j], conversion_dict)
                    temp *= multiplication_converted_to_ion[k]
                    temp *= float(definitions_yaml['Variants'][mediafeedstockother[i]]['Control'][cond])
                    concentration += temp
                    # no break b/c there might be more than one chemical corresponding to organic acid, e.g. sodiumAcetate & potassiumAcetate both become acetate
            condition_array_bigloop.append(cond)
            chemcial_array_bigloop.append(organic_acids[j])
            concentration_array_bigloop.append(concentration)

            # get the variants
            cond = list(definitions_yaml['Variants'][mediafeedstockother[i]]['Variants'].keys())
            for m in range(len(cond)):
                concentration = 0
                for k in range(len(chemicals_in_mediafeedstockother_converted_to_ion)):
                    if organic_acids[j] == chemicals_in_mediafeedstockother_converted_to_ion[k]:
                        temp = get_concentration(solvation_array_small_loop[k], organic_acids[j], conversion_dict)
                        temp *= multiplication_converted_to_ion[k]
                        temp *= float(definitions_yaml['Variants'][mediafeedstockother[i]]['Variants'][cond[m]])
                        concentration += temp
                        # no break b/c there might be more than one chemical corresponding to organic acid, e.g. sodiumAcetate & potassiumAcetate both become acetate
                condition_array_bigloop.append(cond[m])
                chemcial_array_bigloop.append(organic_acids[j])
                concentration_array_bigloop.append(concentration)

    df = pd.DataFrame({'Condition': condition_array_bigloop, 'Chemicals': chemcial_array_bigloop,
                       'Concentration (g/L) From Definitions': concentration_array_bigloop})

    # Add the concentrations from the definitions yaml to the data frame
    conc_from_def = np.zeros(len(hplc_df))
    for i in range(len(hplc_df)):
        for j in range(len(df)):
            if hplc_df['Condition'].iloc[i] == df['Condition'].iloc[j] and hplc_df['Organic Acid'].iloc[i] == df['Chemicals'].iloc[j]:
                conc_from_def[i] += df['Concentration (g/L) From Definitions'].iloc[j]
                break

    hplc_df.insert(11, 'Concentration (g/L) from Definitions', conc_from_def)

    return hplc_df

def get_conc_from_definitions_dilution_chemicals(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict):
    chemicals = list(definitions_yaml['Variants'].keys())
    chemicals_converted_to_ion, multiplication = convert_to_ion(chemicals, match_chemical_w_ion)
    organic_acids = hplc_df['Organic Acid'].unique()

    chemicals_array = []
    concentration_control = []
    condition_array = []
    for i in range(len(chemicals)):
        for j in range(len(organic_acids)):
            if chemicals_converted_to_ion[i] == organic_acids[j]:
                cond = list(definitions_yaml['Variants'][chemicals[i]]['Control'].keys())[0]
                c = definitions_yaml['Variants'][chemicals[i]]['Control'][cond]
                c = get_concentration(c, organic_acids[j], conversion_dict)
                c *= multiplication[i]

                condition_array.append(cond)
                concentration_control.append(c)
                chemicals_array.append(chemicals_converted_to_ion[i])

                cond = list(definitions_yaml['Variants'][chemicals[i]]['Variants'].keys())
                for k in range(len(cond)):
                    c = definitions_yaml['Variants'][chemicals[i]]['Variants'][cond[k]]
                    c = get_concentration(c, organic_acids[j], conversion_dict)
                    c *= multiplication[i]

                    condition_array.append(cond[k])
                    concentration_control.append(c)
                    chemicals_array.append(chemicals_converted_to_ion[i])

                break

    df = pd.DataFrame({'Condition': condition_array, 'Chemicals': chemicals_array,
                        'Concentration (g/L) From Definitions': concentration_control})

    # Add the concentrations from the definitions yaml to the data frame
    conc_from_def = np.zeros(len(hplc_df))
    for i in range(len(hplc_df)):
        for j in range(len(df)):
            if hplc_df['Condition'].iloc[i] == df['Condition'].iloc[j] and hplc_df['Organic Acid'].iloc[i] == df['Chemicals'].iloc[j]:
                conc_from_def[i] += df['Concentration (g/L) From Definitions'].iloc[j]
                break

    hplc_df.insert(11, 'Concentration (g/L) from Definitions', conc_from_def)

    return hplc_df

def get_cond_from_other_sources(hplc_df, other_sources, conversion_dict):
    conc_other_sources = np.zeros(len(hplc_df))
    for i in range(len(hplc_df)):
        for j in range(len(other_sources)):
            if hplc_df['Organic Acid'].iloc[i] == other_sources[j][0]:
                conc_other_sources[i] = other_sources[j][1]
                conc_other_sources[i] *= conversion_dict[other_sources[j][0]]
                break

    hplc_df.insert(12, 'Concentration (g/L) from Other Sources', conc_other_sources)

    return hplc_df

def get_initial_concentration(hplc_df):
    hplc_df.insert(13, 'Initial Concentration (g/L)', hplc_df['Concentration (g/L) from Definitions'] + hplc_df['Concentration (g/L) from Fillstock'] + hplc_df['Concentration (g/L) from Other Sources'])

    return hplc_df

def get_final_minus_initial(hplc_df):
    hplc_df.insert(14, 'Final - Initial Concentration (g/L)', hplc_df['HPLC Measured Concentration (g/L) - Blanked'] - hplc_df['Initial Concentration (g/L)'])

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

def make_color_dictionary(values, hplc_df):
    keys = pd.unique(hplc_df['Detection Method'])
    color_dictionary = dict(zip(keys, values))

    return color_dictionary

def grid_plot(hplc_df, color_dictionary):
    # get list of organic acids
    oa_tested = pd.unique(hplc_df['Organic Acid'])
    oa_tested.sort()

    # find number of rows for plot
    oa_number = len(oa_tested)
    num_y = math.floor(oa_number/5)
    remainder = oa_number % 5
    if remainder != 0 :
        num_y += 1

    # set up the figure
    fig, axs = plt.subplots(num_y, 5, sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0.4)
    fig.set_figwidth(12)
    fig.set_figheight(num_y * 3)

    x = 0; y= 0

    for i in range(oa_number):
        subset = hplc_df.loc[(hplc_df['Organic Acid'] == oa_tested[i])]
        detection_method = pd.unique(subset['Detection Method'])

        if len(detection_method) == 1:
            multiplier = 1
            width = 0.30
        elif len(detection_method) == 2:
            multiplier = 0.5
            width = 0.25
        elif len(detection_method) == 3:
            multiplier = 0
            width = 0.2
        elif len(detection_method) == 4:
            multiplier = -0.5
            width = 0.15
        else:
            multiplier = -1
            width = 0.1

        for j in range(len(detection_method)):
            subset_of_subset = subset.loc[(subset['Detection Method'] == detection_method[j])]

            subset_mean = subset_of_subset.groupby(['Experimental Conditions'], as_index=False).mean(numeric_only=True)
            subset_mean = add_condition_for_sorting(hplc_df, subset_mean)
            subset_mean.sort_values(by=['Condition'], inplace=True)

            subset_std = subset_of_subset.groupby(['Experimental Conditions'], as_index=False).std(numeric_only=True)
            subset_std = add_condition_for_sorting(hplc_df, subset_std)
            subset_std.sort_values(by=['Condition'], inplace=True)

            xx = np.arange(len(subset_mean['Experimental Conditions']))
            subset_mean.insert(0, 'x-axis', xx)

            offset = width * multiplier

            if subset_mean.empty == False:
                if pd.isna(subset_std['HPLC Measured Concentration (g/L) - Blanked'].iloc[0]) == False:
                    axs[y, x].scatter(subset_mean['x-axis'] + offset, subset_mean['HPLC Measured Concentration (g/L) - Blanked'], marker='s', color=color_dictionary[detection_method[j]])
                    axs[y, x].scatter(subset_mean['x-axis'] + offset, subset_mean['Initial Concentration (g/L)'], marker='o', color=color_dictionary[detection_method[j]])
                    axs[y, x].errorbar(subset_mean['x-axis'] + offset, subset_mean['HPLC Measured Concentration (g/L) - Blanked'], yerr=subset_std['HPLC Measured Concentration (g/L) - Blanked'], color=color_dictionary[detection_method[j]], linestyle='None')
                    axs[y, x].errorbar(subset_mean['x-axis'] + offset, subset_mean['Initial Concentration (g/L)'], yerr=subset_std['Initial Concentration (g/L)'], color=color_dictionary[detection_method[j]], linestyle='None')
                else:
                    axs[y, x].scatter(xx + offset, subset_mean['HPLC Measured Concentration (g/L) - Blanked'], marker='s', color=color_dictionary[detection_method[j]])
                    axs[y, x].scatter(xx + offset, subset_mean['Initial Concentration (g/L)'], marker='o', color=color_dictionary[detection_method[j]])

            multiplier += 1

        axs[y, x].set_title(oa_tested[i])
        axs[y, x].tick_params(axis='x', labelrotation = 90)
        axs[y, x].set_xticks(subset_mean['x-axis'] + width, subset_mean['Experimental Conditions'])

        x += 1
        if x > 4:
            x = 0
            y += 1

    if x > 0:
        for i in range(5 - x):
            axs[y, x].tick_params(axis='x', labelrotation=90)

            x += 1

    fig.supylabel('Concentration (g/L)')
    fig.suptitle('Circle: Initial Concentration (g/L) \n Square: Concentration Measured by HPLC (g/L) \n ')
    plt.subplots_adjust(top=0.9, bottom=0.32, left=0.09, right=0.9)

    patchList = []
    for key in color_dictionary:
        data_key = mpatches.Patch(color=color_dictionary[key], label=key)
        patchList.append(data_key)
    axs[0, 4].legend(handles=patchList, bbox_to_anchor=(1, 1), loc='upper left')

    plt.show()

def grid_plot_final_minus_initial(hplc_df, color_dictionary):
    # get list of organic acids
    oa_tested = pd.unique(hplc_df['Organic Acid'])
    oa_tested.sort()

    # find number of rows for plot
    oa_number = len(oa_tested)
    num_y = math.floor(oa_number/5)
    remainder = oa_number % 5
    if remainder != 0 :
        num_y += 1

    # set up the figure
    fig, axs = plt.subplots(num_y, 5, sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0.4)
    fig.set_figwidth(12)
    fig.set_figheight(num_y * 3)

    x = 0; y= 0

    for i in range(oa_number):
        subset = hplc_df.loc[(hplc_df['Organic Acid'] == oa_tested[i])]
        detection_method = pd.unique(subset['Detection Method'])

        if len(detection_method) == 1:
            multiplier = 1
            width = 0.30
        elif len(detection_method) == 2:
            multiplier = 0.5
            width = 0.25
        elif len(detection_method) == 3:
            multiplier = 0
            width = 0.2
        elif len(detection_method) == 4:
            multiplier = -0.5
            width = 0.15
        else:
            multiplier = -1
            width = 0.1
        for j in range(len(detection_method)):
            subset_of_subset = subset.loc[(subset['Detection Method'] == detection_method[j])]

            subset_mean = subset_of_subset.groupby(['Experimental Conditions'], as_index=False).mean(numeric_only=True)
            subset_mean = add_condition_for_sorting(hplc_df, subset_mean)
            subset_mean.sort_values(by=['Condition'], inplace=True)

            subset_std = subset_of_subset.groupby(['Experimental Conditions'], as_index=False).std(numeric_only=True)
            subset_std = add_condition_for_sorting(hplc_df, subset_std)
            subset_std.sort_values(by=['Condition'], inplace=True)

            xx = np.arange(len(subset_mean['Experimental Conditions']))
            subset_mean.insert(0, 'x-axis', xx)

            offset = width * multiplier

            if subset_mean.empty == False:
                if pd.isna(subset_std['Final - Initial Concentration (g/L)'].iloc[0]) == False:
                    axs[y, x].bar(subset_mean['x-axis'] + offset, subset_mean['Final - Initial Concentration (g/L)'], width, yerr=subset_std['Final - Initial Concentration (g/L)'], color=color_dictionary[detection_method[j]])
                else:
                    axs[y, x].bar(subset_mean['x-axis'] + offset, subset_mean['Final - Initial Concentration (g/L)'], width, color=color_dictionary[detection_method[j]])

            multiplier += 1

        axs[y, x].set_title(oa_tested[i])
        axs[y, x].tick_params(axis='x', labelrotation = 90)
        axs[y, x].set_xticks(subset_mean['x-axis'] + width, subset_mean['Experimental Conditions'])
        axs[y, x].axhline(0, color='black', linestyle='dotted')

        x += 1
        if x > 4:
            x = 0
            y += 1

    if x > 0:
        for i in range(5 - x):
            axs[y, x].tick_params(axis='x', labelrotation=90)

            x += 1

    fig.supylabel('Concentration (g/L): Final - Initial')
    plt.subplots_adjust(top=0.9, bottom=0.32, left=0.09, right=0.9)

    patchList = []
    for key in color_dictionary:
        data_key = mpatches.Patch(color=color_dictionary[key], label=key)
        patchList.append(data_key)
    axs[0, 4].legend(handles=patchList, bbox_to_anchor=(1, 1), loc='upper left')

    plt.show()

def individual_plots(hplc_df, color_dictionary):
    # get list of organic acids
    oa_tested = pd.unique(hplc_df['Organic Acid'])
    oa_tested.sort()
    oa_number = len(oa_tested)

    for i in range(oa_number):
        # set up the figure
        fig, axs = plt.subplots(1, 1)
        fig.set_figwidth(9)
        fig.set_figheight(6)

        subset = hplc_df.loc[(hplc_df['Organic Acid'] == oa_tested[i])]
        detection_method = pd.unique(subset['Detection Method'])

        if len(detection_method) == 1:
            multiplier = 1
            width = 0.30
        elif len(detection_method) == 2:
            multiplier = 0.5
            width = 0.25
        elif len(detection_method) == 3:
            multiplier = 0
            width = 0.2
        elif len(detection_method) == 4:
            multiplier = -0.5
            width = 0.15
        else:
            multiplier = -1
            width = 0.1
        for j in range(len(detection_method)):
            subset_of_subset = subset.loc[(subset['Detection Method'] == detection_method[j])]

            subset_mean = subset_of_subset.groupby(['Experimental Conditions'], as_index=False).mean(numeric_only=True)
            subset_mean = add_condition_for_sorting(hplc_df, subset_mean)
            subset_mean.sort_values(by=['Condition'], inplace=True)

            subset_std = subset_of_subset.groupby(['Experimental Conditions'], as_index=False).std(numeric_only=True)
            subset_std = add_condition_for_sorting(hplc_df, subset_std)
            subset_std.sort_values(by=['Condition'], inplace=True)

            xx = np.arange(len(subset_mean['Experimental Conditions']))
            subset_mean.insert(0, 'x-axis', xx)

            offset = width * multiplier

            if subset_mean.empty == False:
                if pd.isna(subset_std['HPLC Measured Concentration (g/L) - Blanked'].iloc[0]) == False:
                    axs.scatter(subset_mean['x-axis'] + offset, subset_mean['HPLC Measured Concentration (g/L) - Blanked'], marker='s', color=color_dictionary[detection_method[j]])
                    axs.scatter(subset_mean['x-axis'] + offset, subset_mean['Initial Concentration (g/L)'], marker='o', color=color_dictionary[detection_method[j]])
                    axs.errorbar(subset_mean['x-axis'] + offset, subset_mean['HPLC Measured Concentration (g/L) - Blanked'], yerr=subset_std['HPLC Measured Concentration (g/L) - Blanked'], color=color_dictionary[detection_method[j]], linestyle='None')
                    axs.errorbar(subset_mean['x-axis'] + offset, subset_mean['Initial Concentration (g/L)'], yerr=subset_std['Initial Concentration (g/L)'], color=color_dictionary[detection_method[j]], linestyle='None')
                else:
                    axs.scatter(subset_mean['x-axis'] + offset, subset_mean['HPLC Measured Concentration (g/L) - Blanked'], marker='s', color=color_dictionary[detection_method[j]])
                    axs.scatter(subset_mean['x-axis'] + offset, subset_mean['Initial Concentration (g/L)'], marker='o', color=color_dictionary[detection_method[j]])

            multiplier += 1

        axs.set_title('Circle: Initial Concentration (g/L) \n Square: Concentration Measured by HPLC (g/L) \n \n' + oa_tested[i])
        axs.tick_params(axis='x', labelrotation=90)
        axs.set_xticks(subset_mean['x-axis'] + width, subset_mean['Experimental Conditions'])

        fig.supylabel('Concentration (g/L)')
        plt.subplots_adjust(top=0.85, bottom=0.45, left=0.2, right=0.8)

        patchList = []
        for key in color_dictionary:
            data_key = mpatches.Patch(color=color_dictionary[key], label=key)
            patchList.append(data_key)
        axs.legend(handles=patchList, bbox_to_anchor=(1, 1), loc='upper left')

        plt.show()
        # plt.savefig('/Users/jennifer/Desktop/' + str(i))

def individual_plots_final_minus_initial(hplc_df, color_dictionary):
    color_dictionary['Initial Conc.'] = 'red'

    # get list of organic acids
    oa_tested = pd.unique(hplc_df['Organic Acid'])
    oa_tested.sort()
    oa_number = len(oa_tested)

    for i in range(oa_number):
        # set up the figure
        fig, axs = plt.subplots(1, 1)
        fig.set_figwidth(9)
        fig.set_figheight(6)

        subset = hplc_df.loc[(hplc_df['Organic Acid'] == oa_tested[i])]
        detection_method = pd.unique(subset['Detection Method'])

        if len(detection_method) == 1:
            multiplier = 1
            width = 0.30
        elif len(detection_method) == 2:
            multiplier = 0.5
            width = 0.25
        elif len(detection_method) == 3:
            multiplier = 0
            width = 0.2
        elif len(detection_method) == 4:
            multiplier = -0.5
            width = 0.15
        else:
            multiplier = -1
            width = 0.1
        for j in range(len(detection_method)):
            subset_of_subset = subset.loc[(subset['Detection Method'] == detection_method[j])]

            subset_mean = subset_of_subset.groupby(['Experimental Conditions'], as_index=False).mean(numeric_only=True)
            subset_mean = add_condition_for_sorting(hplc_df, subset_mean)
            subset_mean.sort_values(by=['Condition'], inplace=True)

            subset_std = subset_of_subset.groupby(['Experimental Conditions'], as_index=False).std(numeric_only=True)
            subset_std = add_condition_for_sorting(hplc_df, subset_std)
            subset_std.sort_values(by=['Condition'], inplace=True)

            xx = np.arange(len(subset_mean['Experimental Conditions']))
            subset_mean.insert(0, 'x-axis', xx)

            offset = width * multiplier

            if subset_mean.empty == False:
                if pd.isna(subset_std['Final - Initial Concentration (g/L)'].iloc[0]) == False:
                    axs.bar(subset_mean['x-axis'] + offset, subset_mean['Final - Initial Concentration (g/L)'], width,
                            yerr=subset_std['Final - Initial Concentration (g/L)'],
                            color=color_dictionary[detection_method[j]])
                else:
                    axs.bar(subset_mean['x-axis'] + offset, subset_mean['Final - Initial Concentration (g/L)'], width,
                            color=color_dictionary[detection_method[j]])

            axs.scatter(subset_mean['x-axis'] + offset, subset_mean['Initial Concentration (g/L)'], marker='_', color='red', s=64)

            multiplier += 1

        axs.set_title(oa_tested[i])
        axs.tick_params(axis='x', labelrotation=90)
        axs.set_xticks(subset_mean['x-axis'] + width, subset_mean['Experimental Conditions'])
        axs.axhline(0, color='black', linestyle='dotted')

        fig.supylabel('Concentration (g/L): Final - Initial')
        plt.subplots_adjust(top=0.9, bottom=0.45, left=0.2, right=0.8)

        patchList = []
        for key in color_dictionary:
            data_key = mpatches.Patch(color=color_dictionary[key], label=key)
            patchList.append(data_key)
        axs.legend(handles=patchList, bbox_to_anchor=(1, 1), loc='upper left')

        plt.show()
        # plt.savefig('/Users/jennifer/Desktop/' + str(i))

def save(hplc_df, experimentID):
    hplc_df.to_csv(experimentID + '.csv')

def main():
    # read in inputs
    plot_final_and_initial, plot_final_minus_initial, colors, experimentID, yaml_definitions_path, yaml_layouts_path, yaml_fillstock_path, hplc_data_path, last_characters, num_last_characters, conversion_dict, other_sources, match_chemical_w_ion = inputs()
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
    # get concentrations from fillstock
    hplc_df = get_initial_from_fillstock(hplc_df, experimentID, yaml_definitions_path, yaml_fillstock_path, conversion_dict, match_chemical_w_ion)
    # get concentrations from definitions yaml
    hplc_df = get_conc_from_definitions(hplc_df, experimentID, yaml_definitions_path, match_chemical_w_ion, conversion_dict, yaml_fillstock_path)
    # get concentrations from other sources
    hplc_df = get_cond_from_other_sources(hplc_df, other_sources, conversion_dict)
    # get initial concentrations of organic acids
    hplc_df = get_initial_concentration(hplc_df)
    # make color dictionary
    color_dictionary = make_color_dictionary(colors, hplc_df)
    # Plot stuff
    if plot_final_and_initial == 'Yes':
        # grid plot final and initial
        grid_plot(hplc_df, color_dictionary)
        # individual plots final and initial
        individual_plots(hplc_df, color_dictionary)
    if plot_final_minus_initial == 'Yes':
        # get final minus initial concentration
        hplc_df = get_final_minus_initial(hplc_df)
        # grid plot final minus initial
        grid_plot_final_minus_initial(hplc_df, color_dictionary)
        # individual plots final minus initial
        individual_plots_final_minus_initial(hplc_df, color_dictionary)
    # Save data frame
    save(hplc_df, experimentID)
    #hplc_df.to_csv('/Users/jennifer/Desktop/HPLC_testing/organic_acids/testing.csv')

main()