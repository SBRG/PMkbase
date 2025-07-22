import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from flask import Flask, render_template, url_for, jsonify
from flask import request
import csv
from flask import current_app as app
from models import db, GrowthData, TraitData, KineticData
from sqlalchemy import or_
import os
import json
from flask import Flask, request, render_template, redirect, url_for
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String, Float
from sqlalchemy.orm import sessionmaker
import csv
import openpyxl
from io import StringIO, BytesIO

species = ['ecoli','pputida','saureus']

plate_map = {'PM01':'Carbon Phenotypes','PM02A':'Carbon Phenotypes','PM03B':'Nitrogen Phenotypes','PM04A':'Phosphorus/Sulfur Phenotypes','PM05':'Supplement Phenotypes',
             'PM06':'Nitrogen Phenotypes','PM07':'Nitrogen Phenotypes','PM08':'Nitrogen Phenotypes','PM09':'Stress Phenotypes',
             'PM10':'Stress Phenotypes','PM11C':'Abx Phenotypes','PM12B':'Abx Phenotypes'}


def get_trait_summary():
    """
    Loads or generates trait summary data for all species and plates.

    Returns:
        tuple: (traits, categories)
            traits (list): List of trait dictionaries for each species.
            categories (list): List of phenotype categories.
    """
    json_file_path = os.path.join('static', 'trait_summary.json')
    if os.path.exists(json_file_path):
        with open(json_file_path, 'r') as json_file:
            data = json.load(json_file)
        traits = data['traits']
        categories = data['categories']
    else:
        all_entries = TraitData.query.all()
        entries = [{'Specie': entry.specie, 'Plate': entry.plate} for entry in all_entries]
        entries = pd.DataFrame(entries)
        entries['Description'] = [plate_map[plate] for plate in entries['Plate']]
        categories = entries['Description'].unique().tolist()
        species = entries['Specie'].unique().tolist()
        traits = []
        for specie in species:
            av_traits = []
            for category in categories:
                av_traits.append(entries[(entries['Specie'] == specie) & (entries['Description'] == category)].shape[0])
            traits.append({'name': specie[0].upper() + '. ' + specie[1:], 'data': av_traits})
        data = {'traits': traits, 'categories': categories}
        with open(json_file_path, 'w') as json_file:
            json.dump(data, json_file)
    return traits, categories

def strain_summary():
    """
    Summarizes the number of strains and plates for each species.

    Returns:
        pd.DataFrame: DataFrame with 'Num Strains' and 'Num Plates' for each species.
    """
    strain_summary = pd.DataFrame(index=species)
    total_strains = []
    total_plates = []
    for specie in species:
        comp_strains = []
        temp_summary = pd.read_csv('static/'+specie+'/metadata/updated_summary.csv')
        total_plates.append(temp_summary.shape[0])
        strains = temp_summary['Strain ID']
        mods = temp_summary['Metadata/Modifications']
        for st in range(0,len(strains)):
            comp_strains.append(strains[st]+mods[st])
        total_strains.append(len(list(set(comp_strains))))
    strain_summary['Num Strains'] = total_strains
    strain_summary['Num Plates'] = total_plates
    return strain_summary

def strain_summary_json():
    """
    Converts strain summary to a JSON-like list of dictionaries.

    Returns:
        list: List of dictionaries with species name, number of strains, and number of plates.
    """
    total_strains = strain_summary()
    out2 = []
    for i in total_strains.index:
        num_specie = total_strains.loc[i,'Num Strains']
        num_plates = total_strains.loc[i,'Num Plates']
        specie = i
        temp_dict = {"name":specie,"y":num_specie,"z":num_plates}
        out2.append(temp_dict)
    return out2

def plate_summary():
    """
    Summarizes the number of plates used for each plate type.

    Returns:
        pd.DataFrame: DataFrame indexed by plate type with column 'num_plates'.
    """
    concat_summary = pd.DataFrame()
    for specie in species:
        temp_summary = pd.read_csv('static/'+specie+'/metadata/updated_summary.csv')
        concat_summary = pd.concat([concat_summary,temp_summary])
    plates = concat_summary['Plate'].unique()
    total_plates = []
    for plate in plates:
        total_plates.append(concat_summary.loc[concat_summary['Plate']==plate].shape[0])
    total_plates_used = pd.DataFrame(index=plates)
    total_plates_used['num_plates'] = total_plates
    return total_plates_used

def get_strain_data(plateid,specie):
    """
    Retrieves growth data and compound mapping for a given plate ID and species.

    Args:
        plateid (str): Plate ID.
        specie (str): Species name.

    Returns:
        tuple: (growth_data, well_char, well_num, compound_dict)
    """
    well_char = ['A','B','C','D','E','F','G','H']
    well_num = ['01','02','03','04','05','06','07','08','09','10','11','12']
    growth_frame = pd.read_csv('static/'+specie+'/data/growth_summary.csv',index_col='Plate IDs')
    growth_frame = growth_frame.loc[plateid]
    growth_calls = np.array(growth_frame['Growth'])
    growth_data = []
    growth_calls = growth_calls.reshape((8,12))
    compounds = growth_frame[['Well','Compound']]
    compound_dict = {}
    for i in range (0,compounds.shape[0]):
        compound_dict[compounds.iloc[i,0]] = compounds.iloc[i,1]
    for i in range(0,12):
        for j in range(0,8):
            growth_data.append([i,j,growth_calls[j,i]])
    return growth_data,well_char,well_num,compound_dict

def get_kinetic_parameters(plateid,strain,param='Max Resp'):
    """
    Retrieves kinetic parameters for a given plate ID and strain.

    Args:
        plateid (str): Plate ID.
        strain (str): Strain name.
        param (str): Parameter to extract (default 'Max Resp').

    Returns:
        tuple: (categories, main_data, error_data, param)
    """
    growth_frame = pd.read_csv('static/'+strain+'/data/kinetic_summary.csv')
    growth_calls = pd.read_csv('static/'+strain+'/data/growth_summary.csv')
    growth_frame = growth_frame.loc[growth_frame['Plate IDs']==plateid]
    growth_calls = growth_calls.loc[growth_calls['Plate IDs']==plateid]
    categories = [w['Well']+': '+w['Compound'] for _,w in growth_frame[['Well','Compound']].drop_duplicates().iterrows()]
    numeric_cols = ['Well','Max Resp', 'Max Resp Rate', 'Time till max resp rate', 'AUC']
    params = growth_frame[numeric_cols].set_index('Well')
    avg_df = params.groupby('Well').mean().reset_index()
    max_df = params.groupby('Well').max().reset_index()
    min_df = params.groupby('Well').min().reset_index()
    colors = ["#FFD1DC" if w == 1 else "#C0C0FF" if w == 0.5 else "#808080" for w in growth_calls['Growth']]
    main_data = []
    error_data = []
    for i in range(0,avg_df.shape[0]):
        main_data.append({'y':avg_df[param][i],'color':colors[i]})
    for i in range(0,max_df.shape[0]):
        error_data.append([round(max_df[param][i],2),round(min_df[param][i],2)])
    return categories,main_data,error_data,param

def get_growth_table(plateid,strain):
    """
    Returns a table of growth data for a given plate ID and strain.

    Args:
        plateid (str): Plate ID.
        strain (str): Strain name.

    Returns:
        list: List of lists containing well, compound, growth, description, KEGG ID, CAS ID.
    """
    growth_frame = pd.read_csv('static/'+strain+'/data/growth_summary.csv')
    growth_frame = growth_frame.loc[growth_frame['Plate IDs']==plateid]
    out2 = []
    for i in growth_frame.index:
        well = growth_frame.loc[i,'Well']
        compound = growth_frame.loc[i,'Compound']
        growth = growth_frame.loc[i,'Growth']
        if(growth==1):
            growth='Yes'
        elif(growth==0):
            growth='No'
        else:
            growth='Uncertain'
        description = growth_frame.loc[i,'Description']
        kegg = growth_frame.loc[i,'KEGG ID']
        cas = growth_frame.loc[i,'CAS ID']
        out2.append([
            str(well),
            str(compound),
            str(growth),
            str(description),
            "<a href=https://www.genome.jp/entry/"+str(kegg)+">"+str(kegg)+"</a>",
            str(cas)])
    return out2

def get_growth_curves(plateid,specie,well='A01'):
    """
    Retrieves growth curves for a given plate ID, species, and well.

    Args:
        plateid (str): Plate ID.
        specie (str): Species name.
        well (str): Well ID (default 'A01').

    Returns:
        tuple: (growth_data, time_scale)
    """
    growth_curves = pd.read_csv('static/'+specie+'/data/plate_summary.csv')
    growth_curves = growth_curves.loc[growth_curves['Plate IDs']==plateid]
    main_growth_curves = growth_curves.loc[growth_curves['Well']==well]
    compound = main_growth_curves['Compound'].tolist()[0]
    plate = main_growth_curves['Plate'].tolist()[0]
    growth_data = []
    if('PM11' in plate or 'PM12' in plate):
        well_num = int(well[1:])
        well_num = well_num - ((well_num%4)-1)
        if(well_num<10):
            well_char = '0'+str(well_num)
        else:
            well_char = str(well_num)
        control_well = well[0]+well_char
        control_growth_curves = growth_curves.loc[growth_curves['Well']==control_well]
        control_compound = control_growth_curves['Compound'].tolist()[0]
        for i in range(0,main_growth_curves.shape[0]):
            temp_dict = {'name':compound+' R'+str(i+1),'data':main_growth_curves.iloc[i,11:-2].tolist()}
            growth_data.append(temp_dict)
        if(well!=control_well):
            for i in range(0,control_growth_curves.shape[0]):
                temp_dict = {'name':control_compound+' R'+str(i+1),'data':control_growth_curves.iloc[i,11:-2].tolist()}
                growth_data.append(temp_dict)
    elif('PM01' in plate or 'PM02' in plate or 'PM03' in plate or 'PM04' in plate or 'PM05' in plate or 'PM06' in plate or 'PM07' in plate or 'PM08' in plate):
        control_well = 'A01'
        control_growth_curves = growth_curves.loc[growth_curves['Well']==control_well]
        control_compound = control_growth_curves['Compound'].tolist()[0]
        for i in range(0,main_growth_curves.shape[0]):
            temp_dict = {'name':compound+' R'+str(i+1),'data':main_growth_curves.iloc[i,11:-2].tolist()}
            growth_data.append(temp_dict)
        if(well!='A01'):
            for i in range(0,control_growth_curves.shape[0]):
                temp_dict = {'name':control_compound+' R'+str(i+1),'data':control_growth_curves.iloc[i,11:-2].tolist()}
                growth_data.append(temp_dict)
    elif('PM09' in plate or 'PM10' in plate):
        for i in range(0,main_growth_curves.shape[0]):
            temp_dict = {'name':compound+' R'+str(i+1),'data':main_growth_curves.iloc[i,11:-2].tolist()}
            growth_data.append(temp_dict)
    time_scale = list(np.arange(0,48.25,0.25))
    chart_data = {'categories':time_scale,'data':growth_data}
    return growth_data,time_scale

def get_all_growth_curves(plateid,specie,well='A01'):
    """
    Retrieves all growth curves for a given plate ID, species, and well.

    Args:
        plateid (str): Plate ID.
        specie (str): Species name.
        well (str): Well ID (default 'A01').

    Returns:
        tuple: (growth_data, time_scale)
    """
    time_scale = list(np.arange(0,48.25,0.25))
    signal_columns = [str(w)+'hrs' for w in time_scale]
    growth_data = []
    growth_curves = pd.read_csv('static/'+specie+'/data/plate_summary.csv')
    growth_curves = growth_curves.loc[growth_curves['Plate IDs']==plateid]
    plate = growth_curves['Plate'].tolist()[0]
    if('PM01' in plate or 'PM02' in plate or 'PM03' in plate or 'PM04' in plate or 'PM05' in plate or 'PM06' in plate or 'PM07' in plate or 'PM08' in plate):
        control_well = 'A01'
        control_growth_curves = growth_curves.loc[growth_curves['Well']==control_well]
        main_growth_curves = growth_curves.loc[growth_curves['Well']==well]
        for _,row in main_growth_curves.iterrows():
            temp_dict = {'name':row['Compound']+' '+row['Replicates'],'data':row[signal_columns].tolist()}
            growth_data.append(temp_dict)
        if(well!='A01'):
            for _,row in control_growth_curves.iterrows():
                temp_dict = {'name':row['Compound']+' '+row['Replicates'],'data':row[signal_columns].tolist()}
                growth_data.append(temp_dict)
    elif('PM09' in plate or 'PM10' in plate):
        main_growth_curves = growth_curves.loc[growth_curves['Well']==well]
        for _,row in main_growth_curves.iterrows():
            temp_dict = {'name':row['Compound']+' '+row['Replicates'],'data':row[signal_columns].tolist()}
            growth_data.append(temp_dict)
    elif('PM11' in plate or 'PM12' in plate):
        control_well = well[0]+'01'
        control_growth_curves = growth_curves.loc[growth_curves['Well']==control_well]
        main_growth_curves = growth_curves.loc[growth_curves['Well']==well]
        for _,row in main_growth_curves.iterrows():
            temp_dict = {'name':row['Compound']+' '+row['Replicates'],'data':row[signal_columns].tolist()}
            growth_data.append(temp_dict)
        if(well!='A01'):
            for _,row in control_growth_curves.iterrows():
                temp_dict = {'name':row['Compound']+' '+row['Replicates'],'data':row[signal_columns].tolist()}
                growth_data.append(temp_dict)
    return growth_data,time_scale

def get_compound_drop_down(plate):
    """
    Returns a list of compounds for a given plate for dropdown selection.

    Args:
        plate (str): Plate ID.

    Returns:
        list: List of strings in the format '<Well>: <Compound>'.
    """
    platedesc = pd.read_csv('static/plate_desc/platedesc.csv')
    platedesc = platedesc[platedesc['Plate']==plate]
    dropdown_compounds = []
    for _,row in platedesc.iterrows():
        dropdown_compounds.append(row['Well']+': '+row['Compound'])
    return dropdown_compounds

def get_control_well_distribution(specie):
    """
    Returns kinetic parameters for control wells for a given species.

    Args:
        specie (str): Species name.

    Returns:
        tuple: Lists of kinetic parameters for growth, no growth, and uncertain growth.
    """
    kinetic_data = pd.read_csv('static/'+specie+'/data/kinetic_summary.csv',index_col='Plate IDs')
    growth = kinetic_data.loc[kinetic_data['Growth']==1]
    no_growth = kinetic_data.loc[kinetic_data['Growth']==0]
    uncertain_growth = kinetic_data.loc[kinetic_data['Growth']==0.5]
    growth_max_resp = growth['Max Resp'].tolist()
    growth_max_resp_rate = growth['Max Resp Rate'].tolist()
    growth_max_time = growth['Time till max resp rate'].tolist()
    growth_max_auc = growth['AUC'].tolist()
    no_growth_max_resp = no_growth['Max Resp'].tolist()
    no_growth_max_resp_rate = no_growth['Max Resp Rate'].tolist()
    no_growth_max_time = no_growth['Time till max resp rate'].tolist()
    no_growth_max_auc = no_growth['AUC'].tolist()
    uncertain_growth_max_resp = uncertain_growth['Max Resp'].tolist()
    return growth_max_resp,growth_max_resp_rate,growth_max_time,growth_max_auc,no_growth_max_resp,no_growth_max_resp_rate,no_growth_max_time,no_growth_max_auc,uncertain_growth_max_resp

def get_control_well_dist(specie):
    """
    Returns control well and growth well data for a given species.

    Args:
        specie (str): Species name.

    Returns:
        tuple: (control_data, growth_data)
    """
    kinetic_data = pd.read_csv('static/'+specie+'/data/kinetic_summary.csv',index_col='Plate IDs')
    control_wells = kinetic_data.loc[kinetic_data['Well']=='A01']
    growth_wells = kinetic_data.loc[kinetic_data['Growth']==1]
    control_data = []
    growth_data = []
    for well in range(0,control_wells.shape[0]):
        plate = control_wells.iloc[well,0]
        if('PM01' in plate or 'PM02' in plate or 'PM03' in plate or 'PM04' in plate or 'PM05' in plate or 'PM06' in plate or 'PM07' in plate or 'PM08' in plate):
            control_data.append(control_wells.iloc[well,5])
            growth_data.append(growth_wells.iloc[well,5])
    return control_data,growth_data

def load_cluster_data(specie):
    """
    Loads cluster data for a given species.

    Args:
        specie (str): Species name.

    Returns:
        dict: Mapping of strain to phylogroup/genome cluster.
    """
    cluster_data = {}
    with open('static/'+specie+'/metadata/tree_clusters.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            cluster_data[row['Strain']] = row['Phylogroup/Genome Cluster']
    return cluster_data

def load_specie_metadata(specie):
    """
    Loads metadata summary for a given species.

    Args:
        specie (str): Species name.

    Returns:
        tuple: (samples, strains, available_plates, clusters, plates)
    """
    specie_summary = pd.read_csv('static/'+specie+'/metadata/updated_summary.csv',index_col='Plate IDs')
    samples = specie_summary.shape[0]
    strains = len(specie_summary['Strain'].unique())
    plates = specie_summary['Plate'].unique()
    available_plates = ''
    for plate in plates[0:-1]:
        available_plates = available_plates+plate+', '
    available_plates = available_plates+plates[-1]
    clusters = len(specie_summary['Phylogroup/Genome Cluster'].unique())
    return samples,strains,available_plates,clusters,plates

def get_compounds_from_plates(plates):
    """
    Returns a list of compound descriptions for the given plate IDs.

    Args:
        plates (list): List of plate IDs.

    Returns:
        list: List of strings describing compounds in the format:
              'Plate: <Plate>, Well: <Well>, Compound: <Compound>, Desc: <Description>'
    """
    platedesc = pd.read_csv('static/plate_desc/platedesc.csv')
    plate_details = platedesc[platedesc['Plate'].isin(plates)]
    plate_details = plate_details[['Plate','Well','Compound','Description']]
    available_comps = ['Plate: '+w['Plate']+', Well: '+w['Well']+', Compound: '+w['Compound']+', Desc: '+w['Description'] 
                       for _,w in plate_details.iterrows()]
    return available_comps

def combine_specie_summaries():
    """
    Combines strain summaries from all species into a single list of dictionaries.

    Returns:
        list: List of dictionaries with keys 'Strain ID', 'Modification', and 'Specie'.
    """
    summary= pd.DataFrame()
    for specie in species:
        temp_dataframe = pd.DataFrame()
        strains = []
        strain_id = []
        mods = []
        sps = []
        temp_summary = pd.read_csv('static/'+specie+'/metadata/updated_summary.csv',index_col='Plate IDs')
        for i in range(0,temp_summary.shape[0]):
            strains.append(temp_summary.iloc[i,1]+'___'+temp_summary.iloc[i,2])
        strains = list(set(strains))
        for strain in strains:
            strain_id.append(strain.split('___')[0])
            mods.append(strain.split('___')[1])
            sps.append(temp_summary['Specie'].tolist().pop(0))
        temp_dataframe['Strain ID'] = strain_id
        temp_dataframe['Modification'] = mods
        temp_dataframe['Specie'] = sps
        summary = pd.concat([summary,temp_dataframe])
    return summary.to_dict('records')

def get_all_compounds_in_all_wells():
    """
    Retrieves all unique compounds and their descriptions from all wells.

    Returns:
        list: List of strings in the format '<Compound>, <Description>'.
    """
    platedesc = pd.read_csv('static/plate_desc/platedesc.csv')
    compounds = []
    for i in range(0,platedesc.shape[0]):
        compounds.append(platedesc.iloc[i,2]+', '+platedesc.iloc[i,3])
    compounds = list(set(compounds))
    return compounds

def get_plate_well_from_compound(compound):
    """
    Finds the plate and well for a given compound description.

    Args:
        compound (str): Compound description in the format '<Compound>, <Description>'.

    Returns:
        tuple: (plate, well) corresponding to the compound.
    """
    platedesc = pd.read_csv('static/plate_desc/platedesc.csv')
    for i in range(platedesc.shape[0]):
        combined_comp = platedesc['Compound'][i]+', '+platedesc['Description'][i]
        if(compound==combined_comp):
            plate = platedesc['Plate'][i]
            well = platedesc['Well'][i]
    return plate,well

def get_plateid_from_strain(strain_list,plate):
    """
    Retrieves plate IDs for a list of strains and a given plate.

    Args:
        strain_list (list): List containing strain information in groups of three.
        plate (str): Plate ID to search for.

    Returns:
        list: List of plate IDs or 'N.A' if not found.
    """
    combined_summary = pd.DataFrame()
    plateids = []
    for specie in species:
        temp_summary = pd.read_csv('static/'+specie+'/metadata/updated_summary.csv',index_col='Plate IDs')
        combined_summary = pd.concat([combined_summary,temp_summary])
    for i in range(0,len(strain_list),3):
        strain = combined_summary[combined_summary['Strain']==strain_list[i]]
        metadata = strain[strain['Modification/Metadata']==strain_list[i+1]]
        plates = metadata['Plate'].tolist()
        if(plate in plates):
            plateids.append((metadata[metadata['Plate']==plate]).index.tolist().pop(0))
        else:
            plateids.append('N.A')
    return plateids

def get_growth_calls_from_plateids(plateids,well,xlabels):
    """
    Retrieves growth calls and signal series for a list of plate IDs and a well.

    Args:
        plateids (list): List of plate IDs.
        well (str): Well ID.
        xlabels (list): List of labels for series.

    Returns:
        tuple: (growth_calls, series, time)
    """
    combined_growth = pd.DataFrame()
    combined_signals = pd.DataFrame()
    growth_calls = []
    series = []
    for specie in species:
        temp_growth = pd.read_csv('static/'+specie+'/data/growth_summary.csv',index_col='Plate IDs')
        temp_signal = pd.read_csv('static/'+specie+'/data/plate_summary.csv',index_col='Plate IDs')
        combined_growth = pd.concat([combined_growth,temp_growth])
        combined_signals = pd.concat([combined_signals,temp_signal])
    i = 0
    for id in plateids:
        if(id=='N.A'):
            growth_calls.append([i,0,0.75])
        else:
            growth = combined_growth.loc[id]
            growth_calls.append([i,0,growth[growth['Well']==well]['Growth'].tolist().pop()])
        i = i+1
    j = 0
    for id in plateids:
        if(id=='N.A'):
            j = j+1
            continue
        else:
            signals = combined_signals.loc[id]
            signals = signals[signals['Well']==well]
            for i in range(0,signals.shape[0]):
                if(i>=2):
                    break
                series.append({'name':xlabels[j]+' R'+str(i+1),'data':signals.iloc[i,7:].tolist()})
        j = j+1
    time = list(np.linspace(0,48,193))
    return growth_calls,series,time

def get_strain_names(strainlist):
    """
    Returns formatted strain names from a list.

    Args:
        strainlist (list): List of strain information.

    Returns:
        list: List of formatted strain names.
    """
    strain_names = []
    for i in range(0,len(strainlist),3):
        strain_names.append(strainlist[i]+'__'+strainlist[i+1]+'__'+strainlist[i+2])
    return strain_names

def get_tracking_growth_data(specie,plate,well):
    """
    Retrieves growth status for all strains for a given species, plate, and well.

    Args:
        specie (str): Species name.
        plate (str): Plate ID.
        well (str): Well ID.

    Returns:
        tuple: (cluster_data, growth_strains, nogrowth_strains)
    """
    all_strains = pd.read_csv('static/'+specie+'/metadata/updated_summary.csv',index_col='Plate IDs')
    all_strains = all_strains['Strain'].unique().tolist()
    growth = TraitData.query.filter(
        TraitData.specie == specie,
        TraitData.growth == 1,
        TraitData.plate == plate,
        TraitData.well == well,
        TraitData.metadata_mods.like('%WT%')
    ).all()
    uncertain = TraitData.query.filter(
        TraitData.specie == specie,
        TraitData.growth == 0.5,
        TraitData.plate == plate,
        TraitData.well == well,
        TraitData.metadata_mods.like('%WT%')
    ).all()
    no_growth = TraitData.query.filter(
        TraitData.specie == specie,
        TraitData.growth == 0,
        TraitData.plate == plate,
        TraitData.well == well,
        TraitData.metadata_mods.like('%WT%')
    ).all()
    growth_strains = [entry.strain for entry in growth]
    uncertain_strains = [entry.strain for entry in uncertain]
    nogrowth_strains = [entry.strain for entry in no_growth]
    cluster_data = {}
    for strain in growth_strains:
        cluster_data[strain] = 'yes'
    for strain in uncertain_strains:
        cluster_data[strain] = 'uncertain/maybe'
    for strain in nogrowth_strains:
        cluster_data[strain] = 'no'
    for strain in all_strains:
        if(strain in growth_strains):
            cluster_data[strain] = 'yes'
        elif(strain in uncertain_strains):
            cluster_data[strain] = 'uncertain/maybe'
        elif(strain in nogrowth_strains):
            cluster_data[strain] = 'no'
        else:
            cluster_data[strain] = 'data N.A'
    return cluster_data,growth_strains,nogrowth_strains

def calculate_specie_inter_cluster_mash_dist(specie):
    """
    Calculates intra-cluster mash distances for a given species.

    Args:
        specie (str): Species name.

    Returns:
        list: List of mash distances within clusters.
    """
    summary = pd.read_csv('static/'+specie+'/metadata/updated_summary.csv',index_col='Plate IDs')
    df_mash = pd.read_csv('static/'+specie+'/mash_distances.tsv',sep='\t')
    clusters = summary[['Strain','Phylogroup/Genome Cluster']].drop_duplicates(keep='first')
    cluster_distance = []
    for cluster in clusters['Phylogroup/Genome Cluster'].unique():
        clst_strains = clusters[clusters['Phylogroup/Genome Cluster']==cluster]['Strain'].unique().tolist()
        cluster_distance = cluster_distance+df_mash[(df_mash['genome1'].isin(clst_strains))&(df_mash['genome2'].isin(clst_strains))]['mash_distance'].unique().tolist()
    return cluster_distance

def calculate_phenotype_median_mash(specie,strains):
    """
    Calculates mash distances for a set of strains in a species.

    Args:
        specie (str): Species name.
        strains (list): List of strain names.

    Returns:
        list: List of mash distances.
    """
    df_mash = pd.read_csv('static/'+specie+'/mash_distances.tsv',sep='\t')
    return df_mash[(df_mash['genome1'].isin(strains))&(df_mash['genome2'].isin(strains))]['mash_distance'].unique().tolist()

def get_tracking_kinetic_params(no_growth_strains,specie,plate,well,growth_strains,param='Max Resp'):
    """
    Retrieves kinetic parameters for growth and no-growth strains.

    Args:
        no_growth_strains (list): List of no-growth strain names.
        specie (str): Species name.
        plate (str): Plate ID.
        well (str): Well ID.
        growth_strains (list): List of growth strain names.
        param (str): Parameter to extract (default 'Max Resp').

    Returns:
        tuple: (kinetic_means, kinetic_errors, categories)
    """
    kinetic_means = []
    kinetic_errors = []
    categories = []
    growth_strain_filter = or_(*[KineticData.strain == strain for strain in growth_strains])
    no_growth_strain_filter = or_(*[KineticData.strain == strain for strain in no_growth_strains])
    growth_kinetics = KineticData.query.filter(
        growth_strain_filter,
        KineticData.specie == specie,
        KineticData.plate == plate,
        KineticData.well == well,
        KineticData.metadata_mods.like('%WT%')
    ).all()
    no_growth_kinetics = KineticData.query.filter(
        no_growth_strain_filter,
        KineticData.specie == specie,
        KineticData.plate == plate,
        KineticData.well == well,
        KineticData.metadata_mods.like('%WT%')
    ).all()
    for strain in growth_strains:
        if(param=='Max Resp'):
            param_list = [w.maxresp for w in growth_kinetics if w.strainid == strain]
            if(param_list):
                kinetic_means.append({'y':round(np.mean(param_list),2),'color':'green'})
                kinetic_errors.append([round(np.max(param_list),2),round(np.min(param_list),2)])
                categories.append(strain)
            else:
                continue
        elif(param=='Max Resp Rate'):
            param_list = [w.maxresprate for w in growth_kinetics if w.strainid == strain]
            if(param_list):
                kinetic_means.append({'y':round(np.mean(param_list),2),'color':'green'})
                kinetic_errors.append([round(np.max(param_list),2),round(np.min(param_list),2)])
                categories.append(strain)
            else:
                continue
        elif(param=='Time till max resp'):
            param_list = [w.timetill for w in growth_kinetics if w.strainid == strain]
            if(param_list):
                kinetic_means.append({'y':round(np.mean(param_list),2),'color':'green'})
                kinetic_errors.append([round(np.max(param_list),2),round(np.min(param_list),2)])
                categories.append(strain)
            else:
                continue
        elif(param=='AUC'):
            param_list = [w.auc for w in growth_kinetics if w.strainid == strain]
            if(param_list):
                kinetic_means.append({'y':round(np.mean(param_list),2),'color':'green'})
                kinetic_errors.append([round(np.max(param_list),2),round(np.min(param_list),2)])
                categories.append(strain)
            else:
                continue
    for strain in no_growth_strains:
        if(param=='Max Resp'):
            param_list = [w.maxresp for w in no_growth_kinetics if w.strainid == strain]
            if(param_list):
                kinetic_means.append({'y':round(np.mean(param_list),2),'color':'red'})
                kinetic_errors.append([round(np.max(param_list),2),round(np.min(param_list),2)])
                categories.append(strain)
            else:
                continue
        elif(param=='Max Resp Rate'):
            param_list = [w.maxresprate for w in no_growth_kinetics if w.strainid == strain]
            if(param_list):
                kinetic_means.append({'y':round(np.mean(param_list),2),'color':'red'})
                kinetic_errors.append([round(np.max(param_list),2),round(np.min(param_list),2)])
                categories.append(strain)
            else:
                continue
        elif(param=='Time till max resp'):
            param_list = [w.timetill for w in no_growth_kinetics if w.strainid == strain]
            if(param_list):
                kinetic_means.append({'y':round(np.mean(param_list),2),'color':'red'})
                kinetic_errors.append([round(np.max(param_list),2),round(np.min(param_list),2)])
                categories.append(strain)
            else:
                continue
        elif(param=='AUC'):
            param_list = [w.auc for w in no_growth_kinetics if w.strainid == strain]
            if(param_list):
                kinetic_means.append({'y':round(np.mean(param_list),2),'color':'red'})
                kinetic_errors.append([round(np.max(param_list),2),round(np.min(param_list),2)])
                categories.append(strain)
            else:
                continue
    return kinetic_means,kinetic_errors,categories

def get_mash_strain_names(specie):
    """
    Returns all unique genome names for a given species from mash distances file.

    Args:
        specie (str): Species name.

    Returns:
        list: List of genome names.
    """
    mash_dist = pd.read_csv('static/'+specie+'/mash_distances.tsv',sep='\t')
    gen1 = mash_dist['genome1'].tolist()
    gen2 = mash_dist['genome2'].tolist()
    genomes = list(set(gen1+gen2))
    genomes = [w.replace('_','') for w in genomes]
    return genomes

def get_example_upload_data():
    """
    Loads example upload data for demonstration.

    Returns:
        list: List of lists representing example rows.
    """
    example_data = pd.read_csv('static/upload_data/example.csv')
    out2 = []
    for i in example_data.index[0:8]:
        out2.append([
            str(example_data.loc[i,'Plate Type']),
            str(example_data.loc[i,'Media']),
            str(example_data.loc[i,'Strain']),
            str(example_data.loc[i,'Specie']),
            str(example_data.loc[i,'Hr']),
            str(round(example_data.loc[i,'A01'],2)),
            str(round(example_data.loc[i,'A02'],2)),
            '...',
            str(round(example_data.loc[i,'H12'],2))])
    return out2

def standardize_columns(df):
    """
    Standardizes column names in uploaded data to expected names.

    Args:
        df (pd.DataFrame): DataFrame with uploaded data.

    Returns:
        pd.DataFrame: DataFrame with standardized column names.
    """
    new_columns = {}
    column_mapping = {
        'Plate Type': ['Plate Type', 'Plate', 'PlateType'],
        'Media': ['Media', 'Medium'],
        'Strain': ['Strain', 'Strain Name', 'StrainName','Name'],
        'Specie': ['Species','Specie','Organism'],
        'Hr': ['Hr', 'Hour', 'Hours']
    }
    for standard_name, variations in column_mapping.items():
        for col in df.columns:
            if col in variations:
                new_columns[col] = standard_name
                break
    df.rename(columns=new_columns, inplace=True)
    return df


