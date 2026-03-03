from flask import Flask, render_template, url_for, jsonify, send_file, Response, session
from flask_cors import CORS, cross_origin
from flask import request
import pandas as pd
import scripts
import json
import random
from models import db, GrowthData, TraitData,KineticData
import numpy as np
from sqlalchemy import asc
import os
from flask import Flask, request, render_template, redirect, url_for
import csv
import openpyxl
from io import StringIO, BytesIO
from werkzeug.utils import secure_filename
import utils
import uuid
from datetime import timedelta
import shutil
from urllib.parse import urlencode
import zipfile
import io
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__,template_folder='templates',static_url_path='/static')
CORS(app, resources={r"/*": {"origins": "*"}})
app.config['CORS_HEADERS'] = 'Content-Type'

app.secret_key = 'sbrg_omnilog'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///growth_data.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['MAX_CONTENT_LENGTH'] = 10000*1024*1024
db.init_app(app)

# with app.app_context():
#     db.create_all()


def ingest_data(specie):
    """
    Ingests growth curve data for a given species from a CSV file and adds it to the database.

    Args:
        specie (str): The species name.

    Returns:
        None

    Side Effects:
        Reads 'static/{specie}/data/plate_summary.csv'.
        Adds GrowthData entries to the database.
    """
    csv_path = f'static/{specie}/data/plate_summary.csv'
    growth_curves = pd.read_csv(csv_path)
    time_scale = list(np.arange(0,48.25,0.25))
    
    for _, row in growth_curves.iterrows():
        signal_columns = [str(w)+'hrs' for w in time_scale]
        signal_data = row[signal_columns].tolist()
        entry = GrowthData(
            plateid=row['Plate IDs'],
            specie=specie,
            well=row['Well'],
            #plate = row['Plate'],
            compound=row['Compound'],
            replicates=row['Replicates'],
            signal_data=signal_data
        )
        db.session.add(entry)
    db.session.commit()


def ingest_trait_data(specie):
    """
    Ingests trait data for a given species from a CSV file and adds it to the database.

    Args:
        specie (str): The species name.

    Returns:
        None

    Side Effects:
        Reads 'static/{specie}/data/growth_summary.csv'.
        Adds TraitData entries to the database.
    """
    csv_path = f'static/{specie}/data/growth_summary.csv'
    growth_calls = pd.read_csv(csv_path)
    
    for _, row in growth_calls.iterrows():
        entry = TraitData(
            plateid=row['Plate IDs'],
            strainid = row['Strain ID'],
            specie=specie,
            metadata_mods = row['Metadata/Modifications'],
            project = row['Project'],
            well=row['Well'],
            plate = row['Plate'],
            media = row['Media'],
            growth = row['Growth'],
            compound=row['Compound'],
            desc = row['Description'],
            strain = row['Strain'],
            phylo = row['Phylogroup/Genome Cluster'],
            mlst = row['MLST']
            )
        
        db.session.add(entry)
    db.session.commit()

def ingest_kinetic_data(specie):
    """
    Ingests kinetic data for a given species from a CSV file and adds it to the database.

    Args:
        specie (str): The species name.

    Returns:
        None

    Side Effects:
        Reads 'static/{specie}/data/kinetic_summary.csv'.
        Adds KineticData entries to the database.
    """
    csv_path = f'static/{specie}/data/kinetic_summary.csv'
    kinetics = pd.read_csv(csv_path)
    for _, row in kinetics.iterrows():
        entry = KineticData(
            plateid=row['Plate IDs'],
            strainid = row['Strain ID'],
            strain = row['Strain'],
            specie=specie,
            metadata_mods = row['Metadata/Modifications'],
            project = row['Project'],
            well=row['Well'],
            plate = row['Plate'],
            media = row['Media'],
            replicates=row['Replicates'],
            compound=row['Compound'],
            keggid=row['KEGG ID'],
            casid=row['CAS ID'],
            maxresp = row['Max Resp'],
            maxresprate = row['Max Resp Rate'],
            timetill = row['Time till max resp rate'],
            auc = row['AUC'],
            growth = row['Growth'],
            mlst = row['MLST'],
            phylo = row['Phylogroup/Genome Cluster']
            )
        
        db.session.add(entry)
    db.session.commit()

@app.route('/download/all_sequences')
def download_all_sequences():
    """
    Downloads a zip file containing all sequence files.

    Returns:
        Response: Sends 'static/all_sequences.zip' as an attachment.
    """
    zip_path = os.path.join('static', 'all_sequences.zip')
    
    return send_file(zip_path, as_attachment=True, download_name='all_sequences.zip')



@app.route('/download/all_pmdata')
def download_all_pmdata():
    """
    Downloads a zip file containing all phenotype microarray data.

    Returns:
        Response: Sends 'static/all_PM_data.zip' as an attachment.
    """
    zip_path = os.path.join('static', 'all_PM_data.zip')
    
    return send_file(zip_path, as_attachment=True, download_name='allPMdata.zip')


@app.route('/download/all_specie_sequences')
def download_all_specie_sequences():
    """
    Downloads a zip file containing all sequence files for a specific species.

    Query Parameters:
        specie (str): The species name.

    Returns:
        Response: Sends '{specie}/sequences.zip' as an attachment.
    """
    specie=request.args.get('specie')
    
    # Define the path for the zip file
    zip_path = os.path.join('static', specie, 'sequences.zip')
    
    return send_file(zip_path, as_attachment=True, download_name=specie+'_sequences.zip')



@app.route('/download/all_specie_pmdata')
def download_all_specie_pmdata():
    """
    Downloads a zip file containing all PM data for a specific species.

    Query Parameters:
        specie (str): The species name.

    Returns:
        Response: Sends '{specie}/data.zip' as an attachment.
    """
    specie=request.args.get('specie')
    
    # Define the path for the zip file
    zip_path = os.path.join('static', specie, 'data.zip')
    
    return send_file(zip_path, as_attachment=True, download_name=specie+'_allPMdata.zip')



@app.route('/download/mainstrain_specie_sequences')
def download_mainstrain_specie_sequences():
    """
    Downloads the sequence file for a specific strain of a given species.

    Query Parameters:
        specie (str): The species name.
        strain (str): The strain name.

    Returns:
        Response: Sends the sequence file (.fna) as an attachment.
        404: If the file is not found.

    Side Effects:
        Reads files from the static directory.

    Raises:
        FileNotFoundError: If the sequence file does not exist.
    """

    specie= request.args.get('specie')
    strain = request.args.get('strain')

    for filename in os.listdir(os.path.join('static', specie, 'sequences')):
        if filename.startswith(strain):
            file_path = os.path.join('static', specie, 'sequences',filename)
            break
    

    return send_file(file_path, as_attachment=True, download_name=specie+'_'+strain+'.fna')



@app.route('/download/mainstrain_specie_growthdata')
def download_mainstrain_specie_growthdata():
    """
    Downloads growth summary data for a specific strain and plate of a species.

    Query Parameters:
        specie (str): The species name.
        plateid (str): The plate ID.
        strain (str): The strain name.

    Returns:
        Response: Sends filtered growth summary CSV as an attachment.
    """

    specie= request.args.get('specie')
    plateid = request.args.get('plateid')
    strain = request.args.get('strain')


    growth = pd.read_csv('static/'+specie+'/data/growth_summary.csv',index_col='Plate IDs')
    growth_filtered = growth.loc[plateid]

    # Convert the filtered dataframe to a CSV
    csv_string = growth_filtered.to_csv()

    # Create a response object and set the appropriate headers
    response = Response(
        csv_string,
        mimetype='text/csv',
        headers={
            'Content-Disposition': f'attachment;filename={specie}_{strain}_growth_data.csv'
        }
    )

    return response




@app.route('/download/mainstrain_specie_kineticdata')
def download_mainstrain_specie_kineticdata():
    """
    Downloads kinetic summary data for a specific strain and plate of a species.

    Query Parameters:
        specie (str): The species name.
        plateid (str): The plate ID.
        strain (str): The strain name.

    Returns:
        Response: Sends filtered kinetic summary CSV as an attachment.
    """

    specie= request.args.get('specie')
    plateid = request.args.get('plateid')
    strain = request.args.get('strain')


    growth = pd.read_csv('static/'+specie+'/data/kinetic_summary.csv',index_col='Plate IDs')
    growth_filtered = growth.loc[plateid]

    # Convert the filtered dataframe to a CSV
    csv_string = growth_filtered.to_csv()

    # Create a response object and set the appropriate headers
    response = Response(
        csv_string,
        mimetype='text/csv',
        headers={
            'Content-Disposition': f'attachment;filename={specie}_{strain}_kinetic_data.csv'
        }
    )

    return response



@app.route('/download/mainstrain_specie_rawdata')
def download_mainstrain_specie_rawdata():
    """
    Downloads raw plate summary data for a specific strain and plate of a species.

    Query Parameters:
        specie (str): The species name.
        plateid (str): The plate ID.
        strain (str): The strain name.

    Returns:
        Response: Sends filtered plate summary CSV as an attachment.
    """

    specie= request.args.get('specie')
    plateid = request.args.get('plateid')
    strain = request.args.get('strain')


    growth = pd.read_csv('static/'+specie+'/data/plate_summary.csv',index_col='Plate IDs')
    growth_filtered = growth.loc[plateid]

    # Convert the filtered dataframe to a CSV
    csv_string = growth_filtered.to_csv()

    # Create a response object and set the appropriate headers
    response = Response(
        csv_string,
        mimetype='text/csv',
        headers={
            'Content-Disposition': f'attachment;filename={specie}_{strain}_raw_data.csv'
        }
    )

    return response



@app.route('/')
@app.route('/index')
def index():
    """
    Renders the main index page with trait summary data.

    Returns:
        Rendered HTML template 'index.html' with trait and category data.
    """
    traits,categories = scripts.get_trait_summary()

    return render_template('index.html',series=traits,categories=categories)

@app.route('/dashboard')
def dashboard():
    """
    Renders the dashboard page.

    Returns:
        Rendered HTML template 'index.html'.
    """
    return render_template('index.html')



@app.route('/tree')
def tree():
    """
    Renders the species tree page with summary statistics.

    Query Parameters:
        specie (str): The species name.

    Returns:
        Rendered HTML template 'tree.html' with tree and compound data.
    """
    specie=request.args.get('specie')
    specie_summary = pd.read_csv('static/'+specie+'/metadata/summary.csv',index_col='Plate IDs')
    plates = specie_summary['Plate'].unique()
    comps = scripts.get_compounds_from_plates(plates)
    spe_mash = scripts.calculate_specie_inter_cluster_mash_dist(specie)
    spe_mash = [w for w in spe_mash if w!=0]
    std_dev = np.var(spe_mash)#np.std(spe_mash, ddof=1)
    specie_inter_mash = np.median(spe_mash)
    lb = specie_inter_mash - std_dev
    ub = specie_inter_mash + std_dev
    specie_mash_min = min(spe_mash, key=lambda x: abs(x - ub))
    specie_mash_max = min(spe_mash, key=lambda x: abs(x - lb))

    return render_template('tree.html',specie=specie,comps=comps,specie_inter_mash=specie_inter_mash,spe_min = specie_mash_min,spe_max = specie_mash_max,
                           specie_individual_points = spe_mash)

@app.route('/tree_json')
def get_tree():
    """
    Returns the species tree data in JSON format, annotated with cluster information.

    Query Parameters:
        specie (str): The species name.

    Returns:
        JSON: Tree data with cluster annotations.
    """
    specie=request.args.get('specie')
    with open("static/"+specie+"/tree.json", "r") as f:
        tree_data = json.load(f)

    cluster_data = scripts.load_cluster_data(specie)

    def add_clusters(node):
        if 'name' in node:
            node['cluster'] = cluster_data.get(node['name'], None)
        if 'children' in node:
            for child in node['children']:
                add_clusters(child)
        return node

    tree_data = add_clusters(tree_data)
    return jsonify(tree_data)



@app.route('/track_tree_json',methods=['GET'])
def get_track_tree():
    """
    Returns tracked tree data and phenotype/kinetic statistics for a given plate and well.

    Query Parameters:
        specie (str): The species name.
        plate (str): The plate name (default 'PM01').
        well (str): The well name (default 'H12').

    Returns:
        JSON: Tree data, phenotype mash statistics, kinetic means/errors, and strain lists.
    """
    from scipy.stats import ttest_ind

    specie=request.args.get('specie')
    plate = request.args.get('plate', 'PM01')
    well = request.args.get('well', 'H12')
    with open("static/"+specie+"/tree.json", "r") as f:
        tree_data = json.load(f)

    cluster_data,growth_strains,nogrowth_strains = scripts.get_tracking_growth_data(specie,plate,well)

    kinetic_means,kinetic_errors,categories = scripts.get_tracking_kinetic_params(nogrowth_strains,specie,plate,well,growth_strains,param='Max Resp')
    spe_mash = scripts.calculate_specie_inter_cluster_mash_dist(specie)
    spe_mash = [w for w in spe_mash if w!=0]

    if(len(growth_strains)>2):
        phen_mash = scripts.calculate_phenotype_median_mash(specie,growth_strains)
        phen_mash = [w for w in phen_mash if w!=0]
        phenotype_mash = np.median(phen_mash)
        std_dev = np.std(phen_mash, ddof=0)#np.var(phen_mash)#np.std(phen_mash, ddof=1)
        lb = phenotype_mash - std_dev
        ub = phenotype_mash + std_dev
        phenotype_mash_min = min(phen_mash, key=lambda x: abs(x - lb))
        phenotype_mash_max = min(phen_mash, key=lambda x: abs(x - ub))
        _, p_value = ttest_ind(spe_mash, phen_mash, equal_var=False)
    
    elif(len(growth_strains)==2):
        phen_mash = scripts.calculate_phenotype_median_mash(specie,growth_strains)
        phen_mash = [w for w in phen_mash if w!=0]
        phenotype_mash = np.median(phen_mash)
        phenotype_mash_min = phenotype_mash
        phenotype_mash_max = phenotype_mash
        p_value = 'null (only 2 strains grow)'#ttest_ind(spe_mash, phen_mash, equal_var=False)

    elif(len(growth_strains)==1):
        phen_mash = scripts.calculate_phenotype_median_mash(specie,growth_strains)
        phenotype_mash = np.median(phen_mash)
        phenotype_mash_min = phenotype_mash
        phenotype_mash_max = phenotype_mash
        p_value = 'null (only 1 strain grows)'#ttest_ind(spe_mash, phen_mash, equal_var=False)

    else:
        phen_mash = []
        phenotype_mash = 'null'
        phenotype_mash_min = 'null'
        phenotype_mash_max = 'null'
        p_value='null (no strains grow)'

    def add_clusters(node):
        if 'name' in node:
            node['cluster'] = cluster_data.get(node['name'], None)
        if 'children' in node:
            for child in node['children']:
                add_clusters(child)
        return node

    tree_data = add_clusters(tree_data)

    return jsonify({
        "tree_data": tree_data,
        "phe_mash": phenotype_mash,
        "phe_min": phenotype_mash_min ,
        "phe_max": phenotype_mash_max,
        "phe_individual_points":phen_mash,
        "kinetic_means":kinetic_means,
        "kinetic_errors":kinetic_errors,
        "all_strains":categories,
        "pval": p_value,
        "growth_strains":growth_strains,
        "nogrowth_strains":nogrowth_strains
    })



@app.route('/signal')
def signal():
    """
    Renders the signal page showing growth curves for a given plate, species, and well.

    Query Parameters:
        pltid (str): Plate ID.
        strn (str): Species name.
        well (str): Well name.

    Returns:
        Rendered HTML template 'signal.html' with growth data and time scale.
    """
    plateid = request.args.get('pltid')
    specie = request.args.get('strn')
    well = request.args.get('well')

    growth_data_entries = GrowthData.query.filter_by(plateid=plateid, specie=specie, well=well).all()
    growth_data = [
        {'name': f"{entry.compound} {entry.replicates}", 'data': entry.signal_data}
        for entry in growth_data_entries
    ]
    time_scale = list(range(0, 49, 1))
    #growth_data,time_scale = scripts.get_all_growth_curves(plateid,specie,well=well)
    return render_template('signal.html',growth_data=growth_data,time_scale=time_scale)


@app.route('/about',methods=['GET', 'POST'])
def about():
    """
    Renders the about page with control and growth well distributions.

    Returns:
        Rendered HTML template 'about.html' with control and growth well data.
    """
    control_wells,growth_wells=scripts.get_control_well_dist('pputida')
    #control_wells = random.sample(control_wells, 100)
    return render_template('about.html',control_wells = control_wells,growth_wells=growth_wells)

@app.route('/plates')
def plates():
    """
    Renders the plates page.

    Returns:
        Rendered HTML template 'plates.html'.
    """
    return render_template('plates.html')

@app.route('/ticket' ,methods=['GET', 'POST'])
def ticket():
    """
    Handles user support ticket submissions.

    POST:
        Receives name, email, and message from form and sends an email.

    GET:
        Renders the ticket submission form.

    Returns:
        Success message or rendered HTML template 'ticket.html'.
    """
    if request.method == 'POST':
        name = request.form['name']
        email = request.form['email']
        message = request.form['message']
        
        scripts.send_email(name, email, message)

        return 'Message sent successfully!'
    return render_template('ticket.html')


@app.route('/explore',methods=['GET', 'POST'])
def explore():
    """
    Renders the comparative analysis page for selected strains and compounds.

    POST:
        Processes selected strains and compound, returns comparative analysis.

    GET:
        Renders the explore page with available entries and options.

    Returns:
        Rendered HTML template 'comparative_analysis.html' or 'explore.html'.
    """
    if request.method == 'POST':
        # selected_entries = request.form.getlist('selected_entries')
        chosen_option = request.form.get('selected_option')
        selected_entries = request.form.getlist('selected_entries[]')
        plate,well = scripts.get_plate_well_from_compound(chosen_option)
        plateids = scripts.get_plateid_from_strain(selected_entries,plate)

        xlabels = scripts.get_strain_names(selected_entries)
        ylabels = [chosen_option]
        growth_calls,series,time = scripts.get_growth_calls_from_plateids(plateids,well,xlabels)
        return render_template('comparative_analysis.html',growth_calls=growth_calls,xlabels=xlabels,ylabels=ylabels,series=series,time=time)
    
    options = scripts.get_all_compounds_in_all_wells()
    entries = scripts.combine_specie_summaries()
    return render_template('explore.html',entries=entries,options=options)

@app.route('/plate_descriptions/json', methods=['GET'])
def plate_descriptions_json():
    """
    Returns plate descriptions in JSON format.

    Query Parameters:
        strain (str): Strain name.

    Returns:
        JSON: Plate description data.
    """
    strain = request.args.get('strain')
    plate_desc = pd.read_csv('./static/'+'plate_desc/platedesc.csv')

    out2 = []

    for i in plate_desc.index:
        plate = plate_desc.loc[i,'Plate']
        well = plate_desc.loc[i,'Well']
        compound = plate_desc.loc[i,'Compound']
        description = plate_desc.loc[i,'Description']
        kegg_id = plate_desc.loc[i,'KEGG ID']
        cas_id = plate_desc.loc[i,'CAS ID']

        out2.append([
            #str(plateid),
            str(plate),
            str(well),
            str(compound),
            str(description),
            # str(kegg_id),
            "<a href=https://www.genome.jp/entry/"+str(kegg_id)+">"+str(kegg_id)+"</a>",
            str(cas_id)])


    #return jsonify(data=out)
    return jsonify(data=out2)

@app.route('/species', methods=['GET'])
def species():
    """
    Renders the species page with metadata and summary information.

    Query Parameters:
        specie (str): Species name.

    Returns:
        Rendered HTML template 'species.html' with metadata.
    """
    specie = request.args.get('specie')
    specie_name = specie[0].upper() +'. '+specie[1:]
    samples,strains,available_plates,clusters,plates = scripts.load_specie_metadata(specie)
    
    return render_template('species.html',specie=specie,specie_name=specie_name,samples=samples,strains=strains,available_plates=available_plates,
                           clusters=clusters,plates=plates)

# Define a sorting key function
def sort_key(compound):
    """
    Sorting key for compounds based on well order.

    Args:
        compound (str): Compound string in format 'Well: Compound'.

    Returns:
        int: Index of the well in plate order.
    """
    wells = []
    for letter in range(ord('A'), ord('H') + 1):
        for num in range(1, 13):
            wells.append(chr(letter) + "{:02d}".format(num))
    well = compound.split(':')[0]
    return wells.index(well)


@app.route('/mainstraindata', methods=['GET'])
def mainstraindata():
    """
    Renders the main strain data page for a given plate and strain.

    Query Parameters:
        pltid (str): Plate ID.
        strn (str): Species name.
        plate (str): Plate name.
        strid (str): Strain ID.
        metadata (str): Metadata.
        media (str): Media.
        strain (str): Strain name.

    Returns:
        Rendered HTML template 'mainstraindata.html' with kinetic and growth data.
    """
    plateid = request.args.get('pltid')
    specie = request.args.get('strn')
    plate = request.args.get('plate')
    strid = request.args.get('strid')
    metadata = request.args.get('metadata')
    media = request.args.get('media')
    strain = request.args.get('strain')
    categories_list, mean_data, error_data,param_name = scripts.get_kinetic_parameters(plateid,specie,param='Max Resp')
    

    wells = []

    for letter in range(ord('A'), ord('H') + 1):
        for num in range(1, 13):
            wells.append(chr(letter) + "{:02d}".format(num))

    #growth_data,time_series,dropdown_names = scripts.get_all_growth_curves(plateid,specie,wells=['A01'])

    growth_data_entries = GrowthData.query.filter_by(plateid=plateid, specie=specie, well='A01').all()
    growth_data = [
        {'name': f"{entry.compound} {entry.replicates}", 'data': entry.signal_data}
        for entry in growth_data_entries
    ]

    # # Corrected code for filtering by a list of wells
    # growth_data_entries = GrowthData.query.filter(GrowthData.plateid == plateid,GrowthData.specie == specie,GrowthData.well.in_(wells)).all()

    # dropdown_compounds = list(set([entry.well +': ' +entry.compound for entry in growth_data_entries]))
    # dropdown_compounds = sorted(dropdown_compounds, key=sort_key)
    dropdown_compounds = scripts.get_compound_drop_down(plate)
    time_series = list(np.arange(0,48.25,0.25))

    return render_template('mainstraindata.html',pltid=plateid,strn=specie,categories = categories_list,mean_data = mean_data,error_data =
                           error_data,param_name=param_name,growth_data=growth_data,time_series=time_series,dropdown_compounds=dropdown_compounds,
                           strid=strid,media=media,metadata=metadata,plate=plate,strain=strain)

@app.route('/update_chart', methods=['GET'])
def update_chart():
    """
    Returns updated kinetic chart data for a given plate, species, and parameter.

    Query Parameters:
        pltid (str): Plate ID.
        strn (str): Species name.
        param (str): Kinetic parameter.

    Returns:
        JSON: Chart categories, mean data, error data, and parameter name.
    """
    plateid = request.args.get('pltid')
    specie = request.args.get('strn')
    param = request.args.get('param')
    categories, mean_data, error_data,param_name = scripts.get_kinetic_parameters(plateid, specie, param=param)
    return jsonify(categories=categories, mean_data=mean_data, error_data=error_data,param_name=param)

@app.route('/update_tracking_kinetics_chart', methods=['POST'])
def update_tracking_kinetics_chart():
    """
    Returns updated tracking kinetics chart data for selected strains and parameter.

    POST Data:
        growth_strains (list): List of strains with growth.
        no_growth_strains (list): List of strains without growth.
        param (str): Kinetic parameter.
        specie (str): Species name.
        plate (str): Plate name.
        well (str): Well name.

    Returns:
        JSON: Chart categories, mean data, error data, and parameter name.
    """
    data = request.get_json()
    growth_strains = data.get('growth_strains')
    no_growth_strains = data.get('no_growth_strains')
    param = data.get('param')
    specie = data.get('specie')
    plate = data.get('plate')
    well = data.get('well')

    mean_data,error_data,categories = scripts.get_tracking_kinetic_params(no_growth_strains, specie, plate, well, growth_strains, param)
    return jsonify(categories=categories, mean_data=mean_data, error_data=error_data,param_name=param)


@app.route('/update_growth_curve', methods=['GET'])
def update_growth_curve():
    """
    Returns updated growth curve data for a given plate, species, and compound.

    Query Parameters:
        pltid (str): Plate ID.
        strn (str): Species name.
        compound (str): Compound string.

    Returns:
        JSON: Growth curve data.
    """
    plateid = request.args.get('pltid')
    specie = request.args.get('strn')
    compound = request.args.get('compound')
    well = compound.split(':')[0]

    growth_data_entries = GrowthData.query.filter_by(plateid=plateid, specie=specie, well=well).all()
    growth_data = [
        {'name': f"{entry.compound} {entry.replicates}", 'data': entry.signal_data}
        for entry in growth_data_entries
    ]

    return jsonify(growth_data=growth_data)

@app.route('/straindata', methods=['GET'])
def straindata():
    """
    Renders the strain data heatmap page for a specific strain.

    Returns:
        Rendered HTML template 'straindata.html' with heatmap and compound data.
    """
    growth_calls,well_char,well_id,compound_dict = scripts.get_strain_data('ECP120')
    chart= {'type': 'heatmap','marginTop': 40,'marginBottom': 80,'plotBorderWidth': 1}
    title= {'text': ''}
    xAxis= {
        'categories': well_id,
        'labels':{'style':{'fontWeight':'bold','fontSize':'2em','fontFamily':'Monospace'}}
    }

    yAxis= {
        'categories': well_char,
        'title': 'null',
        'reversed': 'true',
        'labels':{'style':{'fontWeight':'bold','fontSize':'2em','fontFamily':'Monospace'}}
    }


    legend= {
        'enabled':'false',
        'align': 'right',
        'layout': 'vertical',
        'margin': 0,
        'verticalAlign': 'top',
        'y': 1,
        'symbolHeight': 280
    }

    series= [{
        'name': 'Growth(1)/No Growth(0)/Uncertain(0.5)',
        'borderWidth': 2.5,
        'borderColor':'#0a000f',
        'data': growth_calls,
        'dataLabels': {
            'enabled': 'false',
            'color': '#000000',
        }
    }]

    return render_template('straindata.html',chartID='container', chart=chart, data=growth_calls,
                           title=title,legend = legend,xAxis = xAxis,yAxis=yAxis,compound_dict = compound_dict)

@app.route('/strain_kinetics/json', methods=['GET'])
def strain_kinetics_json():
    """
    Returns kinetic parameters for a given plate and strain in JSON format.

    Query Parameters:
        spec (str): Strain name.
        plate (str): Plate ID.

    Returns:
        JSON: Kinetic parameter data.
    """
    strain = request.args.get('spec')
    plateid = request.args.get('plate')
    out2 = scripts.get_kinetic_parameters(plateid,strain)
    
    return jsonify(data=out2)


@app.route('/strain_growth/json', methods=['GET'])
def strain_growth_json():
    """
    Returns growth table for a given plate and strain in JSON format.

    Query Parameters:
        spec (str): Strain name.
        plate (str): Plate ID.

    Returns:
        JSON: Growth table data.
    """
    strain = request.args.get('spec')
    plateid = request.args.get('plate')
    out2 = scripts.get_growth_table(plateid,strain)
    
    return jsonify(data=out2)

@app.route('/get_growth_curves/json',methods=['POST'])
def get_growth_curves():
    """
    Returns growth curves for a given well, plate, and species.

    POST Data:
        well (str): Well name.
        plateid (str): Plate ID.
        specie (str): Species name.

    Returns:
        JSON: Chart data for growth curves.
    """
    well = request.form['well']
    plateid = request.form['plateid']
    specie = request.form['specie']
    chart_data = scripts.get_growth_curves(well,plateid,specie)
    return jsonify(chart_data) 


@app.route('/strains/json', methods=['GET'])
def strains_json():
    """
    Returns strain metadata in JSON format.

    Query Parameters:
        strain (str): Strain name.

    Returns:
        JSON: Strain metadata.
    """
    strain = request.args.get('strain')
    strain_data = pd.read_csv('./static/'+strain+'/metadata/updated_summary.csv')

    out2 = []

    for i in strain_data.index:
        plateid = strain_data.loc[i,'Plate IDs']
        id = strain_data.loc[i,'Strain ID']
        plate = strain_data.loc[i,'Plate']
        media = strain_data.loc[i,'Media']
        strain = strain_data.loc[i,'Strain']
        metadata = strain_data.loc[i,'Metadata/Modifications']
        phylo = strain_data.loc[i,'Phylogroup/Genome Cluster']
        mlst = strain_data.loc[i,'MLST']
        project = strain_data.loc[i,'Project']
        temperature = strain_data.loc[i,'Temperature']
        respiration = strain_data.loc[i,'Respiration']
        marker = strain_data.loc[i,'Selection Marker']
        reader = strain_data.loc[i,'Plate Reader']
        mode = strain_data.loc[i,'Detection mode']
        replicates = strain_data.loc[i,'Replicates']


        out2.append([
            #"<a href="+url_for('mainstraindata',pltid=str(plateid),strn=request.args.get('strain'))+">"+str(plateid)+"</a>",
            str(plateid),
            #str(plateid),
            str(id),
            str(plate),
            str(media),
            str(strain),
            str(metadata),
            str(phylo),
            str(mlst),
            str(project),
            str(temperature),
            str(respiration),
            str(marker),
            str(reader),
            str(mode),
            str(replicates)])


    #return jsonify(data=out)
    return jsonify(data=out2)

@app.route('/projects/json', methods=['GET'])
def projects():
    """
    Returns project metadata for a given strain in JSON format.

    Query Parameters:
        strain (str): Strain name.

    Returns:
        JSON: Project metadata.
    """
    strain = request.args.get('strain')
    project_data = pd.read_csv('./static/'+strain+'/metadata/project_summary.csv')
    out2 = []

    for i in project_data.index:
        project = project_data.loc[i,'Project']
        description = project_data.loc[i,'Description']
        doi = project_data.loc[i,'DOI']

        out2.append([
            str(project),
            str(description),
            "<a href=https://"+str(doi)+">"+str(doi)+"</a>"])


    #return jsonify(data=out)
    return jsonify(data=out2)
    

@app.route('/dashboard/strains', methods=['GET'])
def dashboard_strains():
    """
    Returns dashboard summary of strains per species in JSON format.

    Returns:
        JSON: Species and number of strains.
    """
    total_strains = strain.strain_summary()
    strain = request.args.get('strain')
    strain_data = pd.read_csv('./static/'+strain+'/metadata/summary.csv')

    out2 = []

    for i in total_strains.index:
        num_specie = total_strains.loc[i,'Num Strains']
        specie = i

        out2.append([
            str(specie),
            str(num_specie),])


    #return jsonify(data=out)
    return jsonify(data=out2)


@app.route('/compound_summary')
def compound_summary():
    """
    Renders the compound summary page for a given plate, well, and compound.

    Query Parameters:
        plate (str): Plate name.
        well (str): Well name.
        compound (str): Compound name.
        desc (str): Description.

    Returns:
        Rendered HTML template 'compound_summary.html'.
    """
    plate = request.args.get('plate')
    well = request.args.get('well')
    compound = request.args.get('compound')
    desc = request.args.get('desc')
    return render_template('compound_summary.html',plate=plate,well=well,compound=compound,desc=desc)

@app.route('/compound_summary_growth/json', methods=['GET'])
def compound_summary_growth_json():
    """
    Returns growth summary for a compound in JSON format.

    Query Parameters:
        plate (str): Plate name.
        well (str): Well name.

    Returns:
        JSON: Growth summary and species percentages.
    """
    plate = request.args.get('plate')
    well = request.args.get('well')

    growth_data_entries = TraitData.query.filter_by(plate=plate,well=well,growth=1).all()
    out2 = [[entry.strain,entry.metadata_mods,entry.specie[0].upper()+'. '+entry.specie[1:],entry.media,
            entry.project,entry.phylo,entry.mlst] for entry in growth_data_entries]
    

    no_growth_data_entries = TraitData.query.filter_by(plate=plate,well=well,growth=0).all()

    
    

    species = [entry.specie for entry in growth_data_entries]
    unique_species = list(set(species))
    
    total_entries = [entry.specie for entry in growth_data_entries] + [entry.specie for entry in no_growth_data_entries]#len(species)
    species_percentage = [(specie, (species.count(specie) / total_entries.count(specie)) * 100) for specie in unique_species]

    unique_species_list = [item[0][0].upper()+'. '+item[0][1:] for item in species_percentage]
    species_percentage_list = [round(item[1],2) for item in species_percentage]

    response = {
        "data": out2,
        "species_percentage_list": species_percentage_list,
        "unique_species_list": unique_species_list
    }



    return jsonify(response)



@app.route('/compound_summary_nogrowth/json', methods=['GET'])
def compound_summary_nogrowth_json():
    """
    Returns no-growth summary for a compound in JSON format.

    Query Parameters:
        plate (str): Plate name.
        well (str): Well name.

    Returns:
        JSON: No-growth summary and species percentages.
    """
    plate = request.args.get('plate')
    well = request.args.get('well')

    no_growth_data_entries = TraitData.query.filter_by(plate=plate,well=well,growth=0).all()
    growth_data_entries = TraitData.query.filter_by(plate=plate,well=well,growth=1).all()

             
    out2 = [[entry.strain,entry.metadata_mods,entry.specie[0].upper()+'. '+entry.specie[1:],entry.media,
            entry.project,entry.phylo,entry.mlst] for entry in no_growth_data_entries]
    

    species = [entry.specie for entry in no_growth_data_entries]
    unique_species = list(set(species))
    
    total_entries = [entry.specie for entry in growth_data_entries] + [entry.specie for entry in no_growth_data_entries]#len(species)
    species_percentage = [(specie, (species.count(specie) / total_entries.count(specie)) * 100) for specie in unique_species]

    unique_species_list = [item[0][0].upper()+'. '+item[0][1:] for item in species_percentage]
    species_percentage_list = [round(item[1],2) for item in species_percentage]

    response = {
        "data": out2,
        "species_percentage_list": species_percentage_list,
        "unique_species_list": unique_species_list
    }


    
    return response


required_columns = [
    'Plate Type', 'Media', 'Strain', 'Specie', 'Hr',
    'A01', 'A02', 'A03', 'A04', 'A05', 'A06', 'A07', 'A08', 'A09', 'A10', 'A11', 'A12',
    'B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09', 'B10', 'B11', 'B12',
    'C01', 'C02', 'C03', 'C04', 'C05', 'C06', 'C07', 'C08', 'C09', 'C10', 'C11', 'C12',
    'D01', 'D02', 'D03', 'D04', 'D05', 'D06', 'D07', 'D08', 'D09', 'D10', 'D11', 'D12',
    'E01', 'E02', 'E03', 'E04', 'E05', 'E06', 'E07', 'E08', 'E09', 'E10', 'E11', 'E12',
    'F01', 'F02', 'F03', 'F04', 'F05', 'F06', 'F07', 'F08', 'F09', 'F10', 'F11', 'F12',
    'G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09', 'G10', 'G11', 'G12',
    'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07', 'H08', 'H09', 'H10', 'H11', 'H12'
]

def allowed_file(filename):
    """
    Checks if the uploaded filename is allowed (CSV or XLSX).

    Args:
        filename (str): The filename to check.

    Returns:
        bool: True if allowed, False otherwise.
    """
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in {'csv', 'xlsx'}

@app.route('/upload')
def upload():
    """
    Renders the upload page for data files.

    Returns:
        Rendered HTML template 'upload.html'.
    """
    return render_template('upload.html')


@app.route('/upload_file', methods=['POST'])
def upload_file():
    """
    Handles file uploads, validates columns and missing data, processes and saves datasets.

    POST:
        Receives files via 'file[]'.

    Returns:
        Rendered HTML template 'processed_uploads.html' on success.
        Error message if missing columns or null entries.
    """
    if 'file[]' not in request.files:
        return 'No file part', 400
    files = request.files.getlist('file[]')

    uploaded_datasets = []
    missing_data_files = {}
    missing_column_list = []
    nan_dataframes = {}

    for file in files:
        if file and allowed_file(file.filename):
            filename = file.filename
            file_bytes = file.read()

            if filename.endswith('.csv'):
                df = pd.read_csv(BytesIO(file_bytes))
            else:
                df = pd.read_excel(BytesIO(file_bytes))

            # Check for required columns

            df = scripts.standardize_columns(df)
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                missing_data_files[filename] = missing_columns
                missing_column_list.append(missing_columns)
                #return f'Missing columns in {filename}: {", ".join(missing_columns)}, pls curate dataset/s', 400
            elif df.isna().any().any():
                nan_dataframes[filename] = 'missing'
            else:
                uploaded_datasets.append(df)


    if missing_data_files:
        missing_files_info = "\n".join([f"{file}: {', '.join(columns)}" for file, columns in missing_data_files.items()])
        return f'Missing columns in the following files:\n{missing_files_info}\nPlease curate the dataset(s).', 400

    if nan_dataframes:
        nan_files_info = "\n".join([f"{file}: {', '}\n" for file, columns in nan_dataframes.items()])
        return f'Missing/Null entries in columns in the following files:\n{nan_files_info}\nPlease curate the dataset(s).', 400


    plate_datatype = utils.make_plate_datatype(uploaded_datasets)
    kinetic_datatype = utils.get_kinetic_dataframe(plate_datatype)
    kinetic_datatype = utils.make_growth_calls(kinetic_datatype,max_resp_threshold=120,alpha=0.05)
    summary_table = utils.make_summary_table(kinetic_datatype)

    ### Add in plateids to kinetic and plate datatypes from summary table
    plate_datatype = pd.merge(plate_datatype,summary_table, on=['Strain','Specie','Plate','Media', 'Replicates'], how='left').set_index('PlateIDs')
    kinetic_datatype = pd.merge(kinetic_datatype,summary_table, on=['Strain','Specie','Plate','Media', 'Replicates'], how='left').set_index('PlateIDs')
    

    unique_key = str(uuid.uuid4())
    # session['unique_key'] = unique_key
    static_dir = os.path.join(app.root_path, 'static','cache', unique_key)
    os.makedirs(static_dir)

    # Save dataframes to the folder
    summary_table_path = os.path.join(static_dir, 'summary_table.csv')
    kinetic_datatype_path = os.path.join(static_dir, 'kinetic_datatype.csv')
    plate_datatype_path = os.path.join(static_dir, 'plate_datatype.csv')

    summary_table.to_csv(summary_table_path)
    kinetic_datatype.to_csv(kinetic_datatype_path)
    plate_datatype.to_csv(plate_datatype_path)

    session['unique_key'] = unique_key

    pca_dict,plate_list,pca_xaxis,pca_yaxis,heatmap_data,heatmap_y,heatmap_x = utils.PCA_analysis(kinetic_datatype)

    
    return render_template('processed_uploads.html',unique_key=unique_key,pca_dict = pca_dict,plate_list=plate_list,def_plate=plate_list[0],
                           pca_xaxis=pca_xaxis,pca_yaxis=pca_yaxis,heatmap_data=heatmap_data,heatmap_y=heatmap_y,heatmap_x=heatmap_x)


@app.route('/delete_session_data')
def delete_session_data():
    """
    Deletes cached session data for the current user.

    Returns:
        str: Success message.
    """
    unique_key = session.get('unique_key')
    if unique_key:
        static_dir = os.path.join(app.root_path, 'static','cache',unique_key)
        if os.path.exists(static_dir):
            shutil.rmtree(static_dir)
        session.pop('unique_key', None)
    return 'Session data deleted', 200


@app.before_request
def session_management():
    """
    Sets session to be permanent and configures session lifetime.

    Returns:
        None
    """
    session.permanent = True
    app.permanent_session_lifetime = timedelta(minutes=30)  # Adjust the lifetime as needed


@app.route('/upload_example/json', methods=['GET'])
def upload_example_json():
    """
    Returns example upload data in JSON format.

    Returns:
        JSON: Example upload data.
    """
    
    out2 = scripts.get_example_upload_data()

    return jsonify(data=out2)

@app.route('/summary_data_upload', methods=['GET'])
def summary_data_upload():
    """
    Returns summary table for uploaded data in JSON format.

    Returns:
        JSON: Summary table data.
    """
    key = session.get('unique_key', [])
    out2 = utils.get_uploaded_summary_table(key)
    return jsonify(data=out2)



@app.route('/uploaded_data_mainstraindata', methods=['GET'])
def uploaded_data_mainstraindata():
    """
    Renders the main strain data page for uploaded datasets.

    Query Parameters:
        plateids (str): Plate IDs.
        strain (str): Strain name.
        plate (str): Plate name.
        specie (str): Species name.
        media (str): Media.
        replicates (str): Replicate count.

    Returns:
        Rendered HTML template 'uploaded_data_mainstraindata.html'.
    """
    plateids = request.args.get('plateids')
    strain = request.args.get('strain')
    plate = request.args.get('plate')
    specie = request.args.get('specie')
    media = request.args.get('media')
    replicate = request.args.get('replicates')

    

    key = session.get('unique_key', [])
    if not key:
        return jsonify(data=[]), 400 
    
    categories_list, mean_data,param_name = utils.get_kinetic_parameters_for_sample(plateids,key,param='Max Resp')

    growth_data = utils.get_growth_curves_for_samples(plateids,key,well='A01')
    time_series = list(np.arange(0,48.25,0.25))
    dropdown_compounds = scripts.get_compound_drop_down(plate)
    return render_template('uploaded_data_mainstraindata.html',plateids=plateids,strain=strain,specie=specie,plate=plate,media=media,replicate=replicate,
                           categories = categories_list,mean_data = mean_data,param_name=param_name,growth_data=growth_data,
                           dropdown_compounds=dropdown_compounds,time_series=time_series)


@app.route('/uploaded_strain_growth/json', methods=['GET'])
def uploaded_strain_growth_json():
    """
    Returns growth table for uploaded datasets in JSON format.

    Query Parameters:
        plateids (str): Plate IDs.

    Returns:
        JSON: Growth table data.
    """
    key = session.get('unique_key', [])
    if not key:
        return jsonify(data=[]), 400  # If key is not found in session, return empty data with bad request

    
    plateids = request.args.get('plateids')
    out2 = utils.get_growth_table(plateids,key)

    return jsonify(data=out2)




@app.route('/update_kinetic_chart_uploaded_data', methods=['GET'])
def update_kinetic_chart_uploaded_data():
    """
    Returns updated kinetic chart data for uploaded datasets.

    Query Parameters:
        param (str): Kinetic parameter.
        plateids (str): Plate IDs.

    Returns:
        JSON: Chart categories, mean data, and parameter name.
    """
    param = request.args.get('param')
    plateids = request.args.get('plateids')
    key = session.get('unique_key', [])
    if not key:
        return jsonify(data=[]), 400 
    
    categories_list, mean_data,param_name = utils.get_kinetic_parameters_for_sample(plateids,key,param)
    
    return jsonify(categories=categories_list, mean_data=mean_data,param_name=param)


@app.route('/update_uploaded_sample_growth_curve', methods=['GET'])
def update_uploaded_sample_growth_curve():
    """
    Returns updated growth curve data for uploaded samples.

    Query Parameters:
        compound (str): Compound string.
        plateids (str): Plate IDs.

    Returns:
        JSON: Growth curve data.
    """
    compound = request.args.get('compound')
    well = compound.split(':')[0]
    plateids = request.args.get('plateids')

    key = session.get('unique_key', [])
    if not key:
        return jsonify(data=[]), 400 

    growth_data = utils.get_growth_curves_for_samples(plateids,key,well)

    return jsonify(growth_data=growth_data)


@app.route('/download_all_processed_upload_data', methods=['GET'])
def download_all_processed_upload_data():
    """
    Downloads all processed upload data as a zip file.

    Returns:
        Response: Sends zip file as an attachment.
    """
    key = session.get('unique_key', [])
    if not key:
        return jsonify(data=[]), 400

    download_directory = os.path.join('static', 'cache', key)
    output_filename = os.path.join('static', 'download_cache', key)

    os.makedirs(output_filename,exist_ok=True)

    # Compress the directory into a zip file
    shutil.make_archive(output_filename, 'zip', download_directory)

    zip_path = output_filename+'.zip'#os.path.join(output_filename,'.zip')

    #return send_file(zip_path, as_attachment=True, download_name='all_sequences.zip')

    # Send the zip file as a downloadable response
    return send_file(zip_path, as_attachment=True,download_name='all_processed_PMdata.zip')


@app.route('/download_sample_processed_upload_data', methods=['GET'])
def download_sample_processed_upload_data():
    """
    Downloads processed upload data for a specific sample as a zip file.

    Query Parameters:
        plateids (str): Plate IDs.
        strain (str): Strain name.
        plate (str): Plate name.

    Returns:
        Response: Sends zip file as an attachment.
    """
    plateids = request.args.get('plateids')
    strain = request.args.get('strain')
    plate = request.args.get('plate')

    key = session.get('unique_key', [])
    if not key:
        return jsonify(data=[]), 400

    signals = pd.read_csv('static/cache/' + key + '/plate_datatype.csv', index_col='PlateIDs')
    signals = signals.loc[plateids]

    kinetics = pd.read_csv('static/cache/' + key + '/kinetic_datatype.csv', index_col='PlateIDs')
    kinetics = kinetics.loc[plateids]

    # Create in-memory string buffers
    signals_csv = io.StringIO()
    kinetics_csv = io.StringIO()

    # Convert the filtered dataframes to CSV
    signals.to_csv(signals_csv)
    kinetics.to_csv(kinetics_csv)

    # Create an in-memory zip file
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w') as z:
        z.writestr(f'{strain}_{plate}_signals_data.csv', signals_csv.getvalue())
        z.writestr(f'{strain}_{plate}_kinetics_data.csv', kinetics_csv.getvalue())

    zip_buffer.seek(0)

    # Create a response object and set the appropriate headers
    response = Response(
        zip_buffer,
        mimetype='application/zip',
        headers={
            'Content-Disposition': f'attachment;filename={strain}_{plate}_data.zip'
        }
    )

    return response

#### Interop-DB Queries
def _parse_ids(data: dict, key: str):
    """Extract a list of IDs from *data* under *key*.

    Raises:
        ValueError: if the key is missing or the value is not a list.
    """
    if key not in data:
        raise ValueError(f"Missing required key '{key}'.")
    if not isinstance(data[key], list):
        raise ValueError(f"'{key}' must be a list – got {type(data[key]).__name__} instead.")
    return data[key]

def _row_to_dict(row):
    """Generic SQLAlchemy → dict (all simple columns)."""
    return {c.name: getattr(row, c.name) for c in row.__table__.columns}

def _growth_row_to_dict(row):
    """GrowthData → dict, making signal_data JSON-safe."""
    d = _row_to_dict(row)
    if isinstance(d.get("signal_data"), list):
        d["signal_data"] = [float(v) for v in d["signal_data"]]
    return d

def _build_strain_entries(strain_id):
    """Return a flat list of dicts — one per (strain, plate) combo — with full params and url."""
    rows = (db.session.query(
                KineticData.plateid,
                KineticData.specie,
                KineticData.plate,
                KineticData.media,
                KineticData.strain,
                KineticData.metadata_mods,
            )
            .filter_by(strainid=strain_id)
            .distinct()
            .all())

    entries = []
    for r in rows:
        params = urlencode({
            "pltid":    r.plateid,
            "strn":     r.specie,
            "plate":    r.plate,
            "media":    r.media or "",
            "strid":    strain_id,
            "metadata": r.metadata_mods or "",
            "strain":   r.strain or "",
        })
        entries.append({
            "strain":   strain_id,
            "plateid":  r.plateid,
            "plate":    r.plate,
            "specie":   r.specie,
            "media":    r.media or "",
            "metadata": r.metadata_mods or "",
            "url":      f"/mainstraindata?{params}",
        })
    return entries

@app.route("/interop-query/query-by-strain", methods=["POST"])
@cross_origin()
def query_by_strain():
    """
    POST body: {"ids": ["S1", "S2", ...]}
    Returns, for each strain ID, all KineticData and TraitData rows
    whose (strainid, plate) match that strain and any plate used
    by that strain.  The plate list is gathered once from KineticData.
    """
    logger.info("query by strain")

    try:
        payload = request.get_json()
        ids = _parse_ids(payload, "ids")

        def _row_to_dict(row):
            return {c.name: getattr(row, c.name) for c in row.__table__.columns}

        entries = []

        unique_plate_ids = set()
        for strain_id in ids:
            plate_ids = {p for (p,) in (
                db.session.query(KineticData.plateid)
                          .filter_by(strainid=strain_id)
                          .distinct()
            )}

            unique_plate_ids |= plate_ids

            # If a strain never occurs, return empty lists instead of error.
            if not plate_ids:
                entries.append({
                    "strainid": strain_id,
                    "urls":     [],
                    "kinethicdata": [],
                    "traitdata":   []
                })
                continue

            kin_rows = (KineticData.query
                        .filter_by(strainid=strain_id)
                        .filter(KineticData.plateid.in_(plate_ids))
                        .all())

            trait_rows = (TraitData.query
                          .filter_by(strainid=strain_id)
                          .filter(TraitData.plateid.in_(plate_ids))
                          .all())

            entries.append({
                "strainid": strain_id,
                "urls":     _build_strain_urls(strain_id),
                "kinethicdata": [_row_to_dict(r) for r in kin_rows],
                "traitdata":    [_row_to_dict(r) for r in trait_rows],
            })


        return jsonify({
            "entries": entries,
        }), 200

    except Exception as exc:
        logger.exception("Error in query_by_strain")
        return jsonify({"error": str(exc)}), 400


@app.route("/interop-query/strains", methods=["GET"])
@cross_origin()
def get_all_strains():
    """
    GET /strains
    Returns a list of all unique strain IDs from the KineticData table.
    """
    logger.info("get all strains")
    
    try:
        strain_ids = db.session.query(KineticData.strainid).distinct().all()

        strains = []
        for (strain_id,) in strain_ids:
            strains.extend(_build_strain_entries(strain_id))

        return jsonify({"strains": strains}), 200

    except Exception as exc:
        logger.exception("Error in get_all_strains")
        return jsonify({"error": str(exc)}), 400
    

if __name__ == "__main__":
    with app.app_context():
         #db.drop_all()
         #db.create_all()
         species_list = ['ecoli', 'pputida', 'saureus']
         for specie in species_list:
             print(specie)
             #ingest_data(specie)
             #ingest_trait_data(specie)
             #ingest_kinetic_data(specie)
    app.run(debug=True)

