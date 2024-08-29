from flask import Flask, render_template, url_for, jsonify, send_file, Response, session
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
import zipfile
import io




app = Flask(__name__,template_folder='templates',static_url_path='/static')
app.secret_key = 'sbrg_omnilog'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///growth_data.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db.init_app(app)

# with app.app_context():
#     db.create_all()


def ingest_data(specie):
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

    zip_path = os.path.join('static', 'all_sequences.zip')
    
    return send_file(zip_path, as_attachment=True, download_name='all_sequences.zip')



@app.route('/download/all_pmdata')
def download_all_pmdata():

    zip_path = os.path.join('static', 'all_PM_data.zip')
    
    return send_file(zip_path, as_attachment=True, download_name='allPMdata.zip')


@app.route('/download/all_specie_sequences')
def download_all_specie_sequences():

    specie=request.args.get('specie')
    
    # Define the path for the zip file
    zip_path = os.path.join('static', specie, 'sequences.zip')
    
    return send_file(zip_path, as_attachment=True, download_name=specie+'_sequences.zip')



@app.route('/download/all_specie_pmdata')
def download_all_specie_pmdata():

    specie=request.args.get('specie')
    
    # Define the path for the zip file
    zip_path = os.path.join('static', specie, 'data.zip')
    
    return send_file(zip_path, as_attachment=True, download_name=specie+'_allPMdata.zip')



@app.route('/download/mainstrain_specie_sequences')
def download_mainstrain_specie_sequences():

    specie= request.args.get('specie')
    strain = request.args.get('strain')

    for filename in os.listdir(os.path.join('static', specie, 'sequences')):
        if filename.startswith(strain):
            file_path = os.path.join('static', specie, 'sequences',filename)
            break
    

    return send_file(file_path, as_attachment=True, download_name=specie+'_'+strain+'.fna')



@app.route('/download/mainstrain_specie_growthdata')
def download_mainstrain_specie_growthdata():

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
    traits,categories = scripts.get_trait_summary()

    return render_template('index.html',series=traits,categories=categories)

@app.route('/dashboard')
def dashboard():
    return render_template('index.html')



@app.route('/tree')
def tree():

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
    control_wells,growth_wells=scripts.get_control_well_dist('pputida')
    #control_wells = random.sample(control_wells, 100)
    return render_template('about.html',control_wells = control_wells,growth_wells=growth_wells)

@app.route('/plates')
def plates():

    return render_template('plates.html')

@app.route('/ticket' ,methods=['GET', 'POST'])
def ticket():
    if request.method == 'POST':
        name = request.form['name']
        email = request.form['email']
        message = request.form['message']
        
        scripts.send_email(name, email, message)

        return 'Message sent successfully!'
    return render_template('ticket.html')


@app.route('/explore',methods=['GET', 'POST'])
def explore():

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
    specie = request.args.get('specie')
    specie_name = specie[0].upper() +'. '+specie[1:]
    samples,strains,available_plates,clusters,plates = scripts.load_specie_metadata(specie)
    
    return render_template('species.html',specie=specie,specie_name=specie_name,samples=samples,strains=strains,available_plates=available_plates,
                           clusters=clusters,plates=plates)

# Define a sorting key function
def sort_key(compound):
    wells = []
    for letter in range(ord('A'), ord('H') + 1):
        for num in range(1, 13):
            wells.append(chr(letter) + "{:02d}".format(num))
    well = compound.split(':')[0]
    return wells.index(well)


@app.route('/mainstraindata', methods=['GET'])
def mainstraindata():
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
    plateid = request.args.get('pltid')
    specie = request.args.get('strn')
    param = request.args.get('param')
    categories, mean_data, error_data,param_name = scripts.get_kinetic_parameters(plateid, specie, param=param)
    return jsonify(categories=categories, mean_data=mean_data, error_data=error_data,param_name=param)

@app.route('/update_tracking_kinetics_chart', methods=['POST'])
def update_tracking_kinetics_chart():

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
    #plateid = 'ECP120'
    #strain = request.args.get('strain')
    strain = request.args.get('spec')
    plateid = request.args.get('plate')
    out2 = scripts.get_kinetic_parameters(plateid,strain)
    
    return jsonify(data=out2)


@app.route('/strain_growth/json', methods=['GET'])
def strain_growth_json():
    #plateid = 'ECP120'
    #strain = request.args.get('strain')
    strain = request.args.get('spec')
    plateid = request.args.get('plate')
    out2 = scripts.get_growth_table(plateid,strain)
    
    return jsonify(data=out2)

@app.route('/get_growth_curves/json',methods=['POST'])
def get_growth_curves():
    well = request.form['well']
    plateid = request.form['plateid']
    specie = request.form['specie']
    chart_data = scripts.get_growth_curves(well,plateid,specie)
    return jsonify(chart_data) 


@app.route('/strains/json', methods=['GET'])
def strains_json():
    strain = request.args.get('strain')
    strain_data = pd.read_csv('./static/'+strain+'/metadata/summary.csv')

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
            str(project),])


    #return jsonify(data=out)
    return jsonify(data=out2)

@app.route('/projects/json', methods=['GET'])
def projects():
    strain = request.args.get('strain')
    project_data = pd.read_csv('./static/'+strain+'/metadata/project_summary.csv')
    out2 = []

    for i in project_data.index:
        project = project_data.loc[i,'Project']
        description = project_data.loc[i,'Description']

        out2.append([
            str(project),
            str(description)])


    #return jsonify(data=out)
    return jsonify(data=out2)
    

@app.route('/dashboard/strains', methods=['GET'])
def dashboard_strains():

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
    plate = request.args.get('plate')
    well = request.args.get('well')
    compound = request.args.get('compound')
    desc = request.args.get('desc')
    return render_template('compound_summary.html',plate=plate,well=well,compound=compound,desc=desc)

@app.route('/compound_summary_growth/json', methods=['GET'])
def compound_summary_growth_json():
    #plateid = 'ECP120'
    #strain = request.args.get('strain')
    plate = request.args.get('plate')
    well = request.args.get('well')
    #compound = request.args.get('compound')

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
    #plateid = 'ECP120'
    #strain = request.args.get('strain')
    plate = request.args.get('plate')
    well = request.args.get('well')
    #compound = request.args.get('compound')

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
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in {'csv', 'xlsx'}

@app.route('/upload')
def upload():
    return render_template('upload.html')


@app.route('/upload_file', methods=['POST'])
def upload_file():
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
    unique_key = session.get('unique_key')
    if unique_key:
        static_dir = os.path.join(app.root_path, 'static','cache',unique_key)
        if os.path.exists(static_dir):
            shutil.rmtree(static_dir)
        session.pop('unique_key', None)
    return 'Session data deleted', 200


@app.before_request
def session_management():
    session.permanent = True
    app.permanent_session_lifetime = timedelta(minutes=30)  # Adjust the lifetime as needed


@app.route('/upload_example/json', methods=['GET'])
def upload_example_json():
    
    out2 = scripts.get_example_upload_data()

    return jsonify(data=out2)

@app.route('/summary_data_upload', methods=['GET'])
def summary_data_upload():
    key = session.get('unique_key', [])
    out2 = utils.get_uploaded_summary_table(key)
    return jsonify(data=out2)



@app.route('/uploaded_data_mainstraindata', methods=['GET'])
def uploaded_data_mainstraindata():

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

    key = session.get('unique_key', [])
    if not key:
        return jsonify(data=[]), 400  # If key is not found in session, return empty data with bad request

    
    plateids = request.args.get('plateids')
    out2 = utils.get_growth_table(plateids,key)

    return jsonify(data=out2)




@app.route('/update_kinetic_chart_uploaded_data', methods=['GET'])
def update_kinetic_chart_uploaded_data():

    param = request.args.get('param')
    plateids = request.args.get('plateids')
    key = session.get('unique_key', [])
    if not key:
        return jsonify(data=[]), 400 
    
    categories_list, mean_data,param_name = utils.get_kinetic_parameters_for_sample(plateids,key,param)
    
    return jsonify(categories=categories_list, mean_data=mean_data,param_name=param_name)



@app.route('/update_uploaded_sample_growth_curve', methods=['GET'])
def update_uploaded_sample_growth_curve():

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

if __name__ == "__main__":
    # with app.app_context():
    #     db.drop_all()
    #     db.create_all()
    #     species_list = ['ecoli', 'pputida', 'saureus']
    #     for specie in species_list:
    #         ingest_data(specie)
    #         ingest_trait_data(specie)
    #         ingest_kinetic_data(specie)
    app.run(debug=True)

