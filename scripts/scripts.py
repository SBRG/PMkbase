import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from flask import Flask, render_template, url_for, jsonify
from flask import request

species = ['ecoli','pputida','saureus']

def strain_summary():

    strain_summary = pd.DataFrame(index=species)
    total_strains = []
    for specie in species:
        comp_strains = []
        temp_summary = pd.read_csv(url_for('.static',filename=specie+'/metadata/summary.csv'))
        #temp_summary = pd.read_csv('../static/'+specie+'/metadata/summary.csv')
        strains = temp_summary['Strain ID']
        mods = temp_summary['Metadata/Modification']
        for st in range(0,len(strains)):
            comp_strains.append(strains[st]+mods[st])
        total_strains.append(len(list(set(comp_strains))))

    strain_summary['Num Strains'] = total_strains

    return strain_summary


def strain_summary_json():

    total_strains = strain_summary()
    out2 = []

    for i in total_strains.index:
        num_specie = total_strains.loc[i,'Num Strains']
        specie = i

        out2.append([
            str(specie),
            str(num_specie),])


    #return jsonify(data=out)
    return jsonify(data=out2)
