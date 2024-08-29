import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import re
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from flask import jsonify
from sklearn.decomposition import PCA
import math
from models import db, GrowthData, TraitData,KineticData

# All wells in the biolog plate
wells =     ['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12',
             'B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12',
             'C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12',
             'D01','D02','D03','D04','D05','D06','D07','D08','D09','D10','D11','D12',
             'E01','E02','E03','E04','E05','E06','E07','E08','E09','E10','E11','E12',
             'F01','F02','F03','F04','F05','F06','F07','F08','F09','F10','F11','F12',
             'G01','G02','G03','G04','G05','G06','G07','G08','G09','G10','G11','G12',
             'H01','H02','H03','H04','H05','H06','H07','H08','H09','H10','H11','H12']



def clean_strain(strain):
    # Remove special characters and spaces
    strain = str(strain)
    cleaned_strain = re.sub(r'[^A-Za-z0-9]+', '', strain)
    # Remove specified strings
    cleaned_strain = re.sub(r'R[1-9]', '', cleaned_strain)
    return cleaned_strain

def clean_plate(plate):
    temp_plate = re.sub(r'[^A-Za-z0-9]', '', plate)
    temp_plate = temp_plate.lower()
    if('pm01' in temp_plate or temp_plate=='pm1'):
        clean_plate = 'PM01'
    elif('pm02' in temp_plate or temp_plate=='pm2'):
        clean_plate = 'PM02A'
    elif('pm03' in temp_plate or temp_plate=='pm3'):
        clean_plate = 'PM03B'
    elif('pm04' in temp_plate or temp_plate=='pm4'):
        clean_plate = 'PM04A'
    elif('pm05' in temp_plate or temp_plate=='pm5'):
        clean_plate = 'PM05'
    elif('pm05' in temp_plate or temp_plate=='pm5'):
        clean_plate = 'PM05'
    elif('pm06' in temp_plate or temp_plate=='pm6'):
        clean_plate = 'PM06'
    elif('pm07' in temp_plate or temp_plate=='pm7'):
        clean_plate = 'PM07'
    elif('pm08' in temp_plate or temp_plate=='pm8'):
        clean_plate = 'PM08'
    elif('pm09' in temp_plate or temp_plate=='pm9'):
        clean_plate = 'PM09'
    elif('pm10' in temp_plate or temp_plate=='pm10'):
        clean_plate = 'PM10'
    elif('pm11' in temp_plate or temp_plate=='pm11c'):
        clean_plate = 'PM11C'
    elif('pm12' in temp_plate or temp_plate=='pm12b'):
        clean_plate = 'PM12B'
    return clean_plate

def get_plate_layout(plate):
    if(plate=='PM02'):
        plate='PM02A'
    plate_layout = pd.read_csv('static/plate_desc/platedesc.csv',index_col='Well')
    plate_layout = plate_layout[plate_layout['Plate']==plate]
    return plate_layout

def label_replicates(group):
    group['Replicates'] = 'R' + (group.groupby(['Strain', 'Well', 'Media','Plate']).cumcount() + 1).astype(str)
    return group

def make_plate_datatype(dataset):

    
    time_interval = np.linspace(0,48,193)
    combined_df = pd.DataFrame(columns=['Strain','Specie','Plate','Well', 'Media','Compound','KEGG ID','CAS ID'] + [f'{t}hrs' for t in time_interval])

    for data in dataset:

        strain = str(data['Strain'].unique().tolist().pop())
        media = str(data['Media'].unique().tolist().pop())
        plate = str(data['Plate Type'].unique().tolist().pop())
        specie = str(data['Specie'].unique().tolist().pop())
        
        plate_layout = get_plate_layout(plate)
        for well in wells:
            # Create a dictionary to store the data for the current file
            file_data = {
                'Strain': strain,
                'Specie': specie,
                'Well': well,
                'Media': media,
                'Plate': plate,
                'Compound': plate_layout.loc[well]['Compound'],
                'Description':plate_layout.loc[well]['Description'],
                'KEGG ID':plate_layout.loc[well]['KEGG ID'],
                'CAS ID':plate_layout.loc[well]['CAS ID']
            }

            # Iterate through the columns, which are the time points
            for col in range(0,len(time_interval)):
                # Extract the time value (e.g., '1hr', '1.5hr') and store it in the dictionary
                file_data[f'{time_interval[col]}hrs'] = data[well].values[col]

            # Append the data for the current file to the combined DataFrame
            combined_df = pd.concat([combined_df,pd.DataFrame([file_data])],ignore_index=True)

    combined_df['Strain'] = combined_df['Strain'].apply(clean_strain)
    combined_df['Plate'] = combined_df['Plate'].apply(clean_plate)

    

    combined_df['Replicates'] = ''
    # Sort the DataFrame by 'Strain', 'Well', and 'Media' to ensure consistent labeling
    #combined_df = combined_df.sort_values(['Strain', 'Well','Media'])

    # Group the DataFrame by 'Strain', 'Well', and 'Media' and apply the labeling function
    combined_df = combined_df.groupby(['Strain', 'Well','Media','Plate'], group_keys=False).apply(label_replicates)

    # Set the columns you want at the beginning
    column_order = ['Strain','Plate','Well', 'Media', 'Replicates','Compound','Description','KEGG ID','CAS ID']

    # Extract the remaining columns
    other_columns = [col for col in combined_df.columns if col not in column_order]

    # Combine the columns in the desired order
    new_column_order = column_order + other_columns

    # Reorder the DataFrame
    combined_df = combined_df[new_column_order]
    combined_df = combined_df.reset_index(drop=True)

    return combined_df

def get_kinetic_parameters(signal,time):
    """
    Extract the max resp, rate, time and auc for a respiration signal
    """
    smoothened_signal = savgol_filter(signal, 50, 3)
    dt = []
    for i in range(0,np.shape(smoothened_signal)[0]-1):
        dt.append((smoothened_signal[i+1]-smoothened_signal[i])/(time[i+1]-time[i]))
    max_resp_rate = np.max(dt)
    time_till = time[np.argmax(smoothened_signal)]
    max_val = np.max(smoothened_signal)
    auc = np.trapz(smoothened_signal,time,dx=0.25)   
    return max_val,max_resp_rate,time_till,auc

def get_kinetic_dataframe(plate_dataframe):
    
    kinetic_dataframe = pd.DataFrame(index=plate_dataframe.index,columns=['Strain','Specie','Well','Plate','Media','Replicates','Compound','Description','KEGG ID','CAS ID','Max Resp','Max Resp Rate','Time till max resp rate','AUC'])
    signal_columns = [str(w)+'hrs' for w in np.linspace(0,48,193)]
    for i in range(0,plate_dataframe.shape[0]):
        kinetic_dataframe.iloc[i,0] = plate_dataframe['Strain'].iloc[i]
        kinetic_dataframe.iloc[i,1] = plate_dataframe['Specie'].iloc[i]
        kinetic_dataframe.iloc[i,2] = plate_dataframe['Well'].iloc[i]
        kinetic_dataframe.iloc[i,3] = plate_dataframe['Plate'].iloc[i]
        kinetic_dataframe.iloc[i,4] = plate_dataframe['Media'].iloc[i]
        kinetic_dataframe.iloc[i,5] = plate_dataframe['Replicates'].iloc[i]
        kinetic_dataframe.iloc[i,6] = plate_dataframe['Compound'].iloc[i]
        kinetic_dataframe.iloc[i,7] = plate_dataframe['Description'].iloc[i]
        kinetic_dataframe.iloc[i,8] = plate_dataframe['KEGG ID'].iloc[i]
        kinetic_dataframe.iloc[i,9] = plate_dataframe['CAS ID'].iloc[i]        
        max_val,max_resp_rate,time_till,auc = get_kinetic_parameters(plate_dataframe.iloc[i,:][signal_columns],np.linspace(0,48,193))
        kinetic_dataframe.iloc[i,10] = max_val
        kinetic_dataframe.iloc[i,11] = max_resp_rate
        kinetic_dataframe.iloc[i,12] = time_till
        kinetic_dataframe.iloc[i,13] = auc

    return kinetic_dataframe



def make_growth_calls(kinetic_dataframe,max_resp_threshold=120,alpha=0.05,negative_control=True):

    global_kinetic_data = KineticData.query.filter_by(growth=0.0).all()
    out2 = [entry.maxresp for entry in global_kinetic_data]


    # Create an empty 'Growth' column and initialize with zeros
    #if(control_well_kinetics.shape[0]):
    kinetic_dataframe['Growth'] = 0
    kinetic_dataframe['Control Well Growth'] = 0

    # Create an empty list to store p-values
    p_values = []

    # Iterate through each row and perform a one-sided z-test
    for index, row in kinetic_dataframe.iterrows():
        # Perform a one-sided z-test using control group statistics
        #z_score = (row['Max Resp'] - control_group['Max Resp'].mean()) / (control_group['Max Resp'].std() / (len(control_group) ** 0.5))
        z_score = (row['Max Resp'] - np.mean(out2)) / (np.std(out2) / (len(out2) ** 0.5))
        p_value = 1 - norm.cdf(z_score)
        p_values.append(p_value)

    # Apply Benjamini-Hochberg correction to p-values
    p_adjusted = multipletests(p_values, method='fdr_bh')[1]

    for i, p_adj in enumerate(p_adjusted):
        if (p_adj < alpha and kinetic_dataframe['Max Resp'].iloc[i] > np.mean(out2) and kinetic_dataframe['Max Resp'].iloc[i]>max_resp_threshold):
            kinetic_dataframe.at[i, 'Growth'] = 1

    kinetic_dataframe['Control Well Growth'] = kinetic_dataframe.apply(lambda row: control_well_growth_condition(row, kinetic_dataframe), axis=1)

    return kinetic_dataframe


def control_well_growth_condition(row, df):
    if not df[(df['Compound'] == 'Negative Control') &
              (df['Strain'] == row['Strain']) &
              (df['Specie'] == row['Specie']) &
              (df['Plate'] == row['Plate']) &
              (df['Media'] == row['Media']) &
              (df['Replicates'] == row['Replicates']) &
              (df['Growth'] == 1)].empty:
        return 1
    else:
        return 0

def make_summary_table(kinetic_frame):

    summary = kinetic_frame[['Strain','Specie','Plate','Media','Replicates','Control Well Growth']].drop_duplicates(keep='first')
    plate_ids = ['ECP'+str(i) for i in range(0,summary.shape[0])]
    
    summary['PlateIDs'] = plate_ids
    #summary = summary.set_index('Plate IDs')

    return summary



def get_upload_kinetic_parameters(growth_frame,plateids,strain,plate,media,replicate,param='Max Resp'):

    growth_calls = growth_frame.loc[(growth_frame['Strain']==strain)&(growth_frame['Plate']==plate)&(growth_frame['Media']==media)&(growth_frame['Replicates']==replicate)]
    growth_calls = growth_frame.loc[plateids]
    categories = [w['Well']+': '+w['Compound'] for _,w in growth_frame[['Well','Compound']].drop_duplicates().iterrows()]
    numeric_cols = ['Well','Max Resp', 'Max Resp Rate', 'Time till max resp rate', 'AUC']
    params = growth_frame[numeric_cols].set_index('Well')#.apply(pd.to_numeric, errors='coerce')

    # Calculate average of numeric columns grouped by Well
    avg_df = params.groupby('Well').mean().reset_index()

    # Calculate max of numeric columns grouped by Well
    max_df = params.groupby('Well').max().reset_index()

    # Calculate min of numeric columns grouped by Well
    min_df = params.groupby('Well').min().reset_index()

    colors = ["#FFD1DC" if w == 1 else "#C0C0FF" if w == 0.5 else "#808080" for w in growth_calls['Growth']]

    main_data = []
    error_data = []

    for i in range(0,avg_df.shape[0]):
        main_data.append({'y':avg_df[param][i],'color':colors[i]})

    for i in range(0,max_df.shape[0]):
        error_data.append([round(max_df[param][i],2),round(min_df[param][i],2)])
    
    return categories,main_data,error_data,param



def get_uploaded_summary_table(session_id):
    summary_table = pd.read_csv('static/cache/'+session_id+'/summary_table.csv')
    out2 = []
    for _, row in summary_table.iterrows():
        out2.append({
            'PlateIDs': str(row['PlateIDs']),
            'Strain': str(row['Strain']),
            'Specie': str(row['Specie']),
            'Plate': str(row['Plate']),
            'Media': str(row['Media']),
            'Replicates': str(row['Replicates']),
            'Control Well Growth':str(row['Control Well Growth']),
        })
    print(out2)
    return out2




def get_growth_table(plateids,session_id):

    growth_frame = pd.read_csv('static/cache/'+session_id+'/kinetic_datatype.csv')
    growth_frame = growth_frame[growth_frame['PlateIDs']==plateids]
    
    
    
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



def get_kinetic_parameters_for_sample(plateids,session_id,param='Max Resp'):
    
    growth_frame = pd.read_csv('static/cache/'+session_id+'/kinetic_datatype.csv',index_col='PlateIDs')
    

    growth_frame = growth_frame.loc[plateids]
    
    
    categories = [w['Well']+': '+w['Compound'] for _,w in growth_frame[['Well','Compound']].drop_duplicates().iterrows()]
    numeric_cols = ['Well','Max Resp', 'Max Resp Rate', 'Time till max resp rate', 'AUC']
    params = growth_frame[numeric_cols].set_index('Well')#.apply(pd.to_numeric, errors='coerce')

    # Calculate average of numeric columns grouped by Well
    avg_df = params.groupby('Well').mean().reset_index()

    colors = ["#FFD1DC" if w == 1 else "#C0C0FF" if w == 0.5 else "#808080" for w in growth_frame['Growth']]

    main_data = []

    for i in range(0,avg_df.shape[0]):
        main_data.append({'y':avg_df[param][i],'color':colors[i]})

    
    return categories,main_data,param


def get_growth_curves_for_samples(plateids,session_id,well):

    growth_frame = pd.read_csv('static/cache/'+session_id+'/plate_datatype.csv',index_col='PlateIDs')
    growth_frame = growth_frame.loc[plateids]
    growth_frame = growth_frame[growth_frame['Well']==well]
    
    signal_columns = [str(w)+'hrs' for w in np.linspace(0,48,193)]
    
    #print(growth_frame[signal_columns].values[0].tolist())
    return [{'name':growth_frame['Well'].values[0]+': '+growth_frame['Compound'].values[0],'data':growth_frame[signal_columns].values[0].tolist()}]



def PCA_analysis(kinetic_datatype):

    kinetic_datatype['Strain Replicate'] = [w['Strain']+' '+w['Replicates']+'('+w['Media']+')' for _,w in kinetic_datatype.iterrows()]

    plate_list = kinetic_datatype['Plate'].unique().tolist()

    pca_dict = {}
    pca_xaxis = {}
    pca_yaxis = {}
    heatmap_data_dict = {}
    heatmap_y_dict = {}
    heatmap_x_dict = {}

    for plate in plate_list:

        df = kinetic_datatype[kinetic_datatype['Plate']==plate]
        df = df.pivot_table(index='Compound',columns='Strain Replicate',values='Growth').copy(deep=False)
        # df.to_csv('df.tsv',sep='\t')
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(df.T)
        pca_dict[plate] = [{'x':pca_result[i,0],'y':pca_result[i,1],'z':10,'name':df.columns[i]} for i in range(0,pca_result.shape[0])]
        pca_xaxis[plate] = 'Explained Var: '+str(round(pca.explained_variance_ratio_[0]*100,2))+'%'
        pca_yaxis[plate] = 'Explained Var: '+str(round(pca.explained_variance_ratio_[1]*100,2))+'%'

        loadings = pca.components_.T * np.sqrt(pca.explained_variance_)

        # Convert loadings to DataFrame for easier manipulation
        loadings_df = pd.DataFrame(loadings, index=df.index, columns=['PC1', 'PC2'])

        # Select top features based on absolute loadings for PC1 and PC2
        top_features_pc1 = loadings_df['PC1'].abs().nlargest(5).index.tolist()
        top_features_pc2 = loadings_df['PC2'].abs().nlargest(5).index.tolist()
        top_features = list(set(top_features_pc1+top_features_pc2))


        loadings = df.loc[top_features]
        heatmap_data,heatmap_y,heatmap_x = make_loadings_hchartshmapstyle(loadings)

        heatmap_data_dict[plate] = heatmap_data
        heatmap_y_dict[plate] = heatmap_y
        heatmap_x_dict[plate] = heatmap_x
        

    return pca_dict, plate_list, pca_xaxis, pca_yaxis, heatmap_data_dict,heatmap_y_dict,heatmap_x_dict



def make_loadings_hchartshmapstyle(loadings_df):

    data = []
    for i in range(0,loadings_df.shape[0]):
        for j in range(0,loadings_df.shape[1]):
            data.append([i,j,loadings_df.iloc[i,j]])

    yaxis_cats = loadings_df.columns.tolist()
    xaxis_cats = loadings_df.index.tolist()
    
    return data,yaxis_cats,xaxis_cats