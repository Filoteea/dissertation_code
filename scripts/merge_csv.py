# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 03:00:26 2022

@author: filot
"""
import pandas as pd
import create_model
import numpy as np

# merge csvs of models
def merge_csv(date):
    if date == '2018-01':
        files = ['outputs/models/models_2018-01.csv', 
                 'outputs/models/models_2018-02.csv', 
                 'outputs/models/models_2018-03.csv', 
                 'outputs/models/models_2018-04.csv', 
                 'outputs/models/models_2018-05.csv', 
                 'outputs/models/models_2018-06.csv']
    else:
        files = ['outputs/models/models_2018-07.csv', 
                 'outputs/models/models_2018-08.csv', 
                 'outputs/models/models_2018-09.csv', 
                 'outputs/models/models_2018-10.csv', 
                 'outputs/models/models_2018-11.csv', 
                 'outputs/models/models_2018-12.csv']
    merged = pd.concat([pd.read_csv(f) for f in files])
    merged.to_csv( f"outputs/models/merged_{date}.csv", index=False, encoding='utf-8-sig')
    return ('done')
#merge_csv('2018-01')
#merge_csv('2018-07')

def merge_models(date):
    if date == '2018-01':
        dates = ['2018-01', '2018-02', '2018-03', '2018-04', '2018-05', '2018-06']
    else:
        dates = ['2018-07', '2018-08', '2018-09', '2018-10', '2018-11', '2018-12']  
    merged = np.extend((create_model.run_model('2018-01'), create_model.run_model('2018-02')))
    '''
    for i in dates:
        model = create_model.run_model(i)
        merged=np.append(model)
    '''
    return merged
#mod = merge_models('2018-02')

def run_models(date):
    if date == '2018-01':
        dates = ['2018-01', '2018-02', '2018-03', '2018-04', '2018-05', '2018-06']
    else:
        dates = ['2018-07', '2018-08', '2018-09', '2018-10', '2018-11', '2018-12']
    for i in dates:
        create_model.run_model(i)
    
date = '2018-07'
run_models(date)





# Combine all three CSV files using the concat method
#merged = pd.concat([pd.read_csv(f) for f in files])
# Export to csv
#merged.to_csv( "outputs/models/merged.csv", index=False, encoding='utf-8-sig')


