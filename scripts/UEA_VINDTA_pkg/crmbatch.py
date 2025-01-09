#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 17:44:49 2019

CRM batches: https://www.nodc.noaa.gov/ocads/oceans/Dickson_CRM/batches.html

@author: Elise Droste (e.droste@uea.ac.uk)
"""

print("Imported crmbatch script.")

import numpy as np
import pandas as pd

def findcrmbatch(crm_batch, dickson_location): 
    lookup = pd.DataFrame(pd.read_csv(dickson_location))
    # units: umol kg-1
    
    if crm_batch in lookup["crm_batch"].values:
        crm_dickson_dic = np.float(lookup["dic"][lookup["crm_batch"]==crm_batch])
        crm_dickson_dic_sd = np.float(lookup["dic_sd"][lookup["crm_batch"]==crm_batch])
        crm_dickson_alk = np.float(lookup["alk"][lookup["crm_batch"]==crm_batch])
        crm_dickson_alk_sd = np.float(lookup["alk_sd"][lookup["crm_batch"]==crm_batch])
        crm_dickson_phos = np.float(lookup["phos"][lookup["crm_batch"]==crm_batch])
        crm_dickson_sil = np.float(lookup["sil"][lookup["crm_batch"]==crm_batch])
        crm_dickson_nitrite = np.float(lookup["nitrite"][lookup["crm_batch"]==crm_batch])
        crm_dickson_nitrate = np.float(lookup["nitrate"][lookup["crm_batch"]==crm_batch])
        crm_dickson_salinity = np.float(lookup["sal"][lookup["crm_batch"]==crm_batch])
            
        crm_dickson = {"dic": crm_dickson_dic, "dic_sd": crm_dickson_dic_sd, "alk": crm_dickson_alk, "alk_sd": crm_dickson_alk_sd,  
                           "phos": crm_dickson_phos, "sil": crm_dickson_sil, "nitrite": crm_dickson_nitrite, "nitrate": crm_dickson_nitrate, 
                           "salinity": crm_dickson_salinity}
        
        return crm_dickson
    
    else: 
        print("No info on other CRM batches currently provided in the lookup. Update crm_batch_lookup.csv. Prepare for error message.")

