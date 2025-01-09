#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 16:48:47 2021

@author: Elise Droste (e.droste@uea.ac.uk)

This script contains the column names in the DataFrames that are used in the processing of the raw VINDTA data for 
dissolved inorganic matter and total alkalinity. 




"""

# Existing columns from original summary file (UEA VINDTAs)
sample = "Sample Name"
comments = "Comments"
dicflag = "DIC flag"
taflag = "TA flag"
analysistype = "Analysis Type"
analysisdate = "Date (analysis)"
analysistime = "Time (analysis)"
station = "Station" # 9999 for junk sea water, 8888 for CRM, 2222 for underway 
cast = "Cast"
niskin = "Niskin"
replicate = "Replicate (diff bottle)"
repeat = "Repeat(same bottle)"
numrepeats = "No. Of Repeats"
depth = "Depth (m)"
ctdtemp = "CTD Temp (C)"
salinity = "Salinity"
phosphate = "Phosphate (uMol/Kg)"
silicate = "Silicate (uMol/Kg)"
nitrate = "Nitrate (uMol/Kg)"
dicpipettevol = "DIC pipette volume (ml)"
dicep = "DIC set end point (cpm)"
dicpippettetemp = "DIC pipette temp sample (C)"
dicbackground = "DIC background (cpm)"
dictime = "DIC total time (s)"
diccounts = "DIC total counts"
diccellnum = "DIC cell Number (sequential)"
sampletemp = "Sample Bottle temp (C)"
tapipettevol = "TA pipette vol (ml)"
tacelltemp = "TA cell temp (C)"
taacidbatch = "TA Acid Batch"
taacidconc = "TA Acid conc (Mol)"
taaciddensity = "TA Acid density (Kg/l)"

# Added prior to creating this class
salinity_availability = "salinity_availability"
ctd_temp_availability  = "ctd_temp_availability"

sample_temp_flag = "sample_temp_flag"

sampling_date = "sampling_date"


# Columns added by PICCOLO_calc_DIC.py
crm_batch = "crm_batch"
analysis_datetime = "analysis_datetime"
sample_temperature_alternative = "sample_temperature_alternative"
dic_certified_umolkg = "DIC_certified_umolkg"
ta_certified_umolkg = "TA_certified_umolkg"
netdiccounts = "netdiccounts"
diccounts_volcorr = "diccounts_volcorr"
dic_sw_density = "dic_sw_density" 
dicpipmass = "dicpipmass"
dic_crm_pip_umol = "dic_crm_pip_umol"
dic_crm_counts_per_umol = "dic_crm_counts_per_umol"
dic_crm_counts_per_umol_mean_allcells = "dic_crm_counts_per_umol_mean_allcells"
dic_crm_counts_per_umol_std_allcells = "dic_crm_counts_per_umol_std_allcells"
dic_sample_pip_umol_allcells = "dic_sample_pip_umol_allcells"
dic_sample_umol_per_l_allcells = "dic_sample_umol_per_l_allcells"
dic_sample_umol_per_kg_allcells = "dic_sample_umol_per_kg_allcells"
dic_crm_counts_per_umol_mean_percell = "dic_crm_counts_per_umol_mean_percell"
dic_crm_counts_per_umol_std_percell = "dic_crm_counts_per_umol_std_percell"
dic_sample_pip_umol = "dic_sample_pip_umol"
dic_sample_umol_per_l = "dic_sample_umol_per_l" 
dic_sample_umol_per_kg = "dic_sample_umol_per_kg"
dic_crm_counts_per_umol_driftcorr = "dic_crm_counts_per_umol_driftcorr"
dic_sample_pip_umol_driftcorr = "dic_sample_pip_umol_driftcorr"
dic_sample_umol_per_l_driftcorr = "dic_sample_umol_per_l_driftcorr"
dic_sample_umol_per_kg_driftcorr = "dic_sample_umol_per_kg_driftcorr"
dic_driftcorr_flag = "dic_driftcorr_flag"

# Columns added elsewhere
dic_sample_umol_per_kg_final = "dic_sample_umol_per_kg_final"
ta_sample_umol_per_kg_final = "ta_sample_umol_per_kg_final"


# Columns added by PICCOLO_calc_TA.py
ta_acid_batch_number = "ta_acid_batch_number"
analysis_batch = "analysis_batch"
titrant_molinity_here = "titrant_molinity_here"
titrant_molinity = "titrant_molinity"
titrant_molinity_std_batch = "titrant_molinity_std_batch"
alkalinity_umolkg = "alkalinity_umolkg"

dic_for_ta = "dic_for_ta"# "dic_sample_umol_per_kg" ## 20240724 NO LONGER: Set which DIC values for the samples are to be used for the TA calculation. < 20240724 - ESD - AWI - changed from dic_sample_umol_per_kg to dic_for_ta >

# IDs
crmid = 8888 # CRM ID
junkid = 9999 # Junk sea water ID 
sampleid = "bottle"
# < 20240207 - ED - SDA - added to allow to select stations for profiles
profile_stations_list = []
uw_stations_list = []
# 20240207 - ED - SDA > 

# Styles
interactive = "on" # "on" or "off"; switch to plot interactive plotly plots


# Styles for qc_conc_profiles()
numberofcols = 6 # for number of columns of interactive profile plots 

cell_colours_list = ["firebrick", "darkorange", "blue",
                            "blueviolet", "brown", "burlywood", "cadetblue", "chartreuse", "chocolate", "coral", "darkolivegreen", "cornflowerblue",
                            "crimson", "cyan", "darkblue", "darkcyan","darkgoldenrod", "darkgray", "darkgreen",
                            "darkkhaki", "darkmagenta","darkorchid", "darkred", "darksalmon", "darkseagreen",
                            "darkslateblue", "darkslategray", "darkslategrey","darkturquoise", "darkviolet", "deeppink", "deepskyblue",
                            "dimgrey", "dodgerblue", "forestgreen", "cornsilk"]
cell_colours = "cell_colours"
coloured_var= diccellnum
fig_height = 2500
fig_width = 1000










