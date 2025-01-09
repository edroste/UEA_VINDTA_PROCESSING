#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 11:14:21 2021

@author: Elise Droste (e.droste@uea.ac.uk)

Conventions for this script to work: 
    - Station = 9999 for junk seawater runs, 8888 for CRM runs
    - Cast = yyyymmdd of analysis for junk seawater runs, CRM batch number for CRM runs 
    - Niskin = junk number of the day for junks, CRM bottle numbre for CRMs 
    - Depth = 0 for junks and CRMs
    - Repeat (same bottle) = automatic repeat run (set by software)
    - Replicate (different bottle) = 1 for junk and CRM runs, unique bottle number for samples

"""


print("Hello! Before running the CalcDIC Class, make sure the filenames are correct.")

# Import required built-in packages
# Consider being more specific for packages e.g. from numpy import max as np_max

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from . import crmbatch

# import plotly
import plotly.graph_objs as go


# Create custom errors 
class Error(Exception):
    """Base class for other exceptions"""
    pass

class IncorrectArgumentError(Error):
    """Raised when incorrect arguments are given to my methods"""
    pass



# Generate a class based on the Pandas DataFrame (df = csv file)
class CalcDIC: 
    """
    TO DO 
    """
    
    #######################################################################################################################
    def __init__(self, df):
        """Initialise the Class by setting the fields as instances. Also do some cleaning up. """
        
        # ---------------------------------------------------------------------
        # Some cleaning up 
        # Copy data frame and remove the junk sea water sample runs + reset the index
        self.df = df[df["Station"]!=9999].reset_index(drop=True)
        print("All junk samples (Station = 9999) have been removed.")
        
        # Add some new columns (others will be created along the way) 
        self.df["crm_batch"] = pd.Series([])
        self.df["DIC_certified_umolkg"] = pd.Series([]) # ED - 20210512 - Norwich - added > 
        self.df["TA_certified_umolkg"] = pd.Series([]) # ED - 20210512 - Norwich - added > 
        # < 20210511 - ED - Norwich - commented; columns will be created in one of the NoteBooks 
        #self.df["salinity_availability"] = 0 # default
        #self.df["ctd_temp_availability"] = 0 # default
        # 20210511 - ED - Norwich >
        self.df["sample_temp_flag"] = pd.Series([]) # flag for sample temperature recorded during analysis. Must be set manually. 
        self.df["sampling_date"] = pd.Series([])
        self.df["sample_temperature_alternative"] = pd.Series([])
        
        
        # Populate the crm_batch column for CRMs only (Cast = batch for CRMs)
        self.df["crm_batch"].loc[self.df["Station"] == 8888] = self.df["Cast"]
        print("Using Cast number for CRMs as crm_batch input.")
        
        # Set DIC flag and TA flag to 0 if -999
        self.df["DIC flag"].loc[self.df["DIC flag"] == -999] = 0
        self.df["TA flag"].loc[self.df["TA flag"] == -999] = 0
        
        # Convert Date and Time format, create datetime column 
        self.df["Date (analysis)"] = pd.to_datetime(self.df["Date (analysis)"], format = "%Y%m%d" )
        self.df["analysis_datetime"] = pd.to_datetime(self.df["Date (analysis)"].astype(str) + " " + self.df["Time (analysis)"].astype(str).str.zfill(6), 
               format = "%Y-%m-%d %H%M%S")
        
        
        # Add all columns as instances to the class 
        # ---------------------------------------------------------------------
        # Set instances for all previously existing variables/columns
        self.sample = self.df["Sample Name"]
        self.comments = self.df["Comments"]
        self.dicflag = self.df["DIC flag"]
        self.taflag = self.df["TA flag"]
        self.analysistype = self.df["Analysis Type"]
        self.analysisdate = self.df["Date (analysis)"]
        self.analysistime = self.df["Time (analysis)"]
        self.station = self.df["Station"] # 9999 for junk sea water, 8888 for CRM, 2222 for underway 
        self.cast = self.df["Cast"] #
        self.niskin = self.df["Niskin"] #
        self.replicate = self.df["Replicate (diff bottle)"]
        self.repeat = self.df["Repeat(same bottle)"]
        self.numrepeats = self.df["No. Of Repeats"]
        self.depth = self.df["Depth (m)"]
        self.ctdtemp = self.df["CTD Temp (C)"] # Will be changed separately
        self.salinity = self.df["Salinity"] # Will be changed separately
        self.phosphate = self.df["Phosphate (uMol/Kg)"] # Will be changed separately
        self.silicate = self.df["Silicate (uMol/Kg)"] # Will be changed separately
        self.nitrate = self.df["Nitrate (uMol/Kg)"] # Will be changed separately
        self.dicpipettevol = self.df["DIC pipette volume (ml)"] # Make sure volume is calibrated
        self.dicep = self.df["DIC set end point (cpm)"]
        self.dicpippettetemp = self.df["DIC pipette temp sample (C)"]
        self.dicbackground = self.df["DIC background (cpm)"]
        self.dictime = self.df["DIC total time (s)"]
        self.diccounts = self.df["DIC total counts"]
        self.diccellnum = self.df["DIC cell Number (sequential)"] # CHECK THIS
        self.sampletemp = self.df["Sample Bottle temp (C)"]
        self.tapipettevol = self.df["TA pipette vol (ml)"] # Make sure volume is calibrated
        self.tacelltemp = self.df["TA cell temp (C)"]
        self.taacidbatch = self.df["TA Acid Batch"] 
        self.taacidconc = self.df["TA Acid conc (Mol)"]
        self.taaciddensity = self.df["TA Acid density (Kg/l)"]
        # < ED - 20210511 - Norwich - added
        self.salinity_availability = self.df["salinity_availability"] 
        self.ctd_temp_availability  = self.df["ctd_temp_availability"] 
        # ED - 20210511 - Norwich > 
        
        # <ED - 20210518 - Norwich - added
        # Check if these columns exist and insert if not
        check_columns = ["sampling_date"]
        for c in check_columns: 
            if c in self.df_ta.columns:
                print(c + " found in the columns of the DataFrame. Make sure it's correctly populated, i.e. not empty.")
            else: 
                print(c + " is not found in the columns of the DataFrame. Creating it here, but will remain empty.")
                print("Needs to be populated before proceeding with further calculations within this class!" )
                self.df_ta[c] = pd.Series([])
        # ED - 20210518 - Norwich > 
    
        
        # ---------------------------------------------------------------------
        # Set instances for all added variables/columns
        self.crm_batch = self.df["crm_batch"]
        # < ED - 20210511 - Norwich - commented, moved up 
        #self.salinity_availability = self.df["salinity_availability"]
        #self.ctd_temp_availability = self.df["ctd_temp_availability"]
        # ED - 20210511 - Norwich > 
        #self.sampling_date = self.df["sampling_date"] # < 20210407 - ED - UEA - commented > 
        self.analysis_datetime = self.df["analysis_datetime"]
        self.sample_temperature_alternative = self.df["sample_temperature_alternative"]
        self.dic_certified_umolkg = self.df["DIC_certified_umolkg"]# < ED - 20210514 - Norwich - added > 
        self.ta_certified_umolkg = self.df["TA_certified_umolkg"]# < ED - 20210514 - Norwich - added > 
        self.sampling_date = self.df["sampling_date"] # < ED - 20210518 - Norwich - added >
        
        # Print warnings/reminders
        print("Set an alternative sample temperature in sample_temperature_alternative column, if relevant!")
        
        
        #######################################################################################################################
    def add_dickson(self, dickson_location = "MOSAiC_pkg/crm_batch_lookup.csv"):
        """
        Looks up the certified concentrations for the relevant CRM batch. 
        CRM batch is specified in DataFrame of self as crm_batch (pd.Series)
        
        args: dickson_location  = string with path to csv file containing the CRM certified values in a particular format. Path is relative to location of current py file. 
        returns: no explicit returns, but does insert the certified concentrations into the rows of the DataFrame of self for CRM runs. 
        
        """
        print("Location of CRM Dickson info: " + dickson_location + ". Change if necessary.")
        
        # Find all unique batch numbers
        for b in self.df["crm_batch"][self.df["Station"] == 8888].unique():
            # Insert the Dickson Certified Values for the CRMs
            self.crm_dickson = crmbatch.findcrmbatch(b, dickson_location) # < 20210407 - ED - UEA - b[0] to b >
                
            # For all CRM runs with this batch number, replace salinity, phosphate, silicate, and nitrate values with those of the CRM batch
            cond = (self.station == 8888) & (self.crm_batch == b)# < 20210407 - ED - UEA - b[0] to b >
            self.salinity[cond] = self.crm_dickson["salinity"]
            self.phosphate[cond] = self.crm_dickson["phos"]
            self.silicate[cond] = self.crm_dickson["sil"]
            self.nitrate[cond] = self.crm_dickson["nitrate"]
            
            # < ED - 20210512 - Norwich - added
            # Also populate the certified DIC and TA for the CRMs in the columns already created in the initialisation of the Class
            self.dic_certified_umolkg[cond] = self.crm_dickson["dic"]
            self.ta_certified_umolkg[cond] = self.crm_dickson["alk"]
            # ED - 20210512 - Norwich > 
        
        # Synchronise with the DataFrame of self
        self.df["Salinity"] = self.salinity
        self.df["Phosphate (uMol/Kg)"] = self.phosphate
        self.df["Silicate (uMol/Kg)"] = self.silicate
        self.df["Nitrate (uMol/Kg)"] = self.nitrate
        
        # < ED - 20210512 - Norwich - added
        self.df["DIC_certified_umolkg"] = self.dic_certified_umolkg
        self.df["TA_certified_umolkg"] = self.ta_certified_umolkg
        # ED - 20210512 - Norwich > 
        
        print("Added the CRM information for salinity, phosphate, silicate, and nitrate.")

        #######################################################################################################################
    def calculate_dic(self, use_alternative_temp = False):
        """
        """
        # Calculate the DIC concentration for all samples using the Dickson CRMs
        
        # Set the temperature you want to use as sample temperature.
        sampletemp = self.sampletemp
        print("Using sample temperature as temperature. You can change temperature to use in MOSAiC_calc_DIC.py")
        if use_alternative_temp == True: 
            sampletemp = self.sample_temperature_alternative
            print("Using sample_temperature_alternative as temperature. ")
        
        
        
        # Calculate the net DIC counts by removing the background DIC counts
        # Counts - (blanks [cts/min]/60)*seconds [sec])
        self.netdiccounts = self.diccounts - ((self.dicbackground/60)*self.dictime)
        self.df["netdiccounts"] = self.netdiccounts
        
        #######################################################################################################################
        # Volume correction for addition of mercuric chloride 
        # 100uL per 0.5L or 50uL per 0.25L of saturated HgCl2 solution
        # 0.5001/0.500 = 0.25005/0.250
        self.diccounts_volcorr = self.netdiccounts*(0.5001/0.5000) # < 20210311 - ED - Norwich - changed diccounts --> netdiccounts > 
        self.df["diccounts_volcorr"] = self.diccounts_volcorr
        # --> not used in calculations
        
        #######################################################################################################################
        # Calculate density of seawater in DIC pipette [kg/m3] or [g/1000ml]
        # self.dic_sw_density = gsw_rho_t_exact(self.salinity, self.sampletemp, 10.1325)
        # <<<<<<< !!!!!!!!!! >>>>>>>> THIS IS HOW OLLIE LEGGE CALCULATED IT, BUT NO REFERENCE. DIFFERENT FROM P 2 OF 4 SOP 12
        # (MILLERO AND POISSON, 1981). This is what Ruth uses too. Check. 
        
        
        RohSMOW = (999.842594 + (6.793952e-2*sampletemp) - 
                   (9.09529e-3*sampletemp**2) + (1.001685e-4*sampletemp**3) -
                   (1.120083e-6*sampletemp**4) + (6.536332e-9*sampletemp**5))
                
        A = (8.24493e-1 - (4.0899e-3*sampletemp) + (7.6438e-5*(sampletemp**2)) -
             (8.2467e-7*(sampletemp**3)) + (5.3875e-9*(sampletemp**4)))
        
        B = -5.72466e-3 + (1.0227e-4*sampletemp) - (1.6546e-6*(sampletemp**2))
        
        C = 4.8314e-4
        
        self.dic_sw_density = RohSMOW + A*self.salinity + B*self.salinity**1.5 + C*self.salinity**2
        self.df["dic_sw_density"] = self.dic_sw_density
                
        #######################################################################################################################
        # Calculate mass of pipette volume [g/pipette]
        # [g / pip] = SWDensity [g / 1000 ml] / 1000 * Pipette Volume [ml]
        self.dicpipmass = (self.dic_sw_density/1000)*self.dicpipettevol
        self.df["dicpipmass"] = self.dicpipmass
        
        # To correct for glass expansion at measurement temperature, correct soelf.dicpipmass with aa3 
        # aa3 = Glass expansion factor relative to calibration at 25C
        # aa3 = 1 + 1e-5 * (self.sampletemp - 25)
        # self.dicpipmass = (self.dicswdensity/1000)*self.dicpipettevol*aa3
        # <<<<<< !!!!!!!!!!!!! >>>>>> CHECK THIS: why has this not been applied? 
        
        #######################################################################################################################        
        # Calculate umol DIC in the pipette volume titrated [umol/pipette] for each CRM, using CRM batch values to get umol per pipette
        # [umol / pip] = CRM DIC [umol/kg] / 1000 * pip mass [g/pip]
        self.dic_crm_pip_umol = (self.crm_dickson["dic"]/1000) * self.dicpipmass[self.sample.str.contains("8888")] # maintains same index!
        self.df["dic_crm_pip_umol"] = self.dic_crm_pip_umol
        
        #######################################################################################################################
        # Calculate counts per mol in pipette volume titrated [CTS/umol] for each CRM, using CRM values to get counts per umol
        #[CTS / umol] = pip CTS [CTS / pip] / pip umol [umol / pip]
        # < 20210311 - ED - Norwich - changed diccounts --> netdiccounts 
        self.dic_crm_counts_per_umol = self.netdiccounts[self.sample.str.contains("8888")]/self.dic_crm_pip_umol[self.sample.str.contains("8888")] 
        # 20210311 - ED - Norwich > 
        self.df["dic_crm_counts_per_umol"] = self.dic_crm_counts_per_umol
        
        #######################################################################################################################
        # Extract the sampling date from the sample name
        
        # Following was causing problems when number of digits didn't fit, so I changed it (see below)
        # self.date_all = pd.to_datetime(self.sample.str.extract(pat="(\d{8})", expand=False) , format="%Y%m%d") # expand=False creates a Series instead of DataFrame
        # self.sampling_date = self.date_all[self.sample.str.contains("9999|8888")==False] # 9999 = junk, 8888 = CRM; could also have used the "Station" variable for this... 
        #self.date_all = pd.to_datetime(self.cast[self.station != 8888] , format="%Y%m%d")
        #self.sampling_date = self.date_all[(self.cast != 9999) or (self.cast != 8888)]
        
        # < 20210407 - ED - UEA - commented 
        #self.sampling_date = pd.to_datetime(self.cast[(self.station != 8888)&(self.station != 0000)] , format="%Y%m%d")
        #self.df["sampling_date"] = self.sampling_date
        # 20210407 - ED - UEA > 

        #######################################################################################################################
        # Calibrate DIC samples / Calculate sample DIC [umol/kg]
        
        # PRELIM VERSION: NO QC ON CRMs DONE YET! 
        # Calculate mean and standard deviation of counts/umol for CRMs
        # "perday" -> DIC concentration calculated using CRM averages per analysis day i.e. cell. 
        # "all" -> DIC concentration calculated using CRM average of ALL CRM runs in the summary file
        # Better would be to use the cell# instead of analysis day column in case analyses cover >1 day !!!! 
        
        # Remember: numpy ignores the nans automatically
        self.dic_crm_counts_per_umol_mean_allcells = np.nanmean(self.dic_crm_counts_per_umol)
        self.dic_crm_counts_per_umol_std_allcells = np.nanstd(self.dic_crm_counts_per_umol)
        
        # WARNING: these function are also applied to the CRMs. The latter are thus divided by their own collective mean. 
        # calculating calibrated moles titrated per pipette [umol]
        # [umol] = counts per pipette [CTS/pip] / counts per umol [CTS/umol]
        
        # DIC concentration based on average of ALL CRM runs in summary file
        self.dic_sample_pip_umol_allcells = self.netdiccounts/self.dic_crm_counts_per_umol_mean_allcells
        
        # Calculate sample concentration [umol/l]
        self.dic_sample_umol_per_l_allcells = (self.dic_sample_pip_umol_allcells/self.dicpipettevol)*1000 # dicpipettevol is in mL

        # Calculate sample concentration [umol/kg]
        self.dic_sample_umol_per_kg_allcells = (self.dic_sample_umol_per_l_allcells/self.dic_sw_density)*1000 # dic_sw_density is in [kg/m3] or [g/1000ml]
        
        # Add to given data frame self.df
        self.df["dic_sample_pip_umol_allcells"] = self.dic_sample_pip_umol_allcells
        self.df["dic_sample_umol_per_l_allcells"] = self.dic_sample_umol_per_l_allcells
        self.df["dic_sample_umol_per_kg_allcells"] = self.dic_sample_umol_per_kg_allcells
        
        ###############################################################
        
        # Create a list of all unique analysis days (i.e. cells)
        index_days = self.dic_crm_counts_per_umol.groupby([self.diccellnum]).size().index #Int64Index
        print("There are "+str(len(index_days))+" unique analysis days (i.e. DIC cells) during which samples were run.")

        # For each of the unique analysis days (i.e. cells), determine the mean of the CRMs' measured DIC counts per umol
        crm_mean_per_cell = self.dic_crm_counts_per_umol.groupby([self.diccellnum]).transform("mean").unique()
        crm_std_per_cell = self.dic_crm_counts_per_umol.groupby([self.diccellnum]).transform("std").unique()

        # Create a dictionary that links the unique analysis days (cells) with their respective mean for CRM DIC counts per umol
        crm_mean_dict = {}
        crm_std_dict = {}
        for i in np.arange(len(index_days)):
            crm_mean_dict[index_days[i]] = crm_mean_per_cell[i]
            crm_std_dict[index_days[i]] = crm_std_per_cell[i]

        # For DIC concentration calculation purposes, we need a Series with values that are the mean CRM values per analysis day
        # NOTE: This process can probably be made more efficient, but am leaving it like this for now. Get on with it. Lambda x?
        crm_mean_serieslist = [] # List created to append crm mean values to; will be turned into a Pandas Series and added to the data frame
        crm_std_serieslist = []
        for i in self.diccellnum:
            if i in crm_mean_dict:
                crm_mean_serieslist.append(crm_mean_dict[i])
                crm_std_serieslist.append(crm_std_dict[i])
        self.df["dic_crm_counts_per_umol_mean_percell"] = crm_mean_serieslist # Add list to data frame
        self.df["dic_crm_counts_per_umol_std_percell"] = crm_std_serieslist
        self.dic_crm_counts_per_umol_mean_percell = self.df["dic_crm_counts_per_umol_mean_percell"] # Add Series as an instance to Class
        self.dic_crm_counts_per_umol_std_percell = self.df["dic_crm_counts_per_umol_std_percell"]
        
        # Do exactly the same as above for "all", but now using the CRM average per analysis day (cell)
        self.dic_sample_pip_umol = self.netdiccounts/self.dic_crm_counts_per_umol_mean_percell
        self.dic_sample_umol_per_l = (self.dic_sample_pip_umol/self.dicpipettevol)*1000 # dicpipettevol is in mL
        self.dic_sample_umol_per_kg = (self.dic_sample_umol_per_l/self.dic_sw_density)*1000 # dic_sw_density is in [kg/m3] or [g/1000ml]
                
        # Add to given data frame self.df
        self.df["dic_sample_pip_umol"] = self.dic_sample_pip_umol
        self.df["dic_sample_umol_per_l"] = self.dic_sample_umol_per_l
        self.df["dic_sample_umol_per_kg"] = self.dic_sample_umol_per_kg
        
        
        self.mean_diff = self.dic_sample_umol_per_kg - self.dic_sample_umol_per_kg_allcells
    
    
    #######################################################################################################################
    def calculate_dic_driftcorr(self, tobe_driftcorr_list, time_frame): 
        """Apply drift corrections to given dates or periods"""
        # If you want to overwrite any previously calculated drift corrections, 
        # set the dic_driftcorr_flag to 0 and dic_crm_counts_per_umol_driftcorr to np.float("NaN") BEFORE applying this method! 

    
        # Check if the flag for DIC drift correction exists
        if "dic_driftcorr_flag" in self.df.columns:
            # If it exists, it means that this method has been applied before, in which case we're NOT going to overwrite what has been done previously! 
            print("Not the first time drift corrections have been done on this dataset!")
    
        else: 
            # Create a column (or overwrite the existing one) that contains the flag for whether samples have been drift corrected
            # 0 = not drift corrected, 1 = drift corrected. Default is set to 0. 
            self.df["dic_driftcorr_flag"] = 0 
    
        if "dic_crm_counts_per_umol_driftcorr" in self.df.columns:
            pass
        else:
            # Create a new column that will contain the DIC counts per umol that will have been drift corrected
            self.df["dic_crm_counts_per_umol_driftcorr"] = np.float("NaN")
    
        if not tobe_driftcorr_list: 
            print("List is empty, no drift corrections done.")
    
        else: 
        
            fig = plt.figure(figsize = (15, 10))
            x, y = 2, 3 # subplot dimensions
            z = 0
        
            try: 
                if time_frame == "single_days":
                    for corr_day in tobe_driftcorr_list:
                        z = z + 1

                        # Change the flag for DIC drift corrected from 0 to 1 
                        self.df["dic_driftcorr_flag"][self.df["Date (analysis)"] == corr_day] = 1

                        # Extract the relevant variables and convert them to an array and reshape so it can be used by the LinearRegression class
                        crm_time_array = self.df["Time (analysis)"][(self.df["Date (analysis)"] == corr_day) & (self.df["Station"] == 8888)].values.reshape(-1,1)
                        samples_time_array = self.df["Time (analysis)"][self.df["Date (analysis)"] == corr_day].values.reshape(-1,1)
                        crm_ctsperumol_array = self.df["dic_crm_counts_per_umol"][(self.df["Date (analysis)"] == corr_day) & (self.df["Station"] == 8888)].values.reshape(-1,1)
                        

                        # Use linear regression to interpolate the CRM counts/umol at the time of each sample
                        linear_regressor = LinearRegression()
                        linear_regressor.fit(crm_time_array, crm_ctsperumol_array)
                        Y_pred = linear_regressor.predict(crm_time_array)

                        # Create a list containing the predicted values; this is used for plotting afterwards. 
                        crm_drift_corr_plt = []

                        # Per sample, determine the counts/umol based on the linear regression
                        for i in samples_time_array:
                            crm_drift_corr = linear_regressor.predict(i.reshape(-1,1))[0][0]

                            # Save the drift corrected value
                            ix = self.df[(self.df["Date (analysis)"] == corr_day) & (self.df["Time (analysis)"] == i[0])].index
                            self.df.loc[ix[0], "dic_crm_counts_per_umol_driftcorr"] = crm_drift_corr

                            # Save the predicted value to the list (used for plotting only)
                            crm_drift_corr_plt.append(crm_drift_corr)


                        # Create a simple plot just to show what has been done 
                        ax = fig.add_subplot(x, y, z)

                        ax.scatter(crm_time_array, crm_ctsperumol_array)
                        ax.plot(crm_time_array, Y_pred, color="red")
                        ax.scatter(samples_time_array, crm_drift_corr_plt)
                        ax.set_xlabel("Time")
                        ax.set_ylabel("counts per umol")
                        ax.set_title(str(corr_day))

                elif time_frame == "multiple_days": 
                
                    z = 1
                
                    # tobe_driftcorr_list needs to have exaclty 2 elements i.e. dates between which linear interpolation/drift correction will be applied
                    date_min = tobe_driftcorr_list[0]
                    date_max = tobe_driftcorr_list[-1]
                
                    # Flag the days that are drift corrected
                    self.df["dic_driftcorr_flag"][(self.df["Date (analysis)"]>date_min) & (self.df["Date (analysis)"] < date_max)] = 1 
                
                    # Extract the relevant variables and convert to array + reshape for use in LinearRegression
                    crm_time_array = self.df["Time (analysis)"][(self.df["Date (analysis)"] > date_min) & (self.df["Date (analysis)"] < date_max) & (self.df["Station"] == 8888)].values.reshape(-1,1)
                    samples_time_array = self.df["Time (analysis)"][(self.df["Date (analysis)"] > date_min) & (self.df["Date (analysis)"] < date_max) & (self.df["Station"] != 8888)].values.reshape(-1,1)
                    crm_ctsperumol_array = self.df["dic_crm_counts_per_umol"][(self.df["Date (analysis)"] > date_min) & (self.df["Date (analysis)"] < date_max) & (self.df["Station"] == 8888)].values.reshape(-1,1)
                
                    # Use linear regression to interpolate CRM counts/umol at time of each sample
                    linear_regressor = LinearRegression()
                    linear_regressor.fit(crm_time_array, crm_ctsperumol_array)
                    Y_pred = linear_regressor.predict(crm_time_array)
                
                    # Create a list containing the predicted values; this is used for plotting afterwards. 
                    crm_drift_corr_plt = []
                
                    # Per sample, determine the counts/umol based on the linear regression
                    for i in self.df[(self.df["Date (analysis)"]>date_min) & (self.df["Date (analysis)"] < date_max) & (self.df["Station"] != 8888)].index: 
                    
                        crm_drift_corr = linear_regressor.predict(self.df["Time (analysis)"][i].reshape(-1,1))[0][0]
                    
                        # Save drift corrected value 
                        self.df.loc[i, "crm_dic_counts_per_umol_driftcorr"] = crm_drift_corr
                    
                        # Save the predicted value to the list (used for plotting only)
                        crm_drift_corr_plt.append(crm_drift_corr)
                    
                    # Create a simple plot just to show what has been done 
                    ax = fig.add_subplot(x, y, z)

                    ax.scatter(crm_time_array, crm_ctsperumol_array)
                    ax.plot(crm_time_array, Y_pred, color="red")
                    ax.scatter(samples_time_array, crm_drift_corr_plt)
                    ax.set_xlabel("Time")
                    ax.set_ylabel("counts per umol")
                    ax.set_title(str(date_min) + " to " + str(date_max))
            
                else:
                    raise IncorrectArgumentError 
        
                # After all the necessary drift corrections have been completed, re-calculate the drift-corrected DIC concentrations 
                # We do not need to check if the new columns already exist. If they don't, they'll automatically be created. If they do, then values are recalculated. 
                # But previously calculated drift corrected values are not overwritten above, so they should still end up with the same calculated value (umol/kg)
                self.df["dic_sample_pip_umol_driftcorr"] = self.df["netdiccounts"]/self.df["dic_crm_counts_per_umol_driftcorr"]
                self.df["dic_sample_umol_per_l_driftcorr"] = (self.df["dic_sample_pip_umol_driftcorr"]/self.df["DIC pipette volume (ml)"])*1000
                self.df["dic_sample_umol_per_kg_driftcorr"] = (self.df["dic_sample_umol_per_l_driftcorr"]/self.df["dic_sw_density"])*1000
        
            except IncorrectArgumentError: 
                print("IncorrectArguentError: time_frame argument can only be single_days or multiple_days. Correct and try again.")
        
        
        
        
        #######################################################################################################################    
    def extract_crms(self):
        # Extract the CRM data into a separate data frame
        self.crm = self.df[self.df["Sample Name"].str.contains("8888")] # this is a dataframe for just the CRM measurements
        self.crm.sample = self.crm["Sample Name"]
        self.crm.comments = self.crm["Comments"]
        self.crm.dicflag = self.crm["DIC flag"]
        self.crm.taflag = self.crm["TA flag"]
        self.crm.analysistype = self.crm["Analysis Type"]
        self.crm.analysisdate = self.crm["Date (analysis)"]
        self.crm.analysistime = self.crm["Time (analysis)"]
        self.crm.station = self.crm["Station"] #
        self.crm.cast = self.crm["Cast"] #
        self.crm.niskin = self.crm["Niskin"] #
        self.crm.replicate = self.crm["Replicate (diff bottle)"]
        self.crm.repeat = self.crm["Repeat(same bottle)"]
        self.crm.numrepeats = self.crm["No. Of Repeats"]
        self.crm.depth = self.crm["Depth (m)"]
        self.crm.ctdtemp = self.crm["CTD Temp (C)"]
        self.crm.salinity = self.crm["Salinity"]
        self.crm.phosphate = self.crm["Phosphate (uMol/Kg)"]
        self.crm.silicate = self.crm["Silicate (uMol/Kg)"]
        self.crm.nitrate = self.crm["Nitrate (uMol/Kg)"]
        self.crm.dicpipettevol = self.crm["DIC pipette volume (ml)"]
        self.crm.dicep = self.crm["DIC set end point (cpm)"]
        self.crm.dicpippettetemp = self.crm["DIC pipette temp sample (C)"]
        self.crm.dicbackground = self.crm["DIC background (cpm)"]
        self.crm.dictime = self.crm["DIC total time (s)"]
        self.crm.diccounts = self.crm["DIC total counts"]
        self.crm.diccellnum = self.crm["DIC cell Number (sequential)"]
        self.crm.sampletemp = self.crm["Sample Bottle temp (C)"]
        self.crm.tapipettevol = self.crm["TA pipette vol (ml)"]
        self.crm.tacelltemp = self.crm["TA cell temp (C)"]
        self.crm.taacidbatch = self.crm["TA Acid Batch"] 
        self.crm.taacidconc = self.crm["TA Acid conc (Mol)"]
        self.crm.taaciddensity = self.crm["TA Acid density (Kg/l)"]
        
        
        self.crm.dic_sample_pip_umol_allcells = self.crm["dic_sample_pip_umol_allcells"]
        self.crm.dic_sample_umol_per_l_allcells = self.crm["dic_sample_umol_per_l_allcells"]
        self.crm.dic_sample_umol_per_kg_allcells = self.crm["dic_sample_umol_per_kg_allcells"]
        self.crm.dic_sample_pip_umol = self.crm["dic_sample_pip_umol"]
        self.crm.dic_sample_umol_per_l = self.crm["dic_sample_umol_per_l"]
        self.crm.dic_sample_umol_per_kg = self.crm["dic_sample_umol_per_kg"]
        # self.crm.date_all = self.crm["date_all"]
        # self.crm.sampling_date = self.crm["sampling_date"] # < 20210407 - ED - UEA - commented > 
        
        
    #######################################################################################################################
    # Quality check the CRMs
    # This method will plot either the counts or the calculated concentrations for the CRMs per analysis day
    def qc_crms(self, use="counts"):
        """
        TO DO: 
            - change analysis days to cell number 
            - check PS117 module for updates 
        """
        
        # WARNING: this function works when analyses have been done within 24 hours, i.e. on the same date (or at least labelled like that)
        # Better would be to distinguish analysis days on cell number, instead of date. Maybe change this at a later point! 
        
        # 'use' argument is used to specify whether to plot the CRM "counts" or the "concentration". Default is counts. 
        
        # Determine the number of analysis days according to cell number
        analysis_days = self.crm.analysisdate.unique()
        x = 4 # Number of rows of subplots
        y = 4 # Number of columns of subplots

        # Create empty list for analysis days that will need drift checking later
        check_drift = []

        # Start figure
        # Plot the CRMs measured for each subplot vs analysis time
        fig = plt.figure(figsize = (20, 20))
        

        # For each cell number, generate a subplot
        for i in np.arange(len(analysis_days)): 
            if use == "counts":
                df_mini = self.crm.loc[self.crm.analysisdate == analysis_days[i]].dropna(subset = ["netdiccounts"]) # Create a mini data frame with just the CRM data per analysis day ...
                # AND remove any rows where the netdiccounts for CRMs is NaN!
                crm_dic = df_mini["netdiccounts"]
                ylim_min_text = np.min(self.crm.netdiccounts)
                ylim_max_text = np.max(self.crm.netdiccounts)
                ylim_min = ylim_min_text-100
                ylim_max = ylim_max_text+100
            elif use == "concentration":
                df_mini = self.crm.loc[self.crm.analysisdate == analysis_days[i]].dropna(subset = ["dic_sample_umol_per_kg"]) # Create a mini data frame with just the CRM data per analysis day ...
                # AND remove any rows where the netdiccounts for CRMs is NaN!
                crm_dic = df_mini["dic_sample_umol_per_kg"]
                ylim_min_text = np.min(self.crm.dic_sample_umol_per_kg)
                ylim_max_text = np.max(self.crm.dic_sample_umol_per_kg)
                ylim_min = ylim_min_text-1
                ylim_max = ylim_max_text+1
            else:
                print("Set argument of qc_crms() function to either 'counts' or to 'concentration'")
            
            crm_time = df_mini["Time (analysis)"]
            z = i+1
            ax = fig.add_subplot(x, y, z)
            ax.scatter(crm_time, crm_dic, marker = "o", c = "k", edgecolors = "b")

            # Perform linear regression using the LinearRegression Class from sklearn.linear_model
            # I used this website for inspiration: https://towardsdatascience.com/linear-regression-in-6-lines-of-python-5e1d0cd05b8d
            crm_dic_array = crm_dic.values.reshape(-1, 1) # values converts series into a numpy array, -1 indicates to calculate the dimension of rows, but to have one column
            crm_time_array = crm_time.values.reshape(-1, 1)

            linear_regressor = LinearRegression() # create object for class
            linear_regressor.fit(crm_time_array, crm_dic_array) # linear regression
            Y_pred = linear_regressor.predict(crm_time_array) # make the predictions

            # Plot the linear regression line
            ax.plot(crm_time_array, Y_pred, c = "r")

            # Add titles
            title = str(analysis_days[i])
            ax.set_title(title)
            
            # Calculate the standard deviation of all CRM netdiccounts values per analysis day
            std_crms_day = np.std(crm_dic_array)
            text2 = "std="+str(std_crms_day)

            # Calculate the min and max of the regression line and add it to the graph
            Y_pred_range = np.max(Y_pred) - np.min(Y_pred)
            text1 = "diff="+str(Y_pred_range)+"cts"
            
            # Add the limitations for the y axis
            ax.set_ylim(ylim_min, ylim_max)
            
            if use == "counts":
                if Y_pred_range > 100:
                    textcolour = "r"
                    check_drift.append(analysis_days[i])
                    printtext = "The following analysis days have drift > 100 counts during the day: "
                else:
                    textcolour = "k"
            elif use == "concentration":
                if Y_pred_range > 3:
                    textcolour = "r"
                    check_drift.append(analysis_days[i])
                    printtext = "The following analysis days have drift > 3 umol/kg during the day: "
                else:
                    textcolour = "k"

            ax.text(np.min(crm_time_array), ylim_max_text, text1, color = textcolour)
            ax.text(np.min(crm_time_array), ylim_min_text, text2)

        print("X axis: analysis time -- Y axis: Net DIC counts (background counts corrected) -- Title: analysis date in yyymmdd -- diff = max-min of lin reg")
        if len(check_drift) > 1:
            print(printtext)
            for i in check_drift:
                print(i)
        else: 
            print("No drift to remark.")
        
        

    #######################################################################################################################
    def qc_conc_profiles(self): 
        """
        TO DO: 
            - use PS117 for base
            - use single values, no means 
        """

    #######################################################################################################################
    
    def qc_conc_timeseries(self):
        """
        
        Plot the DIC concentrations calculated within this class for all samples
        Plotted in an interactive plot using plotly so it's easy for the user to identify samples
        
        TO DO: 
            - adjust sampling date -> create? 
            - clean up 
            - Use PS117 module for updates
        """
        
        # STILL NEED TO CORRECT THIS FUNCTION!!! No sampling_date!! 20210407 > 
        
#        << 20200120 - ED - UEA - Using the columns of the dataframe instead of the instances of the class, because it was causing samples to be plotted against incorrect sampling dates... 
        sdates= self.df["sampling_date"][self.df["Sample Name"].str.contains("9999|8888")==False]
        single_measurements= self.df["dic_sample_umol_per_kg"][self.df["Sample Name"].str.contains("9999|8888")==False]
        
        single_measurements_allcrm= self.df["dic_sample_umol_per_kg_allcells"][self.df["Sample Name"].str.contains("9999|8888")==False]
        
        single_text = np.array(self.df["Sample Name"][self.df["Sample Name"].str.contains("9999|8888")==False])
#        20200120 - ED - UEA >> 
        
        
        # Set lay out properties
        layout = go.Layout(hovermode="closest", # to make sure hover text appears for points that share an x value
                           title = "Check DIC Concentrations for Time Series", 
        #                    xaxis_rangeslider_visible = True
                          )

        # Create figure
        fig = go.Figure(layout=layout)
        
        ###
        fig.add_trace(go.Scatter(
            x = sdates, 
            y = single_measurements, 
            mode = "markers", 
            name = "single measurements - crm per day", 
            text = single_text
        ))
        
        ###
        fig.add_trace(go.Scatter(
            x = sdates, 
            y = single_measurements_allcrm, 
            mode = "markers", 
            name = "single measurements - crm total", 
            text = single_text
        ))
        

        # Set axes titles
        fig.update_xaxes(title_text = "Time")
        fig.update_yaxes(title_text = "DIC [umol kg-1]")

        fig.show()
    
    ####################################################################################################################### 
    def dic_x_chart(self):
        """
        X Chart according to SOP22 in Dickson et al. 
        """
        # Set lay out properties
        layout = go.Layout(hovermode="closest", # to make sure hover text appears for points that share an x value
                           title = "$\overline{X}$ Chart DIC Concentrations for CRMs", 
                           showlegend=False, 
                          )
        
        sequence = np.arange(len(self.crm.index))
        crm_conc = self.crm.dic_sample_umol_per_kg
        
        x = np.mean(self.crm.dic_sample_umol_per_kg)
        s = np.std(self.crm.dic_sample_umol_per_kg)
        ucl = x + 3*s
        uwl = x + 2*s
        lwl = x - 2*s
        lcl = x - 3*s
        
        # Create figure
        fig = go.Figure(layout=layout)
        
        fig.layout.update(
                shapes = [
                        # x
                        go.layout.Shape(
                        type="line",
                            x0=np.min(sequence),
                            y0=x,
                            x1=np.max(sequence),
                            y1=x,
                            line=dict(
                                color="RoyalBlue",
                                width=3,                                
                                ),
                            ),
                        
                        # ucl
                        go.layout.Shape(
                            type="line",
                            x0=np.min(sequence),
                            y0=ucl,
                            x1=np.max(sequence),
                            y1=ucl,
                            line=dict(
                                color="LightSeaGreen",
                                width=3,
                                dash="dot",
                                ),
                            ),
                        # uwl
                        go.layout.Shape(
                                type="line",
                                x0=np.min(sequence),
                                y0=uwl,
                                x1=np.max(sequence),
                                y1=uwl,
                                line=dict(
                                    color="MediumPurple",
                                    width=3,
                                    dash="dashdot",
                                    ),
                                ),
                        #lwl
                        go.layout.Shape(
                                type="line",
                                x0=np.min(sequence),
                                y0=lwl,
                                x1=np.max(sequence),
                                y1=lwl,
                                line=dict(
                                    color="MediumPurple",
                                    width=3,
                                    dash="dashdot",
                                    ),
                                ),
                                        
                        # lcl
                        go.layout.Shape(
                                type="line",
                                x0=np.min(sequence),
                                y0=lcl,
                                x1=np.max(sequence),
                                y1=lcl,
                                line=dict(
                                    color="LightSeaGreen",
                                    width=3,
                                    dash="dot",
                                    ),
                                ),          
                    ]
                )
        
        
        ###
        fig.add_trace(go.Scatter(
            x = sequence, 
            y = crm_conc, 
            mode = "markers", 
            hovertext = np.array(self.crm.sample)
        ))
        
        fig.add_trace(go.Scatter(
            x=[1.5, 1.5, 1.5, 1.5, 1.5],
            y=[ucl+1, uwl+1, x+1, lwl+1, lcl+1],
            text=["UCL", "UWL", "$\overline{X}$", "LWL", "LCL"],
            mode="text",
        ))
        

        # Set axes titles
        fig.update_xaxes(title_text = "Sequence")
        fig.update_yaxes(title_text = "CRM DIC [umol kg-1]")

        fig.show()
    
    #######################################################################################################################
    def dic_r_chart(self):
        """
        R chart according to SOP22 in Dickson et al.    
        """
        
        # Set lay out properties
        layout = go.Layout(hovermode="closest", # to make sure hover text appears for points that share an x value
                           title = "$\overline{R}$ Chart DIC Concentrations for CRMs", 
                           showlegend=False, 
                          )
        
        
        
        crm_rdiff = abs(self.crm.dic_sample_umol_per_kg.groupby(self.crm.niskin).diff())
        sequence = np.arange(len(crm_rdiff))
        
        r = np.mean(crm_rdiff)
        s = np.std(crm_rdiff)
        ucl = r + 3*s
        uwl = r + 2*s
        
        # Create figure
        fig = go.Figure(layout=layout)
        
        fig.layout.update(
                shapes = [
                        # r
                        go.layout.Shape(
                        type="line",
                            x0=np.min(sequence),
                            y0=r,
                            x1=np.max(sequence),
                            y1=r,
                            line=dict(
                                color="RoyalBlue",
                                width=3,
                                #dash="dashdot",
                                ),
                            ),
                        
                        # ucl
                        go.layout.Shape(
                            type="line",
                            x0=np.min(sequence),
                            y0=ucl,
                            x1=np.max(sequence),
                            y1=ucl,
                            line=dict(
                                color="LightSeaGreen",
                                width=3,
                                dash="dot",
                                ),
                            ),
                        # uwl
                        go.layout.Shape(
                                type="line",
                                x0=np.min(sequence),
                                y0=uwl,
                                x1=np.max(sequence),
                                y1=uwl,
                                line=dict(
                                    color="MediumPurple",
                                    width=3,
                                    dash="dashdot",
                                    ),
                                ),          
                    ]
                )
        
        
        ###
        fig.add_trace(go.Scatter(
            x = sequence, 
            y = crm_rdiff, 
            mode = "markers", 
            hovertext = np.array(self.crm.sample)
        ))
        
        fig.add_trace(go.Scatter(
            x=[1.5, 1.5, 1.5],
            y=[ucl+0.5, uwl+0.5, r+0.5],
            text=["UCL", "UWL", "$\overline{R}$]"],
            mode="text"
        ))
        

        # Set axes titles
        fig.update_xaxes(title_text = "Sequence")
        fig.update_yaxes(title_text = "CRM R [umol kg-1]")

        fig.show()
        

print("Imported CalcDIC class.")