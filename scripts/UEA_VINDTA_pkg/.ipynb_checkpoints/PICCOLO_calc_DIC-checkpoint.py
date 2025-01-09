#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 11:14:21 2021

@author: Elise Droste (e.droste@uea.ac.uk)

Conventions for this script to work: 
    
    For CRMs and junks:
    - Station = 9999 for junk seawater runs, 8888 for CRM runs
    - Cast = yyyymmdd of analysis for junk seawater runs, CRM batch number for CRM runs 
    - Niskin = junk number of the day for junks, CRM bottle numbre for CRMs 
    - Depth = 0 for junks and CRMs
    - Repeat (same bottle) = automatic repeat run (set by software)
    - Replicate (different bottle) = 1 for junk and CRM runs, unique bottle number for samples
    
    For samples: 
    - Station = station
    - Cast = cast
    - Niksin = niskin
    - Depth = depth
    - Repeat = repeat from the same bottle (automatic)
    - Replicate = replicate/duplicate sample
    
Set column names in the column_names.py style script, as well as identifiers for CRM and junk runs. 
Column_names.py also contains some styling options. 

"""


print("Hello! Before running the CalcDIC Class, make sure the filenames of the runs are correct.")

# Import required built-in packages
# Consider being more specific for packages e.g. from numpy import max as np_max

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from . import crmbatch

# import plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots

# Import list of columns names 
from . import column_names


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
    Class to determine dissolved inorganic carbon from raw VINDTA data. 
    
    Change column names in column_names.py if required. 
    
    
    """
    
    #######################################################################################################################
    def __init__(self, df):
        """Initialise the Class by setting the fields as instances. Also do some cleaning up. """
        
        # ---------------------------------------------------------------------
        # Some cleaning up 
        # Copy data frame and remove the junk sea water sample runs + reset the index
        self.df = df[df[column_names.station]!=column_names.junkid].reset_index(drop=True)
        print("All junk samples (Station = "+ str(column_names.junkid)+ ") have been removed.")
        
        # Add some new columns (others will be created along the way) 
        self.df[column_names.crm_batch] = pd.Series([])
        self.df[column_names.dic_certified_umolkg] = pd.Series([]) # ED - 20210512 - Norwich - added > 
        self.df[column_names.ta_certified_umolkg] = pd.Series([]) # ED - 20210512 - Norwich - added > 
        # < 20210511 - ED - Norwich - commented; columns will be created in one of the NoteBooks 
        #self.df["salinity_availability"] = 0 # default
        #self.df["ctd_temp_availability"] = 0 # default
        # 20210511 - ED - Norwich >
        
        # < ED - 20210922 - UEA - added sampling_date to this list
        # < ED - 20210528 - Norwich - added as back up 
        check_columns = [column_names.salinity_availability, column_names.ctd_temp_availability, column_names.sampling_date, column_names.sample_temp_flag] # < ED - 20211111 - UEA - added sample_temp_flag > 
        for i in check_columns: 
            if i not in self.df.columns:
                self.df[i] = pd.Series([]) # default
                print("Creating " + i + " as a column")
            else: 
                print(i + " exists as column")
        # ED - 20210528 - Norwich > 
        # ED - 20210922 - UEA > 
        
        # self.df[column_names.sample_temp_flag] = pd.Series([]) # flag for sample temperature recorded during analysis. Must be set manually. # < ED - 20211111 - UEA - commented > 
        # self.df[column_names.sampling_date] = pd.Series([]) # < ED - 20210922 - UEA - commented > 
        self.df[column_names.sample_temperature_alternative] = pd.Series([])
        
        
        # Populate the crm_batch column for CRMs only (Cast = batch for CRMs)
        self.df[column_names.crm_batch].loc[self.df[column_names.station] == column_names.crmid] = self.df[column_names.cast]
        print("Using Cast number for CRMs as crm_batch input.")
        
        # Set DIC flag and TA flag to 0 if -999
        self.df[column_names.dicflag].loc[self.df[column_names.dicflag] == -999] = 0
        self.df[column_names.taflag].loc[self.df[column_names.taflag] == -999] = 0
        
        # Convert Date and Time format, create datetime column 
        self.df[column_names.analysisdate] = pd.to_datetime(self.df[column_names.analysisdate], format = "%Y%m%d" )
        #self.df[column_names.analysistime] = pd.to_datetime(self.df[column_names.analysistime], format = "%H%M%S" )
        self.df[column_names.analysis_datetime] = pd.to_datetime(self.df[column_names.analysisdate].astype(str) + " " + self.df[column_names.analysistime].astype(str).str.zfill(6), 
               format = "%Y-%m-%d %H%M%S")
        
        
        # <ED - 20210518 - Norwich - added
        # Check if these columns exist and insert if not
        check_columns = [column_names.sampling_date]
        for c in check_columns: 
            if c in self.df.columns:
                print(c + " found in the columns of the DataFrame. Make sure it's correctly populated, i.e. not empty.")
            else: 
                print(c + " is not found in the columns of the DataFrame. Creating it here, but will remain empty.")
                print("Needs to be populated before proceeding with further calculations within this class!" )
                self.df[c] = pd.Series([])
        # ED - 20210518 - Norwich > 
    
        
        # Print warnings/reminders
        print("Set an alternative sample temperature in sample_temperature_alternative column, if relevant!")
        
        
        #######################################################################################################################
    def add_dickson(self, dickson_location = "PICCOLO_pkg/crm_batch_lookup.csv"):
        """
        Looks up the certified concentrations for the relevant CRM batch. 
        CRM batch is specified in DataFrame of self as crm_batch (pd.Series)
        
        args: dickson_location  = string with path to csv file containing the CRM certified values in a particular format. Path is relative to location of current py file. 
        returns: no explicit returns, but does insert the certified concentrations into the rows of the DataFrame of self for CRM runs. 
        
        """
        print("Location of CRM Dickson info: " + dickson_location + ". Change if necessary.")
        
        # Find all unique batch numbers
        for b in self.df[column_names.crm_batch][self.df[column_names.station] == column_names.crmid].unique():
            # Insert the Dickson Certified Values for the CRMs
            self.crm_dickson = crmbatch.findcrmbatch(b, dickson_location) # < 20210407 - ED - UEA - b[0] to b >
                
            # For all CRM runs with this batch number, replace salinity, phosphate, silicate, and nitrate values with those of the CRM batch
            cond = (self.df[column_names.station] == column_names.crmid) & (self.df[column_names.crm_batch] == b)# < 20210407 - ED - UEA - b[0] to b >
            self.df[column_names.salinity][cond] = self.crm_dickson["salinity"]
            self.df[column_names.phosphate][cond] = self.crm_dickson["phos"]
            self.df[column_names.silicate][cond] = self.crm_dickson["sil"]
            self.df[column_names.nitrate][cond] = self.crm_dickson["nitrate"]
            
            # < ED - 20210512 - Norwich - added
            # Also populate the certified DIC and TA for the CRMs in the columns already created in the initialisation of the Class
            self.df[column_names.dic_certified_umolkg][cond] = self.crm_dickson["dic"]
            self.df[column_names.ta_certified_umolkg][cond] = self.crm_dickson["alk"]
            # ED - 20210512 - Norwich > 
        
        print("Added the CRM information for salinity, phosphate, silicate, and nitrate.")

        #######################################################################################################################
    def calculate_dic(self, use_alternative_temp = False):
        """
        """
        # Calculate the DIC concentration for all samples using the Dickson CRMs
        
        # Set the temperature you want to use as sample temperature.
        sampletemp = self.df[column_names.sampletemp]
        print("Using sample temperature as temperature. You can change temperature to use in PICCOLO_calc_DIC.py")
        if use_alternative_temp == True: 
            sampletemp = self.df[column_names.sample_temperature_alternative]
            print("Using sample_temperature_alternative as temperature. ")
        
        
        
        # Calculate the net DIC counts by removing the background DIC counts
        # Counts - (blanks [cts/min]/60)*seconds [sec])
        self.df[column_names.netdiccounts] = self.df[column_names.diccounts] - ((self.df[column_names.dicbackground]/60)*self.df[column_names.dictime])
        
        #######################################################################################################################
        # Volume correction for addition of mercuric chloride 
        # 100uL per 0.5L or 50uL per 0.25L of saturated HgCl2 solution
        # 0.5001/0.500 = 0.25005/0.250
        self.df[column_names.diccounts_volcorr] = self.df[column_names.netdiccounts]*(0.5001/0.5000) # < 20210311 - ED - Norwich - changed diccounts --> netdiccounts > 
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
        
        self.df[column_names.dic_sw_density] = RohSMOW + A*self.df[column_names.salinity] + B*self.df[column_names.salinity]**1.5 + C*self.df[column_names.salinity]**2
                
        #######################################################################################################################
        # Calculate mass of pipette volume [g/pipette]
        # [g / pip] = SWDensity [g / 1000 ml] / 1000 * Pipette Volume [ml]
        self.df[column_names.dicpipmass] = (self.df[column_names.dic_sw_density]/1000)*self.df[column_names.dicpipettevol]
        
        # To correct for glass expansion at measurement temperature, correct soelf.dicpipmass with aa3 
        # aa3 = Glass expansion factor relative to calibration at 25C
        # aa3 = 1 + 1e-5 * (self.sampletemp - 25)
        # self.dicpipmass = (self.dicswdensity/1000)*self.dicpipettevol*aa3
        # <<<<<< !!!!!!!!!!!!! >>>>>> CHECK THIS: why has this not been applied? 
        
        #######################################################################################################################        
        # Calculate umol DIC in the pipette volume titrated [umol/pipette] for each CRM, using CRM batch values to get umol per pipette
        # [umol / pip] = CRM DIC [umol/kg] / 1000 * pip mass [g/pip]
        self.df[column_names.dic_crm_pip_umol] = (self.crm_dickson["dic"]/1000) * self.df[column_names.dicpipmass][self.df[column_names.station] == column_names.crmid] # maintains same index!
        
        #######################################################################################################################
        # Calculate counts per mol in pipette volume titrated [CTS/umol] for each CRM, using CRM values to get counts per umol
        #[CTS / umol] = pip CTS [CTS / pip] / pip umol [umol / pip]
        # < 20210311 - ED - Norwich - changed diccounts --> netdiccounts 
        self.df[column_names.dic_crm_counts_per_umol] = self.df[column_names.netdiccounts][self.df[column_names.station] == column_names.crmid]/self.df[column_names.dic_crm_pip_umol][self.df[column_names.station] == column_names.crmid] 
        # 20210311 - ED - Norwich > 
        
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
        self.df[column_names.dic_crm_counts_per_umol_mean_allcells] = np.nanmean(self.df[column_names.dic_crm_counts_per_umol])
        self.df[column_names.dic_crm_counts_per_umol_std_allcells] = np.nanstd(self.df[column_names.dic_crm_counts_per_umol])
        
        # WARNING: these function are also applied to the CRMs. The latter are thus divided by their own collective mean. 
        # calculating calibrated moles titrated per pipette [umol]
        # [umol] = counts per pipette [CTS/pip] / counts per umol [CTS/umol]
        
        # DIC concentration based on average of ALL CRM runs in summary file
        self.df[column_names.dic_sample_pip_umol_allcells] = self.df[column_names.netdiccounts]/self.df[column_names.dic_crm_counts_per_umol_mean_allcells]
        
        # Calculate sample concentration [umol/l]
        self.df[column_names.dic_sample_umol_per_l_allcells] = (self.df[column_names.dic_sample_pip_umol_allcells]/self.df[column_names.dicpipettevol])*1000 # dicpipettevol is in mL

        # Calculate sample concentration [umol/kg]
        self.df[column_names.dic_sample_umol_per_kg_allcells] = (self.df[column_names.dic_sample_umol_per_l_allcells]/self.df[column_names.dic_sw_density])*1000 # dic_sw_density is in [kg/m3] or [g/1000ml]
        
        ###############################################################
        
        # Create a list of all unique analysis days (i.e. cells)
        index_days = self.df[column_names.dic_crm_counts_per_umol].groupby([self.df[column_names.diccellnum]]).size().index #Int64Index
        print("There are "+str(len(index_days))+" unique analysis days (i.e. DIC cells) during which samples were run.")

        # For each of the unique analysis days (i.e. cells), determine the mean of the CRMs' measured DIC counts per umol
        crm_mean_per_cell = self.df[column_names.dic_crm_counts_per_umol].groupby([self.df[column_names.diccellnum]]).transform("mean").unique()
        crm_std_per_cell = self.df[column_names.dic_crm_counts_per_umol].groupby([self.df[column_names.diccellnum]]).transform("std").unique()

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
        for i in self.df[column_names.diccellnum]:
            if i in crm_mean_dict:
                crm_mean_serieslist.append(crm_mean_dict[i])
                crm_std_serieslist.append(crm_std_dict[i])
        self.df[column_names.dic_crm_counts_per_umol_mean_percell] = crm_mean_serieslist # Add list to data frame
        self.df[column_names.dic_crm_counts_per_umol_std_percell] = crm_std_serieslist
        
        # Do exactly the same as above for "all", but now using the CRM average per analysis day (cell)
        self.df[column_names.dic_sample_pip_umol] = self.df[column_names.netdiccounts]/self.df[column_names.dic_crm_counts_per_umol_mean_percell]
        self.df[column_names.dic_sample_umol_per_l] = (self.df[column_names.dic_sample_pip_umol]/self.df[column_names.dicpipettevol])*1000 # dicpipettevol is in mL
        self.df[column_names.dic_sample_umol_per_kg] = (self.df[column_names.dic_sample_umol_per_l]/self.df[column_names.dic_sw_density])*1000 # dic_sw_density is in [kg/m3] or [g/1000ml]
                
        self.mean_diff = self.df[column_names.dic_sample_umol_per_kg] - self.df[column_names.dic_sample_umol_per_kg_allcells]
    
    
    #######################################################################################################################
    def calculate_dic_driftcorr(self, tobe_driftcorr_cell_list, time_frame): 
        """
        Apply drift corrections to given cells or range of cells. 

        INPUT: 
        tobe_driftcorr_cell_list: list containing numbers of DIC cells
        time_frame: "single_cells" or "multiple_cells": whether to apply the drift correction on each individual cell (single cells), or over all cells between (and including) two given cells, i.e. a drift over a longer period of time. 

        OUTPUT: 
        New columns with drift corrected values: 
            - dic_driftcorr_flag
            - dic_crm_counts_per_umol_driftcorr
            - dic_sample_pip_umol_driftcorr
            - dic_sample_umol_per_l_driftcorr
            - dic_sample_umol_per_kg_driftcorr

        Note: this method will override any existing values in the following columns: dic_driftcorr_flag, dic_crm_counts_per_umol_driftcorr,  dic_sample_pip_umol_driftcorr, dic_sample_umol_per_l_driftcorr, dic_sample_umol_per_kg_driftcorr

        Note: This method used to be according to analysis dates, i.e. drift corrections were calculated for given analysis dates instead of DIC cell numbers (changed on 20240723). 
        """
       
    
        # Check if the flag for DIC drift correction exists
        if "dic_driftcorr_flag" in self.df.columns:
            # If it exists, it means that this method has been applied before. Previous results will be overwritten
            print("Drift corrections have likely been done on this dataset before! Previous drift-corrected values will be overwritten. ")
        else: 
            print("First drift correction")
        
        # Create a column (or overwrite the existing one) that contains the flag for whether samples have been drift corrected
        # 0 = not drift corrected, 1 = drift corrected. Default is set to 0. 
        self.df[column_names.dic_driftcorr_flag] = 0 
        # Create a new column (or override column) that will contain the DIC counts per umol that will have been drift corrected
        self.df[column_names.dic_crm_counts_per_umol_driftcorr] = np.float("NaN")

        # Check if list (given in argument) is empty
        if not tobe_driftcorr_cell_list: 
            print("List is empty, no drift corrections done.")

        # If not empty, do the drift correction
        else: 
            
            fig = plt.figure(figsize = (15, 10)) # Create a figure. Subplots will be added along the way
            x, y = 2, 3 # subplot dimensions
            z = 0
        
            try: 
                if time_frame == "single_cells": # If drift correction needs to be done per cell given in tobe_driftcorr_cell_list
                    for corr_cell in tobe_driftcorr_cell_list:
                        z = z + 1 # Add a subplot to figure

                        # Change the flag for DIC drift corrected from 0 to 1 for the given cell number
                        self.df[column_names.dic_driftcorr_flag].where(self.df[column_names.diccellnum] != corr_cell, 1)

                        # Extract the relevant variables and convert them to an array and reshape so they can be used by the LinearRegression class
                        # Note that the sklearn package doesn't like datetime type values, so datetime is converted to a timestamp. 
                        # Also make sure that the arrays used to fit the linear regression do not contains NaNs (sklearn doesn't like it)
                        cond_crm_time_array = (self.df[column_names.diccellnum] == corr_cell) & (self.df[column_names.station] == column_names.crmid) & (~self.df[column_names.dic_crm_counts_per_umol].isna())
                        crm_time_array = self.df[column_names.analysis_datetime].loc[cond_crm_time_array].apply(lambda x: x.timestamp()).values.reshape(-1,1)

                        samples_time_array = self.df[column_names.analysis_datetime].loc[self.df[column_names.diccellnum] == corr_cell].apply(lambda x: x.timestamp()).values.reshape(-1,1)

                        cond_crm_ctsperumol_array = (self.df[column_names.diccellnum] == corr_cell) & (self.df[column_names.station] == column_names.crmid) & (~self.df[column_names.dic_crm_counts_per_umol].isna())
                        crm_ctsperumol_array = self.df[column_names.dic_crm_counts_per_umol].loc[cond_crm_ctsperumol_array].values.reshape(-1,1)
                        

                        # Use linear regression to interpolate the CRM counts/umol at the time of each sample
                        linear_regressor = LinearRegression()
                        linear_regressor.fit(crm_time_array, crm_ctsperumol_array)
                        Y_pred = linear_regressor.predict(crm_time_array)

                        # Create a list containing the predicted values; this is only used for plotting afterwards. 
                        crm_drift_corr_plt = []

                        # Per sample, determine the counts/umol based on the linear regression
                        # This is done per sample's index rather than for all samples at the same time to make sure the result ends up in the right place in the dataframe
                        for i in samples_time_array:
                            crm_drift_corr = linear_regressor.predict(i.reshape(-1,1))[0][0]

                            # Save the drift corrected value
                            # Note that i (i.e. the timestamp) needs to be converted back to datetime
                            i_datetime = pd.to_datetime(i[0], unit='s')
                            ix = self.df.loc[(self.df[column_names.diccellnum] == corr_cell) & (self.df[column_names.analysis_datetime] == i_datetime)].index
                            self.df.loc[ix[0], column_names.dic_crm_counts_per_umol_driftcorr] = crm_drift_corr

                            # Also save the predicted value to the list (used for plotting only)
                            crm_drift_corr_plt.append(crm_drift_corr)


                        # Create a simple plot just to show what has been done 
                        ax = fig.add_subplot(x, y, z)

                        ax.scatter(crm_time_array, crm_ctsperumol_array)
                        ax.plot(crm_time_array, Y_pred, color="red")
                        ax.scatter(samples_time_array, crm_drift_corr_plt)
                        ax.set_xlabel("Time")
                        ax.set_ylabel("counts per umol")
                        ax.set_title(str(corr_cell))

                elif time_frame == "multiple_cells": # If cell corrections needs to be done from the first cell until the last cell given in tobe_driftcorr_cell_list, including all cells in between
                
                    z = 1
                
                    # tobe_driftcorr_list needs to have a minimum of 2 elements. The first and the last element indicate the cells between which linear interpolation/drift correction will be applied
                    cell_min = tobe_driftcorr_cell_list[0]
                    cell_max = tobe_driftcorr_cell_list[-1]
                
                    # Flag the days that are drift corrected
                    self.df[column_names.dic_driftcorr_flag][(self.df[column_names.diccellnum] >= cell_min) & (self.df[column_names.diccellnum] <= cell_max)] = 1 
                
                    # Extract the relevant variables and convert to array + reshape for use in LinearRegression
                    # Note that the sklearn package doesn't like datetime type values, so datetime is converted to a timestamp. 
                    # Also make sure that the arrays used to fit the linear regression do not contains NaNs (sklearn doesn't like it)
                    cond_crm_time_array = (self.df[column_names.diccellnum] >= cell_min) & (self.df[column_names.diccellnum] <= cell_max) & (self.df[column_names.station] == column_names.crmid) & (~self.df[column_names.dic_crm_counts_per_umol].isna())
                    crm_time_array = self.df[column_names.analysis_datetime].loc[cond_crm_time_array].apply(lambda x: x.timestamp()).values.reshape(-1,1)
                    
                    samples_time_array = self.df[column_names.analysis_datetime].loc[(self.df[column_names.diccellnum] >= cell_min) & (self.df[column_names.diccellnum] <= cell_max) & (self.df[column_names.station] != column_names.crmid)].apply(lambda x: x.timestamp()).values.reshape(-1,1)

                    cond_crm_ctsperumol_array = (self.df[column_names.diccellnum] >= cell_min) & (self.df[column_names.diccellnum] <= cell_max) & (self.df[column_names.station] == column_names.crmid) & (~self.df[column_names.dic_crm_counts_per_umol].isna())
                    crm_ctsperumol_array = self.df[column_names.dic_crm_counts_per_umol].loc[cond_crm_ctsperumol_array].values.reshape(-1,1)
                
                    # Use linear regression to interpolate CRM counts/umol at time of each sample
                    linear_regressor = LinearRegression()
                    linear_regressor.fit(crm_time_array, crm_ctsperumol_array)
                    Y_pred = linear_regressor.predict(crm_time_array)
                
                    # Create a list containing the predicted values; this is used for plotting afterwards. 
                    crm_drift_corr_plt = []
                
                    # Per sample, determine the counts/umol based on the linear regression
                    for i in samples_time_array: 
                        crm_drift_corr = linear_regressor.predict(i.reshape(-1,1))[0][0]

                        # Save the drift corrected value
                        # Note that i (i.e. the timestamp) needs to be converted back to datetime
                        i_datetime = pd.to_datetime(i[0], unit='s')
                        ix = self.df.loc[(self.df[column_names.diccellnum] >= cell_min) & (self.df[column_names.diccellnum] <= cell_max) & (self.df[column_names.analysis_datetime] == i_datetime)].index
                        self.df.loc[ix[0], column_names.dic_crm_counts_per_umol_driftcorr] = crm_drift_corr
                                            
                        # Save the predicted value to the list (used for plotting only)
                        crm_drift_corr_plt.append(crm_drift_corr)
                    
                    # Create a simple plot just to show what has been done 
                    ax = fig.add_subplot(x, y, z)

                    ax.scatter(crm_time_array, crm_ctsperumol_array)
                    ax.plot(crm_time_array, Y_pred, color="red")
                    ax.scatter(samples_time_array, crm_drift_corr_plt)
                    ax.set_xlabel("Time")
                    ax.set_ylabel("counts per umol")
                    ax.set_title(str(cell_min) + " to " + str(cell_max))
            
                else:
                    raise IncorrectArgumentError 
        
                # After all the necessary drift corrections have been completed, re-calculate the drift-corrected DIC concentrations 
                # We do not need to check if the new columns already exist. If they don't, they'll automatically be created. If they do, then values are recalculated. 
                # But previously calculated drift corrected values are not overwritten above, so they should still end up with the same calculated value (umol/kg)
                self.df[column_names.dic_sample_pip_umol_driftcorr] = self.df[column_names.netdiccounts]/self.df[column_names.dic_crm_counts_per_umol_driftcorr]
                self.df[column_names.dic_sample_umol_per_l_driftcorr] = (self.df[column_names.dic_sample_pip_umol_driftcorr]/self.df[column_names.dicpipettevol])*1000
                self.df[column_names.dic_sample_umol_per_kg_driftcorr] = (self.df[column_names.dic_sample_umol_per_l_driftcorr]/self.df[column_names.dic_sw_density])*1000

                print("DIC content has been re-calculated according to drift corrected CRM values.")
        
            except IncorrectArgumentError: 
                print("IncorrectArguentError: time_frame argument can only be single_days or multiple_days. Correct and try again.")
        
        
        
        
        #######################################################################################################################    
    def extract_crms(self):
        # Extract the CRM data into a separate data frame
        self.crm = self.df[self.df[column_names.station] == column_names.crmid] # this is a dataframe for just the CRM measurements
        # Index has not been re-set. 
        
        
    #######################################################################################################################
    # Quality check the CRMs
    # This method will plot either the counts or the calculated concentrations for the CRMs per analysis day
    def qc_crms(self, use="counts"):
        """
        Plot CRM counts or concentration per DIC cell. 
        """
        
        # WARNING: this function works when analyses have been done within 24 hours, i.e. on the same date (or at least labelled like that)
        # Better would be to distinguish analysis days on cell number, instead of date. Maybe change this at a later point! 
        
        # 'use' argument is used to specify whether to plot the CRM "counts" or the "concentration". Default is counts. 
        
        # Determine the number of analysis days according to cell number
        analysis_cells = self.crm[column_names.diccellnum].unique()
        x = 6 # Number of rows of subplots # < 20240310 - ESD SDA - changed from 4 to 6 >
        y = 4 # Number of columns of subplots

        # Create empty list for analysis days that will need drift checking later
        check_drift = []

        # Start figure
        # Plot the CRMs measured for each subplot vs analysis time
        fig = plt.figure(figsize = (20, 20))
        

        # For each cell number, generate a subplot
        for i in np.arange(len(analysis_cells)): 
            if use == "counts":
                df_mini = self.crm.loc[self.crm[column_names.diccellnum] == analysis_cells[i]].dropna(subset = [column_names.netdiccounts]) # Create a mini data frame with just the CRM data per analysis day ...
                # AND remove any rows where the netdiccounts for CRMs is NaN!
                crm_dic = df_mini[column_names.netdiccounts]
                ylim_min_text = np.min(self.crm[column_names.netdiccounts])
                ylim_max_text = np.max(self.crm[column_names.netdiccounts])
                ylim_min = ylim_min_text-100
                ylim_max = ylim_max_text+100
            elif use == "concentration":
                df_mini = self.crm.loc[self.crm[column_names.diccellnum] == analysis_cells[i]].dropna(subset = [column_names.dic_sample_umol_per_kg]) # Create a mini data frame with just the CRM data per analysis day ...
                # AND remove any rows where the netdiccounts for CRMs is NaN!
                crm_dic = df_mini[column_names.dic_sample_umol_per_kg]
                ylim_min_text = np.min(self.crm[column_names.dic_sample_umol_per_kg])
                ylim_max_text = np.max(self.crm[column_names.dic_sample_umol_per_kg])
                ylim_min = ylim_min_text-1
                ylim_max = ylim_max_text+1
            else:
                print("Set argument of qc_crms() function to either 'counts' or to 'concentration'")
            
            # < ED - 20210528 - Norwich - changed crm_time array to hour fraction to accomodate past midnight runs
            # Create a series with the hours 
            # crm_time = df_mini[column_names.analysistime] 
            if len(df_mini[column_names.analysisdate].loc[df_mini[column_names.diccellnum] == analysis_cells[i]].unique()) == 1: 
                # If all runs of the dic cell were done on the same day, create crm_time series = hours since start of the day
                crm_time = df_mini[column_names.analysis_datetime].dt.hour +(df_mini[column_names.analysis_datetime].dt.minute/60) + (df_mini[column_names.analysis_datetime].dt.minute/3600)
            elif len(df_mini[column_names.analysisdate].loc[df_mini[column_names.diccellnum] == analysis_cells[i]].unique()) > 1: 
                # First calculate the number of hours since the start of each day
                crm_time = df_mini[column_names.analysis_datetime].dt.hour +(df_mini[column_names.analysis_datetime].dt.minute/60) + (df_mini[column_names.analysis_datetime].dt.minute/3600)
                # If some runs of the same dic cell were done on two different dates (e.g. past mid-night), add 24 hours to the total number of hours for the subsequent day
                add_hours = 24
                for d in np.sort(df_mini[column_names.analysisdate].loc[df_mini[column_names.diccellnum] == analysis_cells[i]].unique())[1:]:
                    # Run through each day after the first day (should be max. 1) and add 24 hours to it
                    ix = df_mini[df_mini[column_names.analysisdate] == d].index # find indeces to which this applies 
                    crm_time.loc[ix] = df_mini.loc[ix][column_names.analysis_datetime].dt.hour +(df_mini.loc[ix][column_names.analysis_datetime].dt.minute/60) + (df_mini.loc[ix][column_names.analysis_datetime].dt.minute/3600) + add_hours
                    add_hours = add_hours + 24 # only relevant if dic cell covers more than 2 days... should never happen... 
            else:
                print("No analysis dates for this cell.")
            # ED - 20210528 - Norwich >
            
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
            title = "Cell No. "+ str(analysis_cells[i])
            ax.set_title(title)
            
            # Calculate the standard deviation of all CRM netdiccounts values per analysis day
            std_crms_day = np.std(crm_dic_array)
            text2 = "1$\sigma$="+str(std_crms_day)

            # Calculate the min and max of the regression line and add it to the graph
            Y_pred_range = np.max(Y_pred) - np.min(Y_pred)
            if use == "counts":
                text1 = "diff="+str(Y_pred_range)+" cts"
            if use == "concentration":
                text1 = "diff="+str(Y_pred_range)+" $\mu$mol kg$^{-1}$"
            
            # Add the limitations for the y axis
            ax.set_ylim(ylim_min, ylim_max)
            
            if use == "counts":
                if Y_pred_range > 100:
                    textcolour = "r"
                    check_drift.append(analysis_cells[i])
                    printtext = "The following cell numbers have drift > 100 counts during the day: "
                else:
                    textcolour = "k"
            elif use == "concentration":
                if Y_pred_range > 3:
                    textcolour = "r"
                    check_drift.append(analysis_cells[i])
                    printtext = "The following cell numbers have drift > 3 $\mu$mol kg$^{-1}$ during the day: "
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
        
        if use == "counts":
            self.qc_crms_cts_fig = fig        
            print("Figure saved as self.qc_crms_cts_fig")
        elif use == "concentration":
            self.qc_crms_concentration_fig = fig        
            print("Figure saved as self.qc_crms_concentration_fig")
        

    
    
    ####################################################################################################################### 
    def dic_x_chart(self):
        """
        X Chart according to SOP22 in Dickson et al. 
        
        :returns: creates an instance for plotly figure
        
        Only if interactive mode == "on"
        """
        
        if column_names.interactive == "on":
            
            
            # Set lay out properties
#            layout = go.Layout(hovermode="closest", # to make sure hover text appears for points that share an x value
#                               title = "$\overline{X}$ Chart DIC Concentrations for CRMs", 
#                               showlegend=False, 
#                              )
            
            sequence = np.arange(len(self.crm.index))
            crm_conc = self.crm[column_names.dic_sample_umol_per_kg]
            
            x = np.mean(self.crm[column_names.dic_sample_umol_per_kg])
            s = np.std(self.crm[column_names.dic_sample_umol_per_kg])
            ucl = x + 3*s
            uwl = x + 2*s
            lwl = x - 2*s
            lcl = x - 3*s
            
            # Create figure
            fig = go.Figure()#layout=layout)
            
            fig.add_hline(y=x, line_width=3, line_color="DarkSlateGrey", 
                          annotation_text="mean", annotation_position="bottom right", annotation_font_size=14, annotation_font_color="black")
            fig.add_hline(y=uwl, line_width=3, line_dash="dash", line_color="darkturquoise", 
                          annotation_text="+2sigma", annotation_position="bottom right", annotation_font_size=14, annotation_font_color="black")
            fig.add_hline(y=lwl, line_width=3, line_dash="dash", line_color="darkturquoise", 
                          annotation_text="-2sigma", annotation_position="bottom right", annotation_font_size=14, annotation_font_color="black")
            fig.add_hline(y=ucl, line_width=3, line_dash="dot", line_color="LightSeaGreen", 
                          annotation_text="+3sigma", annotation_position="bottom right", annotation_font_size=14, annotation_font_color="black")
            fig.add_hline(y=lcl, line_width=3, line_dash="dot", line_color="LightSeaGreen", 
                          annotation_text="-3sigma", annotation_position="bottom right", annotation_font_size=14, annotation_font_color="black")
            
#            fig.layout.update(
#                    shapes = [
#                            # x
#                            go.layout.Shape(
#                            type="line",
#                                x0=np.min(sequence),
#                                y0=x,
#                                x1=np.max(sequence),
#                                y1=x,
#                                line=dict(
#                                    color="RoyalBlue",
#                                    width=3,                                
#                                    ),
#                                ),
#                            
#                            # ucl
#                            go.layout.Shape(
#                                type="line",
#                                x0=np.min(sequence),
#                                y0=ucl,
#                                x1=np.max(sequence),
#                                y1=ucl,
#                                line=dict(
#                                    color="LightSeaGreen",
#                                    width=3,
#                                    dash="dot",
#                                    ),
#                                ),
#                            # uwl
#                            go.layout.Shape(
#                                    type="line",
#                                    x0=np.min(sequence),
#                                    y0=uwl,
#                                    x1=np.max(sequence),
#                                    y1=uwl,
#                                    line=dict(
#                                        color="MediumPurple",
#                                        width=3,
#                                        dash="dashdot",
#                                        ),
#                                    ),
#                            #lwl
#                            go.layout.Shape(
#                                    type="line",
#                                    x0=np.min(sequence),
#                                    y0=lwl,
#                                    x1=np.max(sequence),
#                                    y1=lwl,
#                                    line=dict(
#                                        color="MediumPurple",
#                                        width=3,
#                                        dash="dashdot",
#                                        ),
#                                    ),
#                                            
#                            # lcl
#                            go.layout.Shape(
#                                    type="line",
#                                    x0=np.min(sequence),
#                                    y0=lcl,
#                                    x1=np.max(sequence),
#                                    y1=lcl,
#                                    line=dict(
#                                        color="LightSeaGreen",
#                                        width=3,
#                                        dash="dot",
#                                        ),
#                                    ),          
#                        ]
#                    )
            
            
            ###
            fig.add_trace(go.Scatter(
                x = sequence, 
                y = crm_conc, 
                mode = "markers", 
                hovertext = np.array(self.crm[column_names.sample])
            ))
            
#            fig.add_trace(go.Scatter(
#                x=[1.5, 1.5, 1.5, 1.5, 1.5],
#                y=[ucl+1, uwl+1, x+1, lwl+1, lcl+1],
#                text=["UCL", "UWL", "$\overline{X}$", "LWL", "LCL"],
#                mode="text",
#            ))
            
    
            # Set axes titles
            fig.update_xaxes(title_text = "Sequence")
            fig.update_yaxes(title_text = "CRM DIC [umol/kg]")
            
            fig.update_traces(showlegend=False)
            fig.update_layout(hovermode="closest")
            fig.update_layout(font_size=14)
    
            fig.show()
            
            self.dic_x_chart_fig = fig
            print("Figure saved as self.dic_x_chart_fig")
            
        else: 
            print("Interactive plotting is not switched on. Switch on to see plots generated with Plotly. Otherwise, manually create plots. ")
    
    #######################################################################################################################
    def dic_r_chart(self):
        """
        R chart according to SOP22 in Dickson et al.    
        
        Only if interactive mode == "on"
        """
        
        if column_names.interactive == "on":
        
            # Set lay out properties
            layout = go.Layout(hovermode="closest", # to make sure hover text appears for points that share an x value
                               title = "$\overline{R}$ Chart DIC Concentrations for CRMs", 
                               showlegend=False, 
                              )
            
            
            
            crm_rdiff = abs(self.crm[column_names.dic_sample_umol_per_kg].groupby(self.crm[column_names.niskin]).diff())
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
                hovertext = np.array(self.crm[column_names.sample])
            ))
            
            fig.add_trace(go.Scatter(
                x=[1.5, 1.5, 1.5],
                y=[ucl+0.5, uwl+0.5, r+0.5],
                text=["UCL", "UWL", "$\overline{R}$]"],
                mode="text"
            ))
            
    
            # Set axes titles
            fig.update_xaxes(title_text = "Sequence")
            fig.update_yaxes(title_text = "CRM R [umol/kg]")
    
            fig.show()
            
            self.dic_r_chart_fig = fig
            print("Figure saved as self.dic_r_chart_fig")
            
        else: 
            print("Interactive plotting is not switched on. Switch on to see plots generated with Plotly. Otherwise, manually create plots. ")
    
    
    #######################################################################################################################
    def qc_conc_profiles(self, to_plot="dic_sample_umol_per_kg"): 
        """
        Plot the DIC concentrations calculated within this class for all samples. 
        Plotted in an interactive plot using plotly so it's easy for the user to identify samples. 
        1 plot (interactive!) per station/cast, colour markers according to cell id, pressure on y axis, concentration on x axis
        
        :args: use_dic = dic_sample_umol_per_kg or dic_sample_umol_per_kg_final
        
        :returns: creates an instance for plotly figure
        
        Only if interactive mode == "on"
        """
        
        if column_names.interactive == "on":

            num = 0
            r = 0 
            c = 0
            
            numberofcols = 6#5
            
            numstations = self.df[column_names.station].loc[self.df[column_names.station] != column_names.crmid].unique() # find number of stations (station value for CRMs is "nan", so we're excluding those)
            min_depth = np.nanmin(self.df[column_names.depth])
            max_depth = np.nanmax(self.df[column_names.depth])
            
            subplot_titles = []
            station_casts_dict = {}
            total_stationcasts = 0
            
            for i in numstations: 
                station_casts_dict[i] = self.df[column_names.cast].loc[self.df[column_names.station] == i].unique()
                for j in station_casts_dict[i]: 
                    s = "Station " + str(int(i)) + " Cast " + str(int(j))
                    subplot_titles.append(s)
                    total_stationcasts = total_stationcasts + 1
            
            fig = make_subplots(rows=int(len(total_stationcasts)/numberofcols)+1, cols=numberofcols, subplot_titles = subplot_titles, x_title = "DIC [umol kg-1]", y_title = "Depth [m]")
            
            # Create a column for colour scheme according to DIC cell 
            numcells = self.df[column_names.diccellnum].unique()
            cell_colours = ["firebrick", "darkorange", "blue",
                            "blueviolet", "brown", "burlywood", "cadetblue", "chartreuse", "chocolate", "coral", "darkolivegreen", "cornflowerblue",
                            "crimson", "cyan", "darkblue", "darkcyan","darkgoldenrod", "darkgray", "darkgreen",
                            "darkkhaki", "darkmagenta","darkorchid", "darkred", "darksalmon", "darkseagreen",
                            "darkslateblue", "darkslategray", "darkslategrey","darkturquoise", "darkviolet", "deeppink", "deepskyblue",
                            "dimgrey", "dodgerblue", "forestgreen", "cornsilk"]
            n = 0
            self.df["cell_colour"] = "black"
            thelegend = {}
            for i in numcells:
                self.df["cell_colour"][self.df[column_names.diccellnum] == i] = cell_colours[n]
                thelegend[cell_colours[n]] = i
                n = n + 1
            
            
            # Temporarily fill nans FOR SAMPLES ONLY in "dic_sample_umol_per_kg" (or "dic_sample_umol_per_kg_driftcorr") with zeros (plotly doesn't like the "nan"s)
            self.df[to_plot][self.df[column_names.station] != column_names.crmid] = self.df[to_plot].fillna(0)
            
            for station in numstations: 
                
                for cast in station_casts_dict[station]:
                
                    r = int(num/numberofcols)
                    c = num % numberofcols        
                    num = num + 1
                    
                    
            
                    fig.add_trace(
                        go.Scatter(
                            x = self.df[to_plot][(self.df[column_names.station] == station) & (self.df[column_names.cast] == cast)].values, 
                            y = self.df[column_names.depth][(self.df[column_names.station] == station) & (self.df[column_names.cast] == cast)].values, 
                            mode = "markers+text", 
                            text = self.df[column_names.diccellnum][(self.df[column_names.station] == station) & (self.df[column_names.cast] == cast)].tolist(), 
                            textposition = "bottom right", 
                            textfont = dict(size = 4), 
                            marker = dict(color = self.df["cell_colour"][(self.df[column_names.station] == station) & (self.df[column_names.cast] == cast)].values), 
                            hovertext = self.df[column_names.sample][(self.df[column_names.station]==station) & (self.df[column_names.cast] == cast)].values,
                        ),
                        row = r+1, 
                        col = c+1)
        
                    fig.update_yaxes(range = [max_depth, min_depth])
                
            
            fig.layout.update(height=3500, width=1000, showlegend = False, title_text = "Quality Check Samples", hovermode = "closest")
            fig.show()
            
            # Change any 0 back to nans in "dic_sample_umol_per_kg" or "dic_sample_umol_per_kg_driftcorr" column 
            self.df[to_plot][self.df[column_names.station] != column_names.crmid].replace(to_replace = 0, value = np.nan, inplace = True)
            
            # Return the figure
            # return fig
            self.qc_conc_profiles_fig = fig
            print("Figure saved as self.qc_conc_profiles_fig")
        else: 
            print("Interactive plotting is not switched on. Switch on to see plots generated with Plotly. Otherwise, manually create plots. ")

    #######################################################################################################################
    
    def qc_conc_timeseries(self):
        """
        
        Plot the DIC concentrations calculated within this class for all samples
        Plotted in an interactive plot using plotly so it's easy for the user to identify samples
        
        
        :returns: creates an instance for plotly figure
        
        Only if interactive mode == "on"
        
        """
        
        if column_names.interactive == "on":
        
            sdates= self.df[column_names.sampling_date][self.df[column_names.station] != column_names.crmid]
            single_measurements= self.df[column_names.dic_sample_umol_per_kg][self.df[column_names.station] != column_names.crmid]
            
            single_measurements_allcrm= self.df[column_names.dic_sample_umol_per_kg_allcells][self.df[column_names.station] != column_names.crmid]
            
            single_text = np.array(self.df[column_names.sample][self.df[column_names.station] != column_names.crmid])
            
            
            # Set lay out properties
            layout = go.Layout(hovermode="closest", title = "Check DIC Concentrations for Time Series", )
    
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
            
            self.qc_conc_timeseries_fig = fig
            print("Figure saved as self.qc_conc_timeseries_fig")
            
        else: 
            print("Interactive plotting is not switched on. Switch on to see plots generated with Plotly. Otherwise, manually create plots. ")

print("Imported CalcDIC class.")



