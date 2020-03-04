# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 13:16:08 2020

@author: dfeng
"""

import numpy as np
import pandas as pd


def savetxt_Traveltime_branch1(WS, Ttime, density, flow_index, concentrate, txtfile):
    """
    save travel time to a txt file
    output array: 
        branchID, segID, travel time, density, release_arm, solubility, flow_condition, concentration
    density=1 heavy
    density=0 light
    release_arm=1 East
    flow_index - high:3, medium:2, low:1
    """
        
    outarray = np.vstack((np.ones_like(Ttime), WS.segs1, Ttime, \
                          np.ones_like(Ttime)*density, np.ones_like(Ttime), \
                          np.ones_like(Ttime)*1, np.ones_like(Ttime)*flow_index, concentrate)).T
    np.savetxt(txtfile, outarray, fmt='%1.4e')
    


def savetxt_Traveltime_branch5(self, WS, Ttimes, density, flow_index, concentrates, txtfile):
    """
    output array: 
        branchID, segID, travel time, density, release_arm, solubility, flow_condition, concentration
    density=1 heavy
    density=0 light
    release_arm=5 West
    """
        
    Ttime1 = Ttimes[0]    ## travel times at branch 1
    Ttime5 = Ttimes[1]    ## travel times at branch 5
    
    concentrate1 = concentrates[0]
    concentrate5 = concentrates[1]
        
    outarray1 = np.vstack((np.ones_like(Ttime1), WS.segs1, Ttime1, \
                          np.ones_like(Ttime1)*density, np.ones_like(Ttime1)*5, \
                          np.ones_like(Ttime1), np.ones_like(Ttime1)*flow_index, concentrate1))
        
    #### important !! reverse WS.segs5, from 86 to 121
    outarray5 = np.vstack((np.ones_like(Ttime5)*5, WS.segs5[::-1], Ttime5, \
                          np.ones_like(Ttime5)*density, np.ones_like(Ttime5)*5, \
                          np.ones_like(Ttime5), np.ones_like(Ttime5)*flow_index, concentrate5))
        
    outarray = np.hstack((outarray5, outarray1)).T
        
    #np.savetxt(txtfile, outarray, fmt='%d')
    np.savetxt(txtfile, outarray, fmt='%1.4e')
    

    
def save_excel_Traveltime_branch1(WS, Ttime, density, flow_index, concentrate, water_level, excelfile):
    """
    save travel time info to excel .xlsx file: https://stackoverflow.com/questions/51904126/write-a-numpy-ndarray-to-an-xlsx-spreadsheet
    branchID, segID, travel time, density, release_arm, solubility, flow_condition, concentration, water level
    """
    outarray = np.vstack((np.ones_like(Ttime), WS.segs1, Ttime, \
                          np.ones_like(Ttime)*density, np.ones_like(Ttime), \
                          np.ones_like(Ttime)*1, np.ones_like(Ttime)*flow_index, concentrate, water_level)).T
        
    headers = ["branchID", "segID", "travel_time", "density", "release_arm", \
               "solubility", "flow_condition", "concentration", "water_level"]
                      
    df = pd.DataFrame(outarray, columns=headers)
        
    df.to_excel(excelfile, index=False)
    


def save_excel_Traveltime_branch5(WS, Ttimes, density, flow_index, concentrates, water_levels, excelfile):
    """
    save travel time info to excel .xlsx file for spill initiated at branch 5
    """
        
    Ttime1 = Ttimes[0]    ## travel times at branch 1
    Ttime5 = Ttimes[1]    ## travel times at branch 5
    
    concentrate1 = concentrates[0]
    concentrate5 = concentrates[1]
    
    water_level1 = water_levels[0]
    water_level5 = water_levels[1]
    
    outarray1 = np.vstack((np.ones_like(Ttime1), WS.segs1, Ttime1, \
                          np.ones_like(Ttime1)*density, np.ones_like(Ttime1)*5, \
                          np.ones_like(Ttime1), np.ones_like(Ttime1)*flow_index, \
                          concentrate1, water_level1))
        
    #### important !! reverse WS.segs5, from 86 to 121
    outarray5 = np.vstack((np.ones_like(Ttime5)*5, WS.segs5[::-1], Ttime5, \
                          np.ones_like(Ttime5)*density, np.ones_like(Ttime5)*5, \
                          np.ones_like(Ttime5), np.ones_like(Ttime5)*flow_index, \
                          concentrate5, water_level5))
        
    outarray = np.hstack((outarray5, outarray1)).T
        
    headers = ["branchID", "segID", "travel_time", "density", "release_arm", \
                   "solubility", "flow_condition", "concentration", "water_level"]
        
    df = pd.DataFrame(outarray, columns=headers)
        
    df.to_excel(excelfile, index=False)
    