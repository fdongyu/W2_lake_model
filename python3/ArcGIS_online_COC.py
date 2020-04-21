# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 12:59:18 2020

@author: dfeng
"""

import pandas as pd
from datetime import datetime
import numpy as np
import utm
import osgeo.ogr
from decimal import Decimal

from segmentation import W2_Segmentation

import pdb


class ArcGIS_online_map_COC(W2_Segmentation):
    """
    class to write data for each COC type into ArcGIS shapefile
    """
    
    def __init__(self, workdir, **kwargs):
        
        #self.workdir = workdir
        self.bathfile = '%s\\%s'%(workdir, 'Bth_WB1.npt')
        self.VisSeg2()
        
        self.readCOCtable()


    def readCOCtable(self):
        """
        read COC table, contaminant type, solubility and density
        """

        filename = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\treatability\COC_table_04132020.xlsx'
        
        xl = pd.ExcelFile(filename)
        df = xl.parse('Sheet1', skiprows=2)
        
        df = df.drop(df.index[0])
        
        contaminant = {}
        self.contaminant_names = []
        self.solubilities = []
        self.densities = []
        
        for i in range(1, len(df)):
            contaminant.update({df['CONTAMINANT'].iloc[i]: {'Solubility': df['Solubility Classification'].iloc[i], 'Density': df['Sink/Float'].iloc[i]} })
            
            if df['Sink/Float'].iloc[i] != 'Sink':
                self.contaminant_names.append(df['CONTAMINANT'].iloc[i])
                self.solubilities.append(df['Solubility Classification'].iloc[i])
                self.densities.append(df['Sink/Float'].iloc[i])
        
        self.solubilities = np.asarray(self.solubilities)
        self.solubilities[self.solubilities=='Slightly soluble'] = 'Insoluble'
        
        self.solubilities[self.solubilities=='nan'] = 'Soluble'
        
        #pdb.set_trace()
    
    
    def write_shapefile_combined(self, shpname):  
        #### write data for each category 
        #### soluble, insoluble_float, insoluble_sink 
        
    
        #### Create the shapefile
        # Create the projection
        spatialReference = osgeo.osr.SpatialReference()
        spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
        
        # Create the shape file
        outfile = r'ArcGIS_online\%s'%shpname
        driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
        shapeData = driver.CreateDataSource(outfile)
        
        # Create the layer
        layer = shapeData.CreateLayer('Contour', spatialReference, osgeo.ogr.wkbPoint)
        layerDefinition = layer.GetLayerDefn()
        
        # Create fields containing segment infos
        field_def = osgeo.ogr.FieldDefn('COC', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('BranchID', osgeo.ogr.OFTInteger)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('SegID', osgeo.ogr.OFTInteger)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('Lon', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('Lat', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('T (day)', osgeo.ogr.OFTInteger)
        layer.CreateField(field_def)
        
        ## density - type: string, option: light-0, heavey-1 
        field_def = osgeo.ogr.FieldDefn('Density', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        ## initial release location - type: string, option: East-1, West-5
        field_def = osgeo.ogr.FieldDefn('Initial', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        ## solubility
        field_def = osgeo.ogr.FieldDefn('Solubility', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        ## flow condition
        field_def = osgeo.ogr.FieldDefn('Flow', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        ## concentration
        field_def = osgeo.ogr.FieldDefn('C (mg/L)', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        ## water surface elevation
        field_def = osgeo.ogr.FieldDefn('WSE (ft)', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
        
        ## distance to WTP gate
        field_def = osgeo.ogr.FieldDefn('D (ft)', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
        
        
        def add_feature(layer, COCname, branchID, segs, lines, lon, lat, Ttime, density, Initial_loc, solubility, flows, concentration, water_level, dist):
            """
            function that adds feature to layer
            """    
            ctr=0
            for i in range(len(lines)):
                ctr+=1
                point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
                # Add points individually to the line
                #xy = lines[i]
        
                #line.AddPoint_2D(xy[0][0],xy[0][1])
                #line.AddPoint_2D(xy[1][0],xy[1][1])
                point.AddPoint(lon[i], lat[i])
                # Update the feature with the line data
                featureIndex = ctr
                feature = osgeo.ogr.Feature(layerDefinition)
                #feature.SetStyleString("PEN(c:r,w:5px)")   
                feature.SetGeometry(point)
                feature.SetFID(featureIndex)
                feature.SetGeometryDirectly(point)
                
                # Set the attribute table
                feature.SetField('COC', COCname)
                feature.SetField('BranchID', int(branchID[i])) 
                feature.SetField('SegID', int(segs[i]))   # convert to int() is necessary, osgeo cannot recognize numpy int32 type
                feature.SetField('Lon', "{:.3f}".format(lon[i]))
                feature.SetField('Lat', "{:.3f}".format(lat[i]))
                #feature.SetField('Lon_east', "{:.3f}".format(eastlon[i]))
                #feature.SetField('Lat_east', "{:.3f}".format(eastlat[i]))
                feature.SetField('T (day)', int(Ttime[i]))
                feature.SetField('Density', density[i])
                feature.SetField('Initial', Initial_loc[i])
                feature.SetField('Solubility', solubility[i])
                feature.SetField('Flow', flows[i])
                feature.SetField('C (mg/L)', concentration[i])
                feature.SetField('WSE (ft)', water_level[i])
                feature.SetField('D (ft)', dist[i])
                
                layer.CreateFeature(feature)
        
        
        #### loop to write data for each COC type ####
        for i in range(len(self.contaminant_names[:])):
            
            contaminant_name = self.contaminant_names[i]
            solubility = self.solubilities[i]
            density = self.densities[i]
            
            print ("Writing %i: %s\n"%(i, contaminant_name)) 
            
            if solubility == 'Soluble':
                verbose = 'soluble'
            elif solubility == 'Insoluble' and density == 'Float':
                verbose = 'insoluble_float'
            elif solubility == 'Insoluble' and density == 'Sink':
                verbose = 'insoluble_sink'

            branchIDs_combined, SegIDs_combined, lines_combined, \
            lats_combined, lons_combined, \
            Ttimes_combined, Density_combined, Initial_loc_combined, solubility_combined, \
            flow_combined, concentration_combined, water_level_combined, dist_combined \
                = self.readcategory(verbose)
            
            add_feature(layer, contaminant_name, branchIDs_combined, SegIDs_combined, lines_combined, \
                    lons_combined, lats_combined,\
                    Ttimes_combined, Density_combined, Initial_loc_combined, solubility_combined, \
                    flow_combined, concentration_combined, water_level_combined, dist_combined)
            
            
            
        
    def read_traveltime(self):
        """
        read travel times
        """
        
        basedir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\constant_flow_scenarios\W2\spillway_150\excel_test'
        
        #### Particle travel time branch 1
        excelfile_surface_branch1_scenario1 = r'%s\particle_surface_branch1_scenario1.xlsx'%basedir
        inarray_surface_branch1_scenario1 = pd.read_excel(excelfile_surface_branch1_scenario1).to_numpy()    
        
        excelfile_surface_branch1_scenario2 = r'%s\particle_surface_branch1_scenario2.xlsx'%basedir
        inarray_surface_branch1_scenario2 = pd.read_excel(excelfile_surface_branch1_scenario2).to_numpy()    
        
        excelfile_surface_branch1_scenario3 = r'%s\particle_surface_branch1_scenario3.xlsx'%basedir
        inarray_surface_branch1_scenario3 = pd.read_excel(excelfile_surface_branch1_scenario3).to_numpy()
        
        excelfile_surface_branch1_scenario4 = r'%s\particle_surface_branch1_scenario4.xlsx'%basedir
        inarray_surface_branch1_scenario4 = pd.read_excel(excelfile_surface_branch1_scenario4).to_numpy()
        
        
#        excelfile_bottom_branch1_high = r'excel\flow_rate\particle_bottom_branch1_high.xlsx'
#        inarray_bottom_branch1_high = pd.read_excel(excelfile_bottom_branch1_high).to_numpy()
#        
#        excelfile_bottom_branch1_medium = r'excel\flow_rate\particle_bottom_branch1_medium.xlsx'
#        inarray_bottom_branch1_medium = pd.read_excel(excelfile_bottom_branch1_medium).to_numpy()
#        
#        excelfile_bottom_branch1_low = r'excel\flow_rate\particle_bottom_branch1_low.xlsx'
#        inarray_bottom_branch1_low = pd.read_excel(excelfile_bottom_branch1_low).to_numpy()
        
        
        #### Tracer travel time branch 1
        excelfile_tracer_branch1_scenario1 = r'%s\tracer_branch1_scenario1.xlsx'%basedir
        inarray_tracer_branch1_scenario1 = pd.read_excel(excelfile_tracer_branch1_scenario1).to_numpy()
        
        excelfile_tracer_branch1_scenario2 = r'%s\tracer_branch1_scenario2.xlsx'%basedir
        inarray_tracer_branch1_scenario2 = pd.read_excel(excelfile_tracer_branch1_scenario2).to_numpy()
        
        excelfile_tracer_branch1_scenario3 = r'%s\tracer_branch1_scenario3.xlsx'%basedir
        inarray_tracer_branch1_scenario3 = pd.read_excel(excelfile_tracer_branch1_scenario3).to_numpy()
        
        excelfile_tracer_branch1_scenario4 = r'%s\tracer_branch1_scenario4.xlsx'%basedir
        inarray_tracer_branch1_scenario4 = pd.read_excel(excelfile_tracer_branch1_scenario4).to_numpy()
        
        
        self.inarrays_branch1 = [inarray_surface_branch1_scenario1, inarray_surface_branch1_scenario2, inarray_surface_branch1_scenario3, inarray_surface_branch1_scenario4, \
                                 #         inarray_bottom_branch1_high, inarray_bottom_branch1_medium, inarray_bottom_branch1_low, \
                                          inarray_tracer_branch1_scenario1, inarray_tracer_branch1_scenario2, inarray_tracer_branch1_scenario3, inarray_tracer_branch1_scenario4]
        
        
        #### Particle travel time branch 5
        excelfile_surface_branch5_scenario1 = r'%s\particle_surface_branch5_scenario1.xlsx'%basedir
        inarray_surface_branch5_scenario1 = pd.read_excel(excelfile_surface_branch5_scenario1).to_numpy()
        
        excelfile_surface_branch5_scenario2 = r'%s\particle_surface_branch5_scenario2.xlsx'%basedir
        inarray_surface_branch5_scenario2 = pd.read_excel(excelfile_surface_branch5_scenario2).to_numpy()
        
        excelfile_surface_branch5_scenario3 = r'%s\particle_surface_branch5_scenario3.xlsx'%basedir
        inarray_surface_branch5_scenario3 = pd.read_excel(excelfile_surface_branch5_scenario3).to_numpy()
        
        excelfile_surface_branch5_scenario4 = r'%s\particle_surface_branch5_scenario4.xlsx'%basedir
        inarray_surface_branch5_scenario4 = pd.read_excel(excelfile_surface_branch5_scenario4).to_numpy()
        
        
#        excelfile_bottom_branch5_high = r'excel\flow_rate\particle_bottom_branch5_high.xlsx'
#        inarray_bottom_branch5_high = pd.read_excel(excelfile_bottom_branch5_high).to_numpy()
#        
#        excelfile_bottom_branch5_medium = r'excel\flow_rate\particle_bottom_branch5_medium.xlsx'
#        inarray_bottom_branch5_medium = pd.read_excel(excelfile_bottom_branch5_medium).to_numpy()
#        
#        excelfile_bottom_branch5_low = r'excel\flow_rate\particle_bottom_branch5_low.xlsx'
#        inarray_bottom_branch5_low = pd.read_excel(excelfile_bottom_branch5_low).to_numpy()
        
        
        #### Tracer travel time branch 5
        excelfile_tracer_branch5_scenario1 = r'%s\tracer_branch5_scenario1.xlsx'%basedir
        inarray_tracer_branch5_scenario1 = pd.read_excel(excelfile_tracer_branch5_scenario1).to_numpy()
        
        excelfile_tracer_branch5_scenario2 = r'%s\tracer_branch5_scenario2.xlsx'%basedir
        inarray_tracer_branch5_scenario2 = pd.read_excel(excelfile_tracer_branch5_scenario2).to_numpy()
        
        excelfile_tracer_branch5_scenario3 = r'%s\tracer_branch5_scenario3.xlsx'%basedir
        inarray_tracer_branch5_scenario3 = pd.read_excel(excelfile_tracer_branch5_scenario3).to_numpy()
        
        excelfile_tracer_branch5_scenario4 = r'%s\tracer_branch5_scenario4.xlsx'%basedir
        inarray_tracer_branch5_scenario4 = pd.read_excel(excelfile_tracer_branch5_scenario4).to_numpy()
        
        
        self.inarrays_branch5 = [inarray_surface_branch5_scenario1, inarray_surface_branch5_scenario2, inarray_surface_branch5_scenario3, inarray_surface_branch5_scenario4, \
                                          #inarray_bottom_branch5_high, inarray_bottom_branch5_medium, inarray_bottom_branch5_low, \
                                          inarray_tracer_branch5_scenario1, inarray_tracer_branch5_scenario2, inarray_tracer_branch5_scenario3, inarray_tracer_branch5_scenario4]
        
        
        
        
    def readcategory(self, verbose):
        """
        verbose: soluble, insoluble_float, insoluble_sink
        """
        
        #### coordinate information 
        westlats1 = []
        westlons1 = []
        eastlats1 = []
        eastlons1 = []        
        lines1 = []
        #### points ####
        lats1 = []
        lons1 = []
        for i in range(len(self.westPnts1)):
            westlat1, westlon1 = utm.to_latlon(self.westPnts1[i,0], self.westPnts1[i,1], 14, 'U')
            eastlat1, eastlon1 = utm.to_latlon(self.eastPnts1[i,0], self.eastPnts1[i,1], 14, 'U')
            lines1.append([[westlon1, westlat1], [eastlon1, eastlat1]])
            westlats1.append(westlat1)
            westlons1.append(westlon1)
            eastlats1.append(eastlat1)
            eastlons1.append(eastlon1)
            #### points ####
            lats1.append((westlat1+eastlat1)/2.)
            lons1.append((westlon1+eastlon1)/2.)
            
        westlats5 = []
        westlons5 = []
        eastlats5 = []
        eastlons5 = []        
        lines5 = []
        #### points ####
        lats5 = []
        lons5 = []
        for i in range(len(self.westPnts5)):
            westlat5, westlon5 = utm.to_latlon(self.westPnts5[i,0], self.westPnts5[i,1], 14, 'U')
            eastlat5, eastlon5 = utm.to_latlon(self.eastPnts5[i,0], self.eastPnts5[i,1], 14, 'U')
            lines5.append([[westlon5, westlat5], [eastlon5, eastlat5]])
            westlats5.append(westlat5)
            westlons5.append(westlon5)
            eastlats5.append(eastlat5)
            eastlons5.append(eastlon5)
            #### points ####
            lats5.append((westlat5+eastlat5)/2.)
            lons5.append((westlon5+eastlon5)/2.)
            
        
        self.read_traveltime()
        
        
        
        if verbose == "soluble":
            
            inarrays_branch1 = self.inarrays_branch1[0:4]
            inarrays_branch5 = self.inarrays_branch5[0:4]
            
        elif verbose == "insoluble_float":
            
            inarrays_branch1 = self.inarrays_branch1[4:8]
            inarrays_branch5 = self.inarrays_branch5[4:8]
            
        elif verbose == "insoluble_sink":
            
            print ("Not available: wait for simulation done")
            inarrays_branch1 = self.inarrays_branch1[4:8]
            inarrays_branch5 = self.inarrays_branch5[4:8]
        
        
        Narray_branch1 = len(inarrays_branch1)
        Narray_branch5 = len(inarrays_branch5)
        
        
        #### travel time for branch 1
        Ttime = inarrays_branch1[0][:,2]
        ind0 = np.nonzero(Ttime)[0][0]
        ind = np.arange(ind0, Ttime.shape[0])
        
        #### soluble
        branchIDs_branch1 = []
        SegIDs_branch1 = []
        lines_branch1 = []
        westlats_branch1 = []
        westlons_branch1 = []
        eastlats_branch1 = []
        eastlons_branch1 = []
        lats_branch1 = []
        lons_branch1 = []
        Ttimes_branch1 = []
        Density_branch1 = []
        Initial_loc_branch1 = []
        solubility_branch1 = []
        flow_condition_branch1 = []
        concentration_branch1 = []
        water_level_branch1 = []
        dist_branch1 = []
        
        
        for iarray in range(Narray_branch1):
            
            #### find indexes which segment has travel time
            Ttime_tem = inarrays_branch1[iarray][:,2]
            ind0 = np.nonzero(Ttime_tem)[0][0]
            ind = np.arange(ind0, Ttime_tem.shape[0])
            
            for i in range(inarrays_branch1[0].shape[0]):
                
                if i in ind:
                    branchIDs_branch1.append(inarrays_branch1[iarray][i,0])
                    SegIDs_branch1.append(inarrays_branch1[iarray][i,1])
                    lines_branch1.append(lines1[i])
                    westlats_branch1.append(westlats1[i])
                    westlons_branch1.append(westlons1[i])
                    eastlats_branch1.append(eastlats1[i])
                    eastlons_branch1.append(eastlons1[i])
                    lats_branch1.append(lats1[i])
                    lons_branch1.append(lons1[i])
                    
                    Ttimes_branch1.append(inarrays_branch1[iarray][i,2])
                    if inarrays_branch1[iarray][i,3] == 0:
                        Density_branch1.append('Light')
                    elif inarrays_branch1[iarray][i,3] == 1:
                        Density_branch1.append('Heavy')
                    elif inarrays_branch1[iarray][i,3] == 9:
                        Density_branch1.append('None')
                        
                    Initial_loc_branch1.append('East')
                    
                    if inarrays_branch1[iarray][i,5] == 0:
                        solubility_branch1.append('Insoluble')
                    elif inarrays_branch1[iarray][i,5] == 1:
                        solubility_branch1.append('Soluble')
                        
                    if inarrays_branch1[iarray][i,6] == 1:
                        flow_condition_branch1.append('100 cfs')
                    elif inarrays_branch1[iarray][i,6] == 2:
                        flow_condition_branch1.append('1000 cfs')
                    elif inarrays_branch1[iarray][i,6] == 3:
                        flow_condition_branch1.append('2000 cfs')
                    elif inarrays_branch1[iarray][i,6] == 4:
                        flow_condition_branch1.append('3000 cfs')
                    
                    if inarrays_branch1[iarray][i,7] != 0:
                        concentration_branch1.append("{:.3E}".format(Decimal(inarrays_branch1[iarray][i,7])))
                    else:
                        concentration_branch1.append(str(inarrays_branch1[iarray][i,7]))
                        
                    water_level_branch1.append(inarrays_branch1[iarray][i,8])
                    dist_branch1.append(inarrays_branch1[iarray][i,9])
                    
        branchIDs_branch5 = []
        SegIDs_branch5 = []
        lines_branch5 = []
        westlats_branch5 = []
        westlons_branch5 = []
        eastlats_branch5 = []
        eastlons_branch5 = []
        lats_branch5 = []
        lons_branch5 = []
        Ttimes_branch5 = []
        Density_branch5 = []
        Initial_loc_branch5 = []
        solubility_branch5 = []
        flow_condition_branch5 = []
        concentration_branch5 = []
        water_level_branch5 = []
        dist_branch5 = []  
        
        ## loop over all travel time for each array, find which is in branch 1 and which is in branch 5
        for iarray in range(Narray_branch5): 
            
            #### find indexes which segment has travel time
            Ttime_tem = inarrays_branch5[iarray][:,2]
            
            nbr5 = len(lines5)  ## number of elements in branch 5
            ind1 = np.arange(43, 45) -1 + nbr5    #### hard coded, for release in branch 5 len(branch5)+ [43,44] - 1
            ind5 = np.nonzero(Ttime_tem)[0]
            
            for i in range(inarrays_branch5[0].shape[0]):
                #if iarray==6 and i == 44:
                #    pdb.set_trace()
                if inarrays_branch5[iarray][i,0] == 5: ## at branch 5
                    
                    if i in ind5:
                        branchIDs_branch5.append(inarrays_branch5[iarray][i,0])
                        SegIDs_branch5.append(inarrays_branch5[iarray][i,1])
                        lines_branch5.append(lines5[i])
                        westlats_branch5.append(westlats5[i])
                        westlons_branch5.append(westlons5[i])
                        eastlats_branch5.append(eastlats5[i])
                        eastlons_branch5.append(eastlons5[i])
                        lats_branch5.append(lats5[i])
                        lons_branch5.append(lons5[i])
                        
                        Ttimes_branch5.append(inarrays_branch5[iarray][i,2])
                        if inarrays_branch5[iarray][i,3] == 0:
                            Density_branch5.append('Light')
                        elif inarrays_branch5[iarray][i,3] == 1:
                            Density_branch5.append('Heavy')
                        elif inarrays_branch5[iarray][i,3] == 9:
                            Density_branch5.append('None')
                            
                        if inarrays_branch5[iarray][i,4] == 1:
                            Initial_loc_branch5.append('East')
                        elif inarrays_branch5[iarray][i,4] == 5:
                            Initial_loc_branch5.append('West')
                        
                        if inarrays_branch5[iarray][i,5] == 0:
                            solubility_branch5.append('Insoluble')
                        elif inarrays_branch5[iarray][i,5] == 1:
                            solubility_branch5.append('Soluble')
                            
                        if inarrays_branch5[iarray][i,6] == 1:
                            flow_condition_branch5.append('100 cfs')
                        elif inarrays_branch5[iarray][i,6] == 2:
                            flow_condition_branch5.append('1000 cfs')
                        elif inarrays_branch5[iarray][i,6] == 3:
                            flow_condition_branch5.append('2000 cfs')
                        elif inarrays_branch5[iarray][i,6] == 4:
                            flow_condition_branch5.append('3000 cfs')
                        
                        if inarrays_branch5[iarray][i,7] != 0:
                            concentration_branch5.append("{:.3E}".format(Decimal(inarrays_branch5[iarray][i,7])))
                        else: 
                            concentration_branch5.append(str(inarrays_branch5[iarray][i,7]))
                            
                        water_level_branch5.append(inarrays_branch5[iarray][i,8])
                        dist_branch5.append(inarrays_branch5[iarray][i,9])
                
                
                elif inarrays_branch5[iarray][i,0] == 1: ## at branch 1
                    
                    if i in ind1:
                        
                        branchIDs_branch5.append(inarrays_branch5[iarray][i,0])
                        SegIDs_branch5.append(inarrays_branch5[iarray][i,1])
                        lines_branch5.append(lines1[i-nbr5])
                        westlats_branch5.append(westlats1[i-nbr5])
                        westlons_branch5.append(westlons1[i-nbr5])
                        eastlats_branch5.append(eastlats1[i-nbr5])
                        eastlons_branch5.append(eastlons1[i-nbr5])
                        lats_branch5.append(lats1[i-nbr5])
                        lons_branch5.append(lons1[i-nbr5])
                        
                        Ttimes_branch5.append(inarrays_branch5[iarray][i,2])
                        if inarrays_branch5[iarray][i,3] == 0:
                            Density_branch5.append('Light')
                        elif inarrays_branch5[iarray][i,3] == 1:
                            Density_branch5.append('Heavy')
                        elif inarrays_branch5[iarray][i,3] == 9:
                            Density_branch5.append('None')
                            
                        if inarrays_branch5[iarray][i,4] == 1:
                            Initial_loc_branch5.append('East')
                        elif inarrays_branch5[iarray][i,4] == 5:
                            Initial_loc_branch5.append('West')
                            
                        if inarrays_branch5[iarray][i,5] == 0:
                            solubility_branch5.append('Insoluble')
                        elif inarrays_branch5[iarray][i,5] == 1:
                            solubility_branch5.append('Soluble')
                            
                        if inarrays_branch5[iarray][i,6] == 1:
                            flow_condition_branch5.append('100 cfs')
                        elif inarrays_branch5[iarray][i,6] == 2:
                            flow_condition_branch5.append('1000 cfs')
                        elif inarrays_branch5[iarray][i,6] == 3:
                            flow_condition_branch5.append('2000 cfs')
                        elif inarrays_branch5[iarray][i,6] == 4:
                            flow_condition_branch5.append('3000 cfs')
                            
                        if inarrays_branch5[iarray][i,7] != 0:
                            concentration_branch5.append("{:.3E}".format(Decimal(inarrays_branch5[iarray][i,7])))
                        else:
                            concentration_branch5.append(str(inarrays_branch5[iarray][i,7]))
                            
                        water_level_branch5.append(inarrays_branch5[iarray][i,8])
                        dist_branch5.append(inarrays_branch5[iarray][i,9])
                        
                        
        #### combine all data into one big array
        branchIDs_combined = branchIDs_branch1 + branchIDs_branch5
        SegIDs_combined = SegIDs_branch1 + SegIDs_branch5 
        lines_combined = lines_branch1 + lines_branch5

        lats_combined = lats_branch1 + lats_branch5
        lons_combined = lons_branch1 + lons_branch5
        
        Ttimes_combined = Ttimes_branch1 + Ttimes_branch5
        Density_combined = Density_branch1 + Density_branch5
        Initial_loc_combined = Initial_loc_branch1 + Initial_loc_branch5
        solubility_combined = solubility_branch1 + solubility_branch5
        flow_combined = flow_condition_branch1 + flow_condition_branch5
        concentration_combined = concentration_branch1 + concentration_branch5 
        water_level_combined = water_level_branch1 + water_level_branch5 
        dist_combined = dist_branch1 + dist_branch5
        
        return branchIDs_combined, SegIDs_combined, lines_combined, \
                lats_combined, lons_combined, \
                Ttimes_combined, Density_combined, Initial_loc_combined, solubility_combined, \
                flow_combined, concentration_combined, water_level_combined, dist_combined

if __name__ == "__main__": 
    ## any directory that has the bathymetry file
    wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\constant_flow_scenarios\W2\spillway_150\scenarios_test_new\scenario1'
    ArcMAP = ArcGIS_online_map_COC(wdir)
    ArcMAP.write_shapefile_combined(shpname='Travel_time_COC.shp')