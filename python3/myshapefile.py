# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 16:31:19 2020

@author: dfeng
"""

import numpy as np
import osgeo.ogr
import utm



def writeShpLines_one_branch(WS, Ttime, shpname):
        """
        write travel time shapefile for one branch using GDAL
        When using GDAL, remember to delete previous files when there are bugs
        Note the feature name for travel time can be too long, so we choose Travel_T instead of Travel_time
        """
        
        
        latlon = dict()
        latlon.setdefault('westlat1', [])
        latlon.setdefault('westlon1', [])
        latlon.setdefault('eastlat1', [])
        latlon.setdefault('eastlon1', [])
        
        lines1 = []
        for i in range(len(WS.westPnts1)):
            westlat, westlon = utm.to_latlon(WS.westPnts1[i,0], WS.westPnts1[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(WS.eastPnts1[i,0], WS.eastPnts1[i,1], 14, 'U')
            lines1.append([[westlon, westlat], [eastlon, eastlat]])
            latlon['westlat1'].append(westlat)
            latlon['westlon1'].append(westlon)
            latlon['eastlat1'].append(eastlat)
            latlon['eastlon1'].append(eastlon)
        
        #pdb.set_trace()
        # Create the projection
        spatialReference = osgeo.osr.SpatialReference()
        spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
        
        # Create the shape file
        outfile = '%s.shp'%shpname
        driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
        shapeData = driver.CreateDataSource(outfile)
        
        # Create the layer
        layer = shapeData.CreateLayer('Contour', spatialReference, osgeo.ogr.wkbLineString)
        layerDefinition = layer.GetLayerDefn()
        
        # Create fields containing segment infos
        field_def = osgeo.ogr.FieldDefn('branchID', osgeo.ogr.OFTInteger)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('SegID', osgeo.ogr.OFTInteger)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('Lon_west', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('Lat_west', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('Lon_east', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('Lat_east', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('Travel_T', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
        
        def add_feature(branchID, layer, lines, segs, westlon, westlat, eastlon, eastlat, Ttime):
            """
            function that adds feature to layer
            """    
            ctr=0
            for i in range(len(lines)):
                ctr+=1
                line = osgeo.ogr.Geometry(osgeo.ogr.wkbLineString)
                # Add points individually to the line
                xy = lines[i]
        
                line.AddPoint_2D(xy[0][0],xy[0][1])
                line.AddPoint_2D(xy[1][0],xy[1][1])
                # Update the feature with the line data
                featureIndex = ctr
                feature = osgeo.ogr.Feature(layerDefinition)
                #feature.SetStyleString("PEN(c:r,w:5px)")   
                feature.SetGeometry(line)
                feature.SetFID(featureIndex)
                feature.SetGeometryDirectly(line)
                
                # Set the attribute table
                feature.SetField('branchID', int(branchID)) 
                feature.SetField('SegID', int(segs[i]))   # convert to int() is necessary, osgeo cannot recognize numpy int32 type
                feature.SetField('Travel_T', "{:.1f}".format(Ttime[i]))
                feature.SetField('Lon_west', "{:.3f}".format(westlon[i]))
                feature.SetField('Lat_west', "{:.3f}".format(westlat[i]))
                feature.SetField('Lon_east', "{:.3f}".format(eastlon[i]))
                feature.SetField('Lat_east', "{:.3f}".format(eastlat[i]))
                
                layer.CreateFeature(feature)
        
        
        
        ## only plot nonzero travel time
        #ind = np.nonzero(Ttime)[0]
        #ind = np.append(ind, ind[-1]+1) ## add the index with 0 travel time at the downstream side
        ind0 = np.nonzero(Ttime)[0][0]
        ind = np.arange(ind0, Ttime.shape[0])
        
        #add_feature(1, layer, lines1, WS.segs1, latlon['westlon1'], latlon['westlat1'], latlon['eastlon1'], latlon['eastlat1'], Ttime)
        add_feature(1, layer, np.asarray(lines1)[ind], WS.segs1[ind], 
                    np.asarray(latlon['westlon1'])[ind], np.asarray(latlon['westlat1'])[ind], 
                    np.asarray(latlon['eastlon1'])[ind], np.asarray(latlon['eastlat1'])[ind], Ttime[ind])