# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 16:59:54 2019

@author: dfeng
"""

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import pdb

def Contour2Shp(C,outfile):
    """
    Converts a matplotlib contour object to a shapefile
    
    The shapefile has the filed Contour corresponding to the contour levels set
    in matplotlib.pyplot.contour
    
    **The projection is hardwired to WGS84 for now. This needs updating.
    """
    # Convert the contour values to a shapefile
    # see this example:
    #    http://invisibleroads.com/tutorials/gdal-shapefile-points-save.html
    
    import osgeo.ogr
    
    # Create the projection
    spatialReference = osgeo.osr.SpatialReference()
    spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    
    # Create the shape file
    driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    shapeData = driver.CreateDataSource(outfile)
    
    # Create the layer
    layer = shapeData.CreateLayer('Contour', spatialReference, osgeo.ogr.wkbLineString)
    
    # Create a field containing the contour level
    field_def = osgeo.ogr.FieldDefn('Contour', osgeo.ogr.OFTReal)
    layer.CreateField(field_def)
    layerDefinition = layer.GetLayerDefn()
    
    
    # Loop through the contour object to get the coordinates of each layer
    ctr=0
    for coll,lev in zip(C.collections,C.levels):
        #pdb.set_trace()
        coords = coll.get_paths()
        for xyObj in coords:
            ctr+=1
            line = osgeo.ogr.Geometry(osgeo.ogr.wkbLineString)
            # Add points individually to the line
            for xy in xyObj.vertices:
                line.AddPoint_2D(xy[0],xy[1])
            
            # Update the feature with the line data
            featureIndex = ctr
            feature = osgeo.ogr.Feature(layerDefinition)
            feature.SetGeometry(line)
            feature.SetFID(featureIndex)
            feature.SetGeometryDirectly(line)
            # Set the contour level
            feature.SetField('Contour',lev)
            
            layer.CreateFeature(feature)
    
    # Close the shape file
    shapeData.Destroy()
    print ('Complete - shapefile written to:\n      %s'%outfile)


delta = 0.025
x = np.arange(-3.0, 3.0, delta)
y = np.arange(-2.0, 2.0, delta)
X, Y = np.meshgrid(x, y)
Z1 = np.exp(-X**2 - Y**2)
Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
Z = (Z1 - Z2) * 2

fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z)
ax.clabel(CS, inline=1, fontsize=10)
ax.set_title('Simplest default with labels')

Contour2Shp(CS,'Contour.shp')

plt.show()