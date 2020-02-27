# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 13:53:40 2019

@author: dfeng
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib as mpl

from w2_contour import W2_Contour
from bathymetry import W2_Bathymetry
from segmentation import W2_Segmentation

import pdb


class Tracer_Travel_Time(W2_Contour):
    
    #branch2:
    DHS2 = 33
    #branch3:
    DHS3 = 34
    #bracnh4
    DHS4 = 42
    #branch5
    DHS5 = 43
    
    flows = {'high': 3, 'medium': 2, 'low': 1} 
    

    def __init__(self, workdir, **kwargs):
        
        self.workdir = workdir
        
    
    
    def travel_time_full_branch(self, starttime, initialBranch=1, initialSeg=1, flow_condition='high', write2shp=False, write2txt=True):
        """
        call travel time function, and calculate the travel time for all branches: 1 to 5
        """
        self.initialBranch = initialBranch
        self.initialSeg = initialSeg
        
        #### read conservative tracer file for all branches, so don't need to read the file multiple times
        self.Readcpl()
        
        #starttime = 0
        #endtime = 149
        #starttime = 334
        #endtime = 500
        starttime =  starttime
        endtime = starttime + 300
        
        #### calculate the travel time for all branches  endtime=299
        Ttime1 = self.travel_time(starttime=starttime, endtime=endtime, branchID=1)
        #Ttime2 = self.travel_time(endtime=endtime, branchID=2)
        #Ttime3 = self.travel_time(endtime=endtime, branchID=3)
        #Ttime4 = self.travel_time(endtime=endtime, branchID=4)
        Ttime5 = self.travel_time(starttime=starttime, endtime=endtime, branchID=5)
        
        
        if self.initialBranch == 5:
            
            print ('Spill released at branch 5 \n')
            #Ttime1[self.DHS5-2:] += Ttime5[-1]
            #### because branch 1 and branch 5 are calculated separate
            if Ttime1[self.DHS5] < Ttime5[-2]:
                dt_tem = Ttime1[self.DHS5+1] - Ttime1[self.DHS5]
                Ttime1[self.DHS5] = Ttime5[-2] + 1
                Ttime1[self.DHS5+1] = Ttime1[self.DHS5] + dt_tem
                
        elif self.initialBranch == 1:
            
            print ('Spill released at branch 1 \n')
        
        
        ## delete inactive segments
        Ttime1 = Ttime1[1:-1] 
        #Ttime2 = Ttime2[1:-1] 
        #Ttime3 = Ttime3[1:-1] 
        #Ttime4 = Ttime4[1:-1] 
        Ttime5 = Ttime5[1:-1]

        ## branch 2, 3, 4, 5, adding time from the main branch
        #Ttime2[Ttime2!=0] = Ttime2[Ttime2!=0] - Ttime2[-1] + Ttime1[self.DHS2-2]
        #Ttime3[Ttime3!=0] = Ttime3[Ttime3!=0] - Ttime3[-1] + Ttime1[self.DHS3-2]
        #Ttime4[Ttime4!=0] = Ttime4[Ttime4!=0] - Ttime4[-1] + Ttime1[self.DHS4-2]
        #Ttime5[Ttime5!=0] = Ttime5[Ttime5!=0] - Ttime5[-1] + Ttime1[self.DHS5-2]
        
        
        ## also remember to remove inactivate segments from Ttimes
        Ttimes = [Ttime1, Ttime5]

        ## save to txt file
        if write2txt:
            #### call segment class for segment information
            Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
            WS = W2_Segmentation(Bthfile)
            WS.VisSeg2()
            
            #### save travel time data to txt file ####
            if self.initialBranch == 1:
                
                ## conversion to zero travel time at donwstream gate (only on nonzero values)
                Ttime = Ttimes[0]
                Ttime[Ttime!=0] = Ttime[-1] - Ttime[Ttime!=0]
                
                ## save txt
                txtfile=r'txt\tracer_branch%s_%s.txt'%(str(self.initialBranch), flow_condition)
                density = 9    ## density information not useful for non-soluble contanminants
                self.savetxt_Traveltime_branch1(WS, Ttime, density, self.flows[flow_condition], txtfile)
                
                
            elif self.initialBranch == 5:
                
                ## conversion to zero travel time at donwstream gate (only on nonzero values)
                MaxTime = Ttimes[0][-1]
                for Ttime in Ttimes:
                    Ttime[Ttime!=0] = MaxTime - Ttime[Ttime!=0]
                
                ## save txt
                txtfile=r'txt\tracer_branch%s_%s.txt'%(str(self.initialBranch), flow_condition)
                density = 9
                self.savetxt_Traveltime_branch5(WS, Ttimes, density, self.flows[flow_condition], txtfile)

                
        
    
        self.Ttime_plot(Ttimes, plotBranch=True, write2shp=write2shp)    
        
        
    def travel_time(self, starttime=0, endtime=1460, branchID=1):
        """
        calculate the travel time (or median arrival time) for conservative tracers
        timestep - when to calculate the travel time 
        branchID - the travel time of which branch
        """
        print ("Calculate travel time for branch %s ... \n"%str(branchID))
        
        #### read bathymetry information
        Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
        WB = W2_Bathymetry(Bthfile)
        pat = WB.VisBranch2(branchID)
        
        #### create empty array for travel time
        Ttime = np.zeros_like(WB.X)
        
        #### calculate travel time ####
        #### loop over all time steps ####
        #### starttime to endtime, during which period to search for tracer travel time
        for tstep in range(starttime, endtime): 
            print ('Time step = %s \n'%str(tstep))
            
            ## search for index for each branch 
            ## algorithm find the distance between two elements in self.X_flow that are large, 
            ## this is where the two branches separate
            dist = np.diff(self.X_flow[tstep])
            inds = np.where(dist>1200)[0]
            if branchID == 1:
                ind0 = 0
                ind1 = inds[0]
            elif branchID == 2:
                ind0 = inds[0]+1
                ind1 = inds[1]
            elif branchID == 3:
                ind0 = inds[1]+1
                ind1 = inds[2]
            elif branchID == 4:
                ind0 = inds[2]+1
                ind1 = inds[3]
            elif branchID == 5:
                ind0 = inds[3]+1
                ind1 = len(self.X_flow[tstep]) - 1   ## -1 remove the array size mismatch issue

            ## tracer locations
            X_flow = self.X_flow[tstep][ind0:ind1+1]
            Z_flow = self.Z_flow[tstep][ind0:ind1+1]
            X_flow = np.asarray(X_flow)
            Z_flow = np.asarray(Z_flow)
            
            ## align coordinates with the grid
            dx =  WB.X.max() - X_flow.max()
            X_flow += dx
            
            vartem = np.asarray( self.var_output['Tracer']['value'][tstep][ind0:ind1+1] )
            
            #### quality control if X_flow, vartem not in the same shape, resize
            if X_flow.shape != vartem.shape:
                #pdb.set_trace()
                Lmin = np.min([X_flow.shape[0], vartem.shape[0]])
                if X_flow.shape[0] > vartem.shape[0]:
                    X_flow = np.delete(X_flow, np.arange(Lmin, X_flow.shape[0]))
                elif X_flow.shape[0] < vartem.shape[0]:
                    vartem = np.delete(vartem, np.arange(Lmin, vartem.shape[0]))
                    
            
            ## search nonzero X coordinates
            X_nonzero = self.Xsearch(X_flow, Z_flow, vartem)  
            ## search nonzero X index
            ind_nonzero = self.GridSearch(WB.X, X_nonzero, branchID)  ## output is segID
            
            
            if ind_nonzero.size != 0:     #### there is nonzero tracer in the branch
                
                if branchID in [1,5]:
                    
                    #pdb.set_trace()
                    #### if tracer is released at day 385.9 at segment 18, the tracer would travel to segment 22 at day 386
                    #### so it is important to record the initial segment in the travel time array
                    if tstep == starttime and starttime >= 1:
                        #Ttime[np.min(ind_nonzero)-1] = tstep - 1
                        if branchID == 1:
                            Ttime[self.initialSeg-1] = tstep - 1
                        elif branchID == 5:
                            Ttime[self.initialSeg-85] = tstep - 1
                            
                    
                    ind_max = np.max(ind_nonzero)
                    ## add the travel time to the pre-defined 1D array
                    if Ttime[ind_max] == 0:
                        if len(Ttime.nonzero()[0]) != 0:
                            ind_tem = np.where(Ttime==Ttime.max())[0][-1]
                            Ttime[ind_tem+1: ind_max] = tstep
                        else:
                            Ttime[ind_max] = tstep
                            
                            
                elif branchID in [2,3,4]:
                    if len(Ttime.nonzero()[0])==0 and len(ind_nonzero) == 1:
                        Ttime[ind_nonzero[0]] = tstep
                    else:
                        ind_min = np.min(ind_nonzero)
                        ## add the travel time to the pre-defined 1D array
                        if Ttime[ind_min] == 0:
                            if len(Ttime.nonzero()[0]) != 0:
                                ind_tem = np.where(Ttime==Ttime.max())[0][0]
                                Ttime[ind_min:ind_tem] = tstep
                            else:
                                Ttime[ind_min] = tstep
            print (Ttime)
            
        if branchID in [1, 5]:
            Ttime[-2] = Ttime[-3] + 1
        
        return Ttime
    
    
    
    def Ttime_plot(self, Ttimes, plotBranch=True, write2shp=False):
        """
        Visualize the travel time on W2 segments
        multiple lines' color changing the values corresponding to gradient color:
            https://stackoverflow.com/questions/55587674/how-to-make-multiple-lines-color-changing-with-values-corresponding-to-gradient
        Matplotlib - add colorbar to a sequence of line plots:
            https://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots
        """
        
        #### call segment class for plotting
        Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
        WS = W2_Segmentation(Bthfile)
        WS.VisSeg2()
        
        
        if self.initialBranch == 1: 
            ## branch 1
            
            MaxTime = Ttimes[0][-1]
            Ttime = Ttimes[0]
            ## conversion to zero travel time at donwstream gate (only on nonzero values)
            Ttime[Ttime!=0] = Ttime[-1] - Ttime[Ttime!=0]
        
            #### write to ArcGIS
            if write2shp:
                self.writeShpLines_one_branch(WS, Ttime, shpname='traveltime_branch1')
        
            #### line color
            zz = np.asarray(Ttime)
            r = (zz.astype(np.float)-zz.min())/(zz.max()-zz.min())
            g = np.zeros_like(r)
            b = 1 - r
            
            colorlist = []
            for i in range(len(zz)):
                colorlist.append((r[i], g[i], b[i]))
            #### gray color for zero travel time
            for i in range(len(zz)):
                if zz[i] == 0:
                    colorlist[i] = 'gray'
                    
            colorlist = []
            for i in range(len(zz)):
                colorlist.append((r[i], g[i], b[i]))
            #### gray color for zero travel time
            for i in range(len(zz)):
                if zz[i] == 0:
                    colorlist[i] = 'gray'
            #### colorbar: Setting up a colormap that's a simple transtion
            mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
            # Using contourf to provide my colorbar info, then clearing the figure
            ss = 1
            levels = np.arange(zz.min(),zz.max(),ss)
            CS3 = plt.contourf([[0,0],[0,0]], levels, cmap=mymap)
            plt.clf()
            
            #### matplotlib plotting
            plt.rcParams.update({'font.size': 14})
            fig = plt.figure(figsize=(10,12))
            ax = fig.add_subplot(111)
            
            if plotBranch:
                ax.plot(WS.Pnts1[:,0], WS.Pnts1[:,1], '.-', markersize=2, color='0.5')
            
            #### Plot segments
            for i in range(len(WS.westPnts1)):
                ax.plot([WS.westPnts1[i,0], WS.eastPnts1[i,0]], [WS.westPnts1[i,1], WS.eastPnts1[i,1]], '-', color=colorlist[i])
            for i in range(len(WS.segs1)):
                if Ttimes[0][i] != 0:
                    ax.annotate('%s'%str(int(Ttimes[0][i])), (WS.Pnts1[i,0], WS.Pnts1[i,1]), color=colorlist[i], fontsize=8)
        
            
        elif self.initialBranch == 5:
            ## branch 5
            
            nt = []
            for l in Ttimes:
                nt.append(len(l))
            
            MaxTime = Ttimes[0][-1]
            for Ttime in Ttimes:
                Ttime[Ttime!=0] = MaxTime - Ttime[Ttime!=0]
            
            #### write to ArcGIS
            if write2shp:
                self.writeShpLines(WS, Ttimes, shpname='traveltime_branch5')
            
            #### line color
            zz = np.asarray(Ttimes[0].tolist() + Ttimes[1].tolist())
            r = (zz.astype(np.float)-zz.min())/(zz.max()-zz.min())
            g = np.zeros_like(r)
            b = 1 - r    
        
            colorlist = []
            for i in range(len(zz)):
                colorlist.append((r[i], g[i], b[i]))
            #### gray color for zero travel time
            for i in range(len(zz)):
                if zz[i] == 0:
                    colorlist[i] = 'gray'
                    
            #### colorbar: Setting up a colormap that's a simple transtion
            mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
            # Using contourf to provide my colorbar info, then clearing the figure
            ss = 1
            levels = np.arange(zz.min(),zz.max(),ss)
            CS3 = plt.contourf([[0,0],[0,0]], levels, cmap=mymap)
            plt.clf()
            
            #### matplotlib plotting
            plt.rcParams.update({'font.size': 14})
            fig = plt.figure(figsize=(10,12))
            ax = fig.add_subplot(111)
            
            if plotBranch:
                ax.plot(WS.Pnts1[:,0], WS.Pnts1[:,1], '.-', markersize=2, color='0.5')
                ax.plot(WS.Pnts5[:,0], WS.Pnts5[:,1], '.-', markersize=2, color='0.5')
            
            for i in range(len(WS.westPnts1)):
                ax.plot([WS.westPnts1[i,0], WS.eastPnts1[i,0]], [WS.westPnts1[i,1], WS.eastPnts1[i,1]], '-', color=colorlist[i])
            for i in range(len(WS.westPnts5)):
                ax.plot([WS.westPnts5[i,0], WS.eastPnts5[i,0]], [WS.westPnts5[i,1], WS.eastPnts5[i,1]], '-', color=colorlist[i+nt[0]])

            #### Plot segments               
            for i in range(len(WS.segs1)):
                if Ttimes[0][i] != 0:
                    ax.annotate('%s'%str(int(Ttimes[0][i])), (WS.Pnts1[i,0], WS.Pnts1[i,1]), color=colorlist[i], fontsize=8)
            for i in range(len(WS.segs5)):
                if Ttimes[1][::-1][i] != 0:
                    ax.annotate('%s'%str(int(Ttimes[1][::-1][i])), (WS.Pnts5[i,0], WS.Pnts5[i,1]), color=colorlist[sum(nt[:])-i-1], fontsize=8)

        ax.set_xlabel('Easting [m]')
        ax.set_ylabel('Northing [m]')
        fig.colorbar(CS3)
        ax.set_aspect(True)
        fig.tight_layout()
        plt.show()
        

        
            
    def Xsearch(self, xin, zin, varin):
        """
        A searching function, 
        Input X and variable
        Output X with nonzero value (Note the masked value -99)
        Note: do not use the non-zero value at the surface layer
        """
        
        varin = np.ma.masked_array(varin,mask=varin==self.mask_value)
        
        ####quality control the masked array
        varin[(varin.mask==False)&(varin<1e-15)]=0
        
        ind = varin.nonzero()[0]  ## This nonzero assumption to find tracer time is not practical
        
        ## 50 percentile arrival time
        pp = 0.5  #### here I use 75 percentile, this value is subject to changes  0.25 works
        idx = pp * (len(ind) - 1)    
        idx = int(idx + 0.5) * (-1)
        ind_sub = np.argpartition(varin[ind], idx)[idx:]
        ind = ind[ind_sub]
        
        #pdb.set_trace()
        
        #### check if nonzero values is at the surface layer
        if len(ind) == 2:
            ind_tem = np.argpartition(zin[ind], len(ind)-1)[:len(ind)-1]
            ind = ind[ind_tem]
        elif len(ind) > 2:
            ind_tem = np.argpartition(zin[ind], len(ind)-2)[:len(ind)-2]
            ind = ind[ind_tem]
        
        return xin[ind]
    
    def GridSearch(self, xgrid, xxin, branchID):
        """
        Search the grid index of nonzero tracer values
        Input a list of X coordinates with nonzero values
        Output a list of grid index in the horizontal direction
        """
        
        segID = np.searchsorted(xgrid, xxin)

        return np.unique(segID)
    
    
    
    def savetxt_Traveltime_branch1(self, WS, Ttime, density, flow_index, txtfile):
        """
        save travel time to a txt file
        output array: 
            branchID, segID, travel time, density, release_arm, solubility, flow_condition
        density=1 heavy
        density=0 light
        release_arm=1 East
        flow_index - high:3, medium:2, low:1
        """
        
        outarray = np.vstack((np.ones_like(Ttime), WS.segs1, Ttime, \
                              np.ones_like(Ttime)*density, np.ones_like(Ttime), \
                              np.ones_like(Ttime)*1, np.ones_like(Ttime)*flow_index)).T
        np.savetxt(txtfile, outarray, fmt='%d')
        
        
    def savetxt_Traveltime_branch5(self, WS, Ttimes, density, flow_index, txtfile):
        """
        output array: 
            branchID, segID, travel time, density, release_arm, solubility, flow_condition
        density=1 heavy
        density=0 light
        release_arm=5 West
        """
        
        Ttime1 = Ttimes[0]    ## travel times at branch 1
        Ttime5 = Ttimes[1]    ## travel times at branch 5
        
        outarray1 = np.vstack((np.ones_like(Ttime1), WS.segs1, Ttime1, \
                              np.ones_like(Ttime1)*density, np.ones_like(Ttime1)*5, \
                              np.ones_like(Ttime1), np.ones_like(Ttime1)*flow_index))
        
        #### important !! reverse WS.segs5, from 86 to 121
        outarray5 = np.vstack((np.ones_like(Ttime5)*5, WS.segs5[::-1], Ttime5, \
                              np.ones_like(Ttime5)*density, np.ones_like(Ttime5)*5, \
                              np.ones_like(Ttime5), np.ones_like(Ttime5)*flow_index))
        
        outarray = np.hstack((outarray5, outarray1)).T
        
        np.savetxt(txtfile, outarray, fmt='%d')
        
    
    
    def writeShpLines(self, WS, Ttimes, shpname):
        """
        write segment lines into ArcGIS readable shapefile
        Input segment object
        """
        import shapefile as shp
        import utm
        
        w = shp.Writer(shpname)
        
        w.field('branchID','C')
        w.field('SegID','C')
        w.field('Lon_west','C')  #float - needed for coordinates
        w.field('Lat_west','C') 
        w.field('Lon_east','C')
        w.field('Lat_east','C')
        w.field('Travel_T','C')
        
        #### only plot non-zero segments
        ind1 = np.nonzero(Ttimes[0])[0]
        ind1 = np.append(ind1, ind1[-1]+1)
        
        ind5 = np.nonzero(Ttimes[1])[0]
        
        
        #for i in range(len(WS.westPnts1)):
        for i in ind1:
            westlat, westlon = utm.to_latlon(WS.westPnts1[i,0], WS.westPnts1[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(WS.eastPnts1[i,0], WS.eastPnts1[i,1], 14, 'U')
            w.line([[[westlon, westlat], [eastlon, eastlat]]])
            w.record('branch1', str(WS.segs1[i]), "{:.3f}".format(westlon), "{:.3f}".format(westlat), \
                     "{:.3f}".format(eastlon), "{:.3f}".format(eastlat), "{:.1f}".format(Ttimes[0][i]))
            
            
        #for i in range(len(WS.westPnts5)):
        for i in ind5:
            westlat, westlon = utm.to_latlon(WS.westPnts5[i,0], WS.westPnts5[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(WS.eastPnts5[i,0], WS.eastPnts5[i,1], 14, 'U')
            w.line([[[westlon, westlat], [eastlon, eastlat]]])
            w.record('branch5', str(WS.segs5[::-1][i]), "{:.3f}".format(westlon), "{:.3f}".format(westlat), \
                     "{:.3f}".format(eastlon), "{:.3f}".format(eastlat), "{:.1f}".format(Ttimes[1][i]))
        
        w.close()
        
        prj = open("%s.prj"%shpname, "w") 
        epsg = 'GEOGCS["WGS 84",'
        epsg += 'DATUM["WGS_1984",'
        epsg += 'SPHEROID["WGS 84",6378137,298.257223563]]'
        epsg += ',PRIMEM["Greenwich",0],'
        epsg += 'UNIT["degree",0.0174532925199433]]'
        prj.write(epsg)
        prj.close()
        
        #pdb.set_trace()
        
    def writeShpLines_one_branch(self, WS, Ttime, shpname):
        """
        write travel time shapefile for one branch using GDAL
        When using GDAL, remember to delete previous files when there are bugs
        Note the feature name for travel time can be too long, so we choose Travel_T instead of Travel_time
        """
        import osgeo.ogr
        import utm
        
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
        ind = np.nonzero(Ttime)[0]
        ind = np.append(ind, ind[-1]+1) ## add the index with 0 travel time at the downstream side
        
        #add_feature(1, layer, lines1, WS.segs1, latlon['westlon1'], latlon['westlat1'], latlon['eastlon1'], latlon['eastlat1'], Ttime)
        add_feature(1, layer, np.asarray(lines1)[ind], WS.segs1[ind], 
                    np.asarray(latlon['westlon1'])[ind], np.asarray(latlon['westlat1'])[ind], 
                    np.asarray(latlon['eastlon1'])[ind], np.asarray(latlon['eastlat1'])[ind], Ttime[ind])
            
        
        
if __name__ == "__main__":
    
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191213_1533_tracer_test'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20200103_1326_tracer_test_branch5'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20200211_1013_tracer_test_branch1'
    
    
    #### it is important to specify at which segment the tracer is initially released
    #### because for medium flow case, the tracer is released at segment 17 at day 725.9,
    #### but at day 726, the tracer might be found in segment 15. 
    
    
    #### branch 1
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20200221_0930_tracer_high_branch1'
    #TTT = Tracer_Travel_Time(wdir)
    #TTT.travel_time_full_branch(starttime=385, initialBranch=1, initialSeg=18, flow_condition='high', write2shp=False)
    
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20200221_0952_tracer_medium_branch1'
    #TTT = Tracer_Travel_Time(wdir)
    #TTT.travel_time_full_branch(starttime=725, initialBranch=1, initialSeg=17, flow_condition='medium', write2shp=False)
    
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20200226_1129_tracer_low_branch1'
    #TTT = Tracer_Travel_Time(wdir)
    #TTT.travel_time_full_branch(starttime=1085, initialBranch=1, initialSeg=18, flow_condition='low', write2shp=False)
    
    
    #### branch 5
    wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20200226_1142_tracer_high_branch5'  ## timestep=385
    TTT = Tracer_Travel_Time(wdir)
    TTT.travel_time_full_branch(starttime=385, initialBranch=5, initialSeg=105, flow_condition='high', write2shp=False)
    #TTT.travel_time_full_branch(initialBranch=5, write2shp=False)
    