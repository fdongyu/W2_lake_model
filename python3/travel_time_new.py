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

    def __init__(self, workdir, **kwargs):
        
        self.workdir = workdir
        
        
    def travel_time(self, endtime=1460, branchID=1):
        """
        calculate the travel time (or earliest arrival time) for conservative tracers
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
        for tstep in range(endtime): 
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
            ind_nonzero = self.GridSearch(WB.X, X_nonzero, branchID)  
            
            
            if ind_nonzero.size != 0:     #### there is nonzero tracer in the branch
                
                if branchID == 1:
                    ind_max = np.max(ind_nonzero)
                    ## add the travel time to the pre-defined 1D array
                    if Ttime[ind_max] == 0:
                        if len(Ttime.nonzero()[0]) != 0:
                            ind_tem = np.where(Ttime==Ttime.max())[0][-1]
                            Ttime[ind_tem+1: ind_max] = tstep
                        else:
                            Ttime[ind_max] = tstep
                elif branchID in [2,3,4,5]:
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
            #pdb.set_trace()
        
        if branchID == 1:
            Ttime[-2] = Ttime[-3] + 1
        
        return Ttime
        
    
    def travel_time_full_branch(self, write2shp=False):
        """
        call travel time function, and calculate the travel time for all branches: 1 to 5
        """
        
        #### read conservative tracer file for all branches, so don't need to read the file multiple times
        self.Readcpl()
        
        #### calculate the travel time for all branches 
        Ttime1 = self.travel_time(endtime=299, branchID=1)
        Ttime2 = self.travel_time(endtime=299, branchID=2)
        Ttime3 = self.travel_time(endtime=299, branchID=3)
        Ttime4 = self.travel_time(endtime=299, branchID=4)
        Ttime5 = self.travel_time(endtime=299, branchID=5)
        
        ## delete inactive segments
        Ttime1 = Ttime1[1:-1] 
        Ttime2 = Ttime2[1:-1] 
        Ttime3 = Ttime3[1:-1] 
        Ttime4 = Ttime4[1:-1] 
        Ttime5 = Ttime5[1:-1]
        
        ## branch 2, 3, 4, 5, adding time from the main branch
        Ttime2[Ttime2!=0] = Ttime2[Ttime2!=0] - Ttime2[-1] + Ttime1[self.DHS2-2]
        Ttime3[Ttime3!=0] = Ttime3[Ttime3!=0] - Ttime3[-1] + Ttime1[self.DHS3-2]
        Ttime4[Ttime4!=0] = Ttime4[Ttime4!=0] - Ttime4[-1] + Ttime1[self.DHS4-2]
        Ttime5[Ttime5!=0] = Ttime5[Ttime5!=0] - Ttime5[-1] + Ttime1[self.DHS5-2]
        
        ## also remember to remove inactivate segments from Ttimes
        Ttimes = [Ttime1, Ttime2, Ttime3, Ttime4, Ttime5]
        
        #self.Ttime_plot(Ttimes, write2shp=write2shp)
        self.Ttime_plot_contourf(Ttimes, write2shp=write2shp)  ## create contour lines
        
    
    
    def Ttime_plot(self, Ttimes, plotBranch=True, write2shp=False):
        """
        Visualize the travel time on W2 segments
        multiple lines' color changing the values corresponding to gradient color:
            https://stackoverflow.com/questions/55587674/how-to-make-multiple-lines-color-changing-with-values-corresponding-to-gradient
        Matplotlib - add colorbar to a sequence of line plots:
            https://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots
        """
        
        #### calculate travel time first
        #self.travel_time(endtime=1460, branchID=1)
        
        
        #### call segment class for plotting
        Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
        WS = W2_Segmentation(Bthfile)
        WS.VisSeg2()
        
        
        #### create shapefile to integrate with ArcGIS
        if write2shp:
            self.writeShpLines(WS, Ttimes)
        
        pdb.set_trace()
        
        
        nt = []
        for l in Ttimes:
            nt.append(len(l))
            
        #### line color
        zz = Ttimes[0].tolist() + Ttimes[1].tolist() + Ttimes[2].tolist() + Ttimes[3].tolist() + Ttimes[4].tolist()
        zz = np.asarray(zz)
        r = (zz.astype(np.float)-zz.min())/(zz.max()-zz.min())
        g = np.zeros_like(r)
        b = 1 - r
        
        colorlist = []
        for i in range(len(zz)):
            colorlist.append((r[i], g[i], b[i]))
        
#        #### black color for zero travel time
        for i in range(len(zz)):
            if zz[i] == 0:
                colorlist[i] = 'gray'
                
        #### colorbar
        # Setting up a colormap that's a simple transtion
        mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
        # Using contourf to provide my colorbar info, then clearing the figure
        ss = 1
        levels = np.arange(zz.min(),zz.max(),ss)
        CS3 = plt.contourf([[0,0],[0,0]], levels, cmap=mymap)
        plt.clf()
        
        
        plt.rcParams.update({'font.size': 14})
        fig = plt.figure(figsize=(10,12))
        ax = fig.add_subplot(111)
        
        if plotBranch:
            ax.plot(WS.Pnts1[:,0], WS.Pnts1[:,1], '.-', markersize=2, color='0.5')
            ax.plot(WS.Pnts2[:,0], WS.Pnts2[:,1], '.-', markersize=2, color='0.5')
            ax.plot(WS.Pnts3[:,0], WS.Pnts3[:,1], '.-', markersize=2, color='0.5')
            ax.plot(WS.Pnts4[:,0], WS.Pnts4[:,1], '.-', markersize=2, color='0.5')
            ax.plot(WS.Pnts5[:,0], WS.Pnts5[:,1], '.-', markersize=2, color='0.5')
        
        
        #### Plot segments
        for i in range(len(WS.westPnts1)):
            ax.plot([WS.westPnts1[i,0], WS.eastPnts1[i,0]], [WS.westPnts1[i,1], WS.eastPnts1[i,1]], '-', color=colorlist[i])
        for i in range(len(WS.westPnts2)):
            ax.plot([WS.westPnts2[i,0], WS.eastPnts2[i,0]], [WS.westPnts2[i,1], WS.eastPnts2[i,1]], '-', color=colorlist[i+sum(nt[:1])])
        for i in range(len(WS.westPnts3)):
            ax.plot([WS.westPnts3[i,0], WS.eastPnts3[i,0]], [WS.westPnts3[i,1], WS.eastPnts3[i,1]], '-', color=colorlist[i+sum(nt[:2])])
        for i in range(len(WS.westPnts4)):
            ax.plot([WS.westPnts4[i,0], WS.eastPnts4[i,0]], [WS.westPnts4[i,1], WS.eastPnts4[i,1]], '-', color=colorlist[i+sum(nt[:3])])
        for i in range(len(WS.westPnts5)):
            ax.plot([WS.westPnts5[i,0], WS.eastPnts5[i,0]], [WS.westPnts5[i,1], WS.eastPnts5[i,1]], '-', color=colorlist[i+sum(nt[:4])])
        
        
        #### len(WS.segs1) = 44, segments 2-45
#        for i in range(len(WS.segs1)):
#            ax.annotate('%s:%s'%(str(WS.segs1[i]),str(int(Ttimes[0][i]))), (WS.Pnts1[i,0], WS.Pnts1[i,1]), color=(r[i],g, b[i]), fontsize=8)
#            #ax.annotate('Seg %s, Time=%s'%(str(WS.segs1[i]),str(int(self.Ttime[i+1]))), (WS.Pnts1[i,0], WS.Pnts1[i,1]), color=(r[i],g, b[i]), fontsize=18)  ## for zoom plot only
#        for i in range(len(WS.segs2)):
#            ax.annotate('%s:%s'%(str(WS.segs2[i]),str(int(Ttimes[1][::-1][i]))), (WS.Pnts2[i,0], WS.Pnts2[i,1]), color=(r[i],g, b[i]), fontsize=8)
#        for i in range(len(WS.segs3)):
#            ax.annotate('%s:%s'%(str(WS.segs3[i]),str(int(Ttimes[2][::-1][i]))), (WS.Pnts3[i,0], WS.Pnts3[i,1]), color=(r[i],g, b[i]), fontsize=8)
#        for i in range(len(WS.segs4)):
#            ax.annotate('%s:%s'%(str(WS.segs4[i]),str(int(Ttimes[3][::-1][i]))), (WS.Pnts4[i,0], WS.Pnts4[i,1]), color=(r[i],g, b[i]), fontsize=8)
#        for i in range(len(WS.segs5)):
#            ax.annotate('%s:%s'%(str(WS.segs5[i]),str(int(Ttimes[4][::-1][i]))), (WS.Pnts5[i,0], WS.Pnts5[i,1]), color=(r[i],g, b[i]), fontsize=8)
            
        #### label segment without the segment ID
        for i in range(len(WS.segs1)):
            if Ttimes[0][i] != 0:
                ax.annotate('%s'%str(int(Ttimes[0][i])), (WS.Pnts1[i,0], WS.Pnts1[i,1]), color=colorlist[i], fontsize=8)
                #ax.annotate('Seg %s, Time=%s'%(str(WS.segs1[i]),str(int(self.Ttime[i+1]))), (WS.Pnts1[i,0], WS.Pnts1[i,1]), color=(r[i],g, b[i]), fontsize=18)  ## for zoom plot only
        for i in range(len(WS.segs2)):
            if Ttimes[1][::-1][i] != 0:
                ax.annotate('%s'%str(int(Ttimes[1][::-1][i])), (WS.Pnts2[i,0], WS.Pnts2[i,1]), color=colorlist[sum(nt[:2])-i-1], fontsize=8)
        for i in range(len(WS.segs3)):
            if Ttimes[2][::-1][i] != 0:
                ax.annotate('%s'%str(int(Ttimes[2][::-1][i])), (WS.Pnts3[i,0], WS.Pnts3[i,1]), color=colorlist[sum(nt[:3])-i-1], fontsize=8)
        for i in range(len(WS.segs4)):
            if Ttimes[3][::-1][i] != 0:
                ax.annotate('%s'%str(int(Ttimes[3][::-1][i])), (WS.Pnts4[i,0], WS.Pnts4[i,1]), color=colorlist[sum(nt[:4])-i-1], fontsize=8)
        for i in range(len(WS.segs5)):
            if Ttimes[4][::-1][i] != 0:
                ax.annotate('%s'%str(int(Ttimes[4][::-1][i])), (WS.Pnts5[i,0], WS.Pnts5[i,1]), color=colorlist[sum(nt[:])-i-1], fontsize=8)
            
        #pdb.set_trace()
            
        ax.set_xlabel('Easting [m]')
        ax.set_ylabel('Northing [m]')
        fig.colorbar(CS3)
        ax.set_aspect(True)
        fig.tight_layout()
        plt.show()
        
        
    def Ttime_plot_contourf(self, Ttimes, plotBranch=True, write2shp=False):
        """
        function Ttime_plot_contourf is different from function Ttime_plot 
        This function creates a collection of lines like contourf lines, 
        so that we can import the contourf lines to ArcGIS for visualization
        """
        
        #### call segment class for plotting
        Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
        WS = W2_Segmentation(Bthfile)
        WS.VisSeg2()
        
        nt = []
        for l in Ttimes:
            nt.append(len(l))
    
        
        
            
    def Xsearch(self, xin, zin, varin):
        """
        A searching function, 
        Input X and variable
        Output X with nonzero value (Note the masked value -99)
        Note: do not use the non-zero value at the surface layer
        """
        
        #if xin.shape != varin.shape:
        #    raise ValueError("X and Tracer shape not the same, check!")
        
        varin = np.ma.masked_array(varin,mask=varin==self.mask_value)
        
        ####quality control the masked array
        varin[(varin.mask==False)&(varin<1e-15)]=0
        
        ind = varin.nonzero()[0]  ## This nonzero assumption to find tracer time is not practical
        
        ## 50 percentile arrival time
        pp = 0.25  #### here I use 75 percentile, this value is subject to changes  0.25 works
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
#        if branchID == 1:
#            segID = np.searchsorted(xgrid, xxin) ## segment number where nonzero value first occurs
#        elif branchID == 5:
#            raise ValueError("In development.")  ## need to add the segmentID where branch 5 starts
        
        #pdb.set_trace()
        return np.unique(segID)
    
    
    def writeShpLines(self, WS, Ttimes):
        """
        write segment lines into ArcGIS readable shapefile
        Input segment object
        """
        import shapefile as shp
        import utm
        
        w = shp.Writer('traveltime')
        
        w.field('branch_ID','C')
        w.field('segment_ID','C')
        w.field('Lon_west','C')  #float - needed for coordinates
        w.field('Lat_west','C') 
        w.field('Lon_east','C')
        w.field('Lat_east','C')
        w.field('Travel_time','C')
        
        for i in range(len(WS.westPnts1)):
            westlat, westlon = utm.to_latlon(WS.westPnts1[i,0], WS.westPnts1[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(WS.eastPnts1[i,0], WS.eastPnts1[i,1], 14, 'U')
            #lines.append([[westlat, westlon], [eastlat, eastlon]])
            w.line([[[westlon, westlat], [eastlon, eastlat]]])
            w.record('branch1', str(WS.segs1[i]), "{:.3f}".format(westlon), "{:.3f}".format(westlat), \
                     "{:.3f}".format(eastlon), "{:.3f}".format(eastlat), "{:.1f}".format(Ttimes[0][i]))
            
        for i in range(len(WS.westPnts2)):
            westlat, westlon = utm.to_latlon(WS.westPnts2[i,0], WS.westPnts2[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(WS.eastPnts2[i,0], WS.eastPnts2[i,1], 14, 'U')
            w.line([[[westlon, westlat], [eastlon, eastlat]]])
            w.record('branch2', str(WS.segs2[::-1][i]), "{:.3f}".format(westlon), "{:.3f}".format(westlat), \
                     "{:.3f}".format(eastlon), "{:.3f}".format(eastlat), "{:.1f}".format(Ttimes[1][i]))
            
        for i in range(len(WS.westPnts3)):
            westlat, westlon = utm.to_latlon(WS.westPnts3[i,0], WS.westPnts3[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(WS.eastPnts3[i,0], WS.eastPnts3[i,1], 14, 'U')
            w.line([[[westlon, westlat], [eastlon, eastlat]]])
            w.record('branch3', str(WS.segs3[::-1][i]), "{:.3f}".format(westlon), "{:.3f}".format(westlat), \
                     "{:.3f}".format(eastlon), "{:.3f}".format(eastlat), "{:.1f}".format(Ttimes[2][i]))
        
        for i in range(len(WS.westPnts4)):
            westlat, westlon = utm.to_latlon(WS.westPnts4[i,0], WS.westPnts4[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(WS.eastPnts4[i,0], WS.eastPnts4[i,1], 14, 'U')
            w.line([[[westlon, westlat], [eastlon, eastlat]]])
            w.record('branch4', str(WS.segs4[::-1][i]), "{:.3f}".format(westlon), "{:.3f}".format(westlat), \
                     "{:.3f}".format(eastlon), "{:.3f}".format(eastlat), "{:.1f}".format(Ttimes[3][i]))
            
        for i in range(len(WS.westPnts5)):
            westlat, westlon = utm.to_latlon(WS.westPnts5[i,0], WS.westPnts5[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(WS.eastPnts5[i,0], WS.eastPnts5[i,1], 14, 'U')
            w.line([[[westlon, westlat], [eastlon, eastlat]]])
            w.record('branch5', str(WS.segs5[::-1][i]), "{:.3f}".format(westlon), "{:.3f}".format(westlat), \
                     "{:.3f}".format(eastlon), "{:.3f}".format(eastlat), "{:.1f}".format(Ttimes[4][i]))
        
        w.close()
        
        prj = open("traveltime.prj", "w") 
        epsg = 'GEOGCS["WGS 84",'
        epsg += 'DATUM["WGS_1984",'
        epsg += 'SPHEROID["WGS 84",6378137,298.257223563]]'
        epsg += ',PRIMEM["Greenwich",0],'
        epsg += 'UNIT["degree",0.0174532925199433]]'
        prj.write(epsg)
        prj.close()
        
        pdb.set_trace()
        
        
if __name__ == "__main__":
    
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191127_1115_tracer_test'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191127_1336_tracer_test'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191127_1631_tracer_test'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191202_1100_tracer_test'
    wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191213_1533_tracer_test'
    TTT = Tracer_Travel_Time(wdir)
    #TTT.travel_time()
    #TTT.Ttime_plot()
    TTT.travel_time_full_branch()
    