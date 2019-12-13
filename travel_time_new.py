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
    
    
    def __init__(self, workdir, **kwargs):
        
        self.workdir = workdir
        
        
    def travel_time(self, endtime=1460, branchID=1):
        """
        calculate the travel time (or earliest arrival time) for conservative tracers
        timestep - when to calculate the travel time 
        branchID - the travel time of which branch
        """
        
        #### read bathymetry information
        Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
        WB = W2_Bathymetry(Bthfile)
        pat = WB.VisBranch2(branchID)
        
        
        #### create empty array for travel time
        self.Ttime = np.zeros_like(WB.X)
        
        ### read conservative tracer file
        self.Readcpl()
        
        #### calculate travel time ####
        #### loop over all time steps ####
        for tstep in range(endtime): 
            print 'Time step = %s \n'%str(tstep)
            if self.Ttime[-2] == 0:
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
                    ind1 = len(self.X_flow[tstep])
    
                ## tracer locations
                X_flow = TTT.X_flow[tstep][ind0:ind1+1]
                Z_flow = TTT.Z_flow[tstep][ind0:ind1+1]
                X_flow = np.asarray(X_flow)
                Z_flow = np.asarray(Z_flow)
            
                vartem = np.asarray( self.var_output['Tracer']['value'][tstep][ind0:ind1+1] )
                ## search nonzero X coordinates
                X_nonzero = self.Xsearch(X_flow, Z_flow, vartem)  
                ## search nonzero X index
                ind_nonzero = self.GridSearch(WB.X, X_nonzero, branchID)  
                ind_max = np.max(ind_nonzero)
                ## add the travel time to the pre-defined 1D array
                if self.Ttime[ind_max] == 0:
                    if len(self.Ttime.nonzero()[0]) != 0:
                        ind_tem = np.where(self.Ttime==self.Ttime.max())[0][-1]
                        self.Ttime[ind_tem+1: ind_max] = tstep
                    else:
                        self.Ttime[ind_max] = tstep
                print self.Ttime
            
            else:
                break    
        #pdb.set_trace()
        
        
    
    def Ttime_plot(self, plotBranch=True):
        """
        Visualize the travel time on W2 segments
        multiple lines' color changing the values corresponding to gradient color:
            https://stackoverflow.com/questions/55587674/how-to-make-multiple-lines-color-changing-with-values-corresponding-to-gradient
        Matplotlib - add colorbar to a sequence of line plots:
            https://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots
        """
        
        #### calculate travel time first
        self.travel_time(endtime=1460, branchID=1)
        
        
        #### call segment class for plotting
        Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
        WS = W2_Segmentation(Bthfile)
        WS.VisSeg2()
        
        
        #### line color
        self.Ttime[-2] = self.Ttime[-3] + 1
        zz = self.Ttime[1:-1]
        r = (zz.astype(np.float)-zz.min())/(zz.max()-zz.min())
        g = 0
        b = 1 - r
        
        #### colorbar
        #pdb.set_trace()
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
            ax.plot(WS.Pnts1[:,0], WS.Pnts1[:,1], '.-k')
            ax.plot(WS.Pnts2[:,0], WS.Pnts2[:,1], '.-b')
            ax.plot(WS.Pnts3[:,0], WS.Pnts3[:,1], '.-b')
            ax.plot(WS.Pnts4[:,0], WS.Pnts4[:,1], '.-b')
            ax.plot(WS.Pnts5[:,0], WS.Pnts5[:,1], '.-b')
        
        
        for i in range(len(WS.westPnts1)):
            ax.plot([WS.westPnts1[i,0], WS.eastPnts1[i,0]], [WS.westPnts1[i,1], WS.eastPnts1[i,1]], '-', color=(r[i], g, b[i]))
#        for i in range(len(WS.westPnts2)):
#            ax.plot([WS.westPnts2[i,0], WS.eastPnts2[i,0]], [WS.westPnts2[i,1], WS.eastPnts2[i,1]], '-k')
#        for i in range(len(WS.westPnts3)):
#            ax.plot([WS.westPnts3[i,0], WS.eastPnts3[i,0]], [WS.westPnts3[i,1], WS.eastPnts3[i,1]], '-k')
#        for i in range(len(WS.westPnts4)):
#            ax.plot([WS.westPnts4[i,0], WS.eastPnts4[i,0]], [WS.westPnts4[i,1], WS.eastPnts4[i,1]], '-k')
#        for i in range(len(WS.westPnts5)):
#            ax.plot([WS.westPnts5[i,0], WS.eastPnts5[i,0]], [WS.westPnts5[i,1], WS.eastPnts5[i,1]], '-k')
        
        #### len(WS.segs1) = 44, segments 2-45
        for i in range(len(WS.segs1)):
            #ax.annotate('%s:%s'%(str(WS.segs1[i]),str(int(self.Ttime[i+1]))), (WS.Pnts1[i,0], WS.Pnts1[i,1]), color=(r[i],g, b[i]), fontsize=8)
            ax.annotate('Seg %s, Time=%s'%(str(WS.segs1[i]),str(int(self.Ttime[i+1]))), (WS.Pnts1[i,0], WS.Pnts1[i,1]), color=(r[i],g, b[i]), fontsize=18)
            
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
        
        if xin.shape != varin.shape:
            raise ValueError("X and Tracer shape not the same, check!")
        
        varin = np.ma.masked_array(varin,mask=varin==self.mask_value)
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
        
        if branchID == 1:
            segID = np.searchsorted(xgrid, xxin) ## segment number where nonzero value first occurs
        elif branchID == 5:
            raise ValueError("In development.")  ## need to add the segmentID where branch 5 starts
        
        #pdb.set_trace()
        return np.unique(segID)
        
        
        
        
if __name__ == "__main__":
    #pdb.set_trace()
    #wdir = r'C:\Users\dfeng\Downloads\v42\Tests\20191112_tracer'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191113_tracer_test'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191127_1115_tracer_test'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191127_1336_tracer_test'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191127_1631_tracer_test'
    wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191202_1100_tracer_test'
    TTT = Tracer_Travel_Time(wdir)
    #TTT.travel_time()
    TTT.Ttime_plot()
    