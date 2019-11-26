# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 08:18:03 2019

@author: dfeng
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from segmentation import W2_Segmentation

import pdb



class Particle_Travel_Time(W2_Segmentation):
    """
    general class for visualizing particle travel time
    """
    
    def __init__(self, workdir, **kwargs):
        self.__dict__.update(kwargs)
        
        self.workdir = workdir
        
        self.bathfile = '%s\\Bth_WB1.npt'%self.workdir
        
        self.readBathymetry()
        
    
    def TravelTime(self, timelist):
        """
        read a list of travel times
        """
        
        XX = []
        YY = []
        Ttravel = []
        
        
        for i in range(len(timelist)):
            
            ttem = timelist[i]
            filename = '%s\\finalparticle\\finalparticle_%s.csv'%(self.workdir, str(ttem))
            print "Reading file %s ... \n"%filename
            Part, Seg, Layer, Branch, Xlocation, LaterDist = self.ReadFinalcsv(filename)
            
            x, y = self.CalcXY(Part, Seg, Layer, Branch, Xlocation, LaterDist)
            tt = np.ones_like(x) * ttem   # travel time
            
            for j in range(len(x)):
                XX.append(x[j])
                YY.append(y[j])
                Ttravel.append(tt[j])
            
            
        XX = np.asarray(XX)
        YY = np.asarray(YY)
        Ttravel = np.asarray(Ttravel)
        
        outarray = np.vstack((Ttravel, XX, YY)).T
        np.savetxt('contour_data.txt', outarray)
        
        
#        ## save to txt file
#        f = open('contour_data.txt', 'w')
#        for i in range(len(XX)):
#            f.write("%s  %s  %s\n"%str(XX[i]), str(YY[i]), str(Ttravel[i]))
#        f.close()
        
        #pdb.set_trace()
        
    
    def ContourPlot(self, PlotGrid=True):
        """
        create contour plot for travel time
        """
        
        inarray = np.loadtxt('contour_data.txt')
        
        T = inarray[:,0]
        
        T = T - T.min()
        X = inarray[:,1]
        Y = inarray[:,2]
        
        plt.rcParams.update({'font.size': 18})
        fig = plt.figure(figsize=(11.5,8))
        ax = fig.add_subplot(111)
        
        if PlotGrid==True:
            
            self.VisSeg2()
            
            ax.plot(self.Pnts1[:,0], self.Pnts1[:,1], '.-k')
            ax.plot(self.Pnts2[:,0], self.Pnts2[:,1], '.-b')
            ax.plot(self.Pnts3[:,0], self.Pnts3[:,1], '.-b')
            ax.plot(self.Pnts4[:,0], self.Pnts4[:,1], '.-b')
            ax.plot(self.Pnts5[:,0], self.Pnts5[:,1], '.-b')
            
            for i in range(len(self.westPnts1)):
                ax.plot([self.westPnts1[i,0], self.eastPnts1[i,0]], [self.westPnts1[i,1], self.eastPnts1[i,1]], '-k')
            for i in range(len(self.westPnts2)):
                ax.plot([self.westPnts2[i,0], self.eastPnts2[i,0]], [self.westPnts2[i,1], self.eastPnts2[i,1]], '-k')
            for i in range(len(self.westPnts3)):
                ax.plot([self.westPnts3[i,0], self.eastPnts3[i,0]], [self.westPnts3[i,1], self.eastPnts3[i,1]], '-k')
            for i in range(len(self.westPnts4)):
                ax.plot([self.westPnts4[i,0], self.eastPnts4[i,0]], [self.westPnts4[i,1], self.eastPnts4[i,1]], '-k')
            for i in range(len(self.westPnts5)):
                ax.plot([self.westPnts5[i,0], self.eastPnts5[i,0]], [self.westPnts5[i,1], self.eastPnts5[i,1]], '-k')
            
            for i in range(len(self.segs1)):
                ax.annotate('%s'%str(self.segs1[i]), (self.Pnts1[i,0], self.Pnts1[i,1]), color='b', fontsize=8)
            for i in range(len(self.segs2)):
                ax.annotate('%s'%str(self.segs2[i]), (self.Pnts2[i,0], self.Pnts2[i,1]), color='b', fontsize=8)
            for i in range(len(self.segs3)):
                ax.annotate('%s'%str(self.segs3[i]), (self.Pnts3[i,0], self.Pnts3[i,1]), color='b', fontsize=8)
            for i in range(len(self.segs4)):
                ax.annotate('%s'%str(self.segs4[i]), (self.Pnts4[i,0], self.Pnts4[i,1]), color='b', fontsize=8)
            for i in range(len(self.segs5)):
                ax.annotate('%s'%str(self.segs5[i]), (self.Pnts5[i,0], self.Pnts5[i,1]), color='b', fontsize=8)
        ax.set_aspect(True)
        
        cmap = plt.set_cmap('bone_r')
        cs = ax.tricontourf(X, Y, T, cmap=cmap)
        
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(cs, cax=cax, orientation='vertical')
        cb.ax.tick_params(labelsize=18)
        cb.ax.yaxis.offsetText.set_fontsize(14)
        cb.set_label('Travel Time', fontsize=18)
        
        plt.show()
        
        #pdb.set_trace()
        
        
    
    def CalcXY(self, Part, Seg, Layer, Branch, Xlocation, LaterDist):
        """
        calculate x and y coordinates for each particle 
        """
        nn = len(Part)
        
        x = []
        y = []
        
        for i in range(nn):
            particleID = Part[i]
            segID = Seg[i]
            branchID = Branch[i]
            layerID = Layer[i]
            laterDist = LaterDist[i]
            print "Reading particle %s ... \n"%str(particleID)
            
            segs = self.lyrs_segID[self.lyrs_branchID==branchID]
            seg_length = np.asarray(self.seg_length)[self.lyrs_branchID==branchID]
            seg_ori = np.asarray(self.seg_ori)[self.lyrs_branchID==branchID]
            lyrs = self.lyrs[self.lyrs_branchID==branchID]
            ## truncate boundary cells 
            segs = segs[1:-1]
            seg_length = seg_length[1:-1]
            seg_ori = seg_ori[1:-1]
            lyrs = lyrs[1:-1]
            
            ## read corresponding seg points info
            segs, Pnts= self.BranchPnt(branchID=branchID)
            
            ## goes longitudinal direction
            xtem = Pnts[segID-2][0] - Xlocation[i]*np.sin(seg_ori[segID-2])
            ytem = Pnts[segID-2][1] - Xlocation[i]*np.cos(seg_ori[segID-2])
            
            ## goes Lateral direction starting from left bank
            x.append( xtem - lyrs[segID-2][layerID]/2.*np.cos(np.pi*2-seg_ori[segID-2]) +  laterDist*np.cos(np.pi*2-seg_ori[segID-2]) )
            y.append( ytem - lyrs[segID-2][layerID]/2.*np.sin(np.pi*2-seg_ori[segID-2]) +  laterDist*np.sin(np.pi*2-seg_ori[segID-2]) )
        
        
        return x, y
        
        
    
    def ReadFinalcsv(self, filename):
        """
        Read the csv file, final particles location
        """
        
        column_names = ['Part', 'Seg', 'Xlocation', 'Layer', 'VerticalDistfromTop', \
                        'LateralDistfromLeftBank', 'Branch', 'ParticleInModel', \
                        'JDAYleftsystem', 'DetentionTime', 'RemovalMechanism', \
                        'SedVelocity', 'DateStart']
        
        #df = pd.read_fwf(filename, skiprows=1,names = column_names)
        df = pd.read_csv(filename, skiprows=1,names = column_names, delimiter=',')
        
        Part = [int(p) for p in df['Part'].values]
        Seg = [int(s) for s in df['Seg'].values]
        Layer = [int(y) for y in df['Layer'].values]
        Branch = [int(b) for b in df['Branch'].values]
        Xlocation = [float(x) for x in df['Xlocation'].values]
        LaterDist = [float(l) for l in df['LateralDistfromLeftBank'].values]
        
        return Part, Seg, Layer, Branch, Xlocation, LaterDist
        
        
        

#### For testing ####       
        
if __name__ == "__main__":
    
    wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\particle_tracking_test\20191125_1229_test5'
    PTT = Particle_Travel_Time(wdir)
    
    timelist = [12,13,14,15,16,17,18,19,20,21,22,23,24]
    #timelist = [12,13]
    
    #PTT.TravelTime(timelist)
    PTT.ContourPlot()