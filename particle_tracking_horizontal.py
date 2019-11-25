# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:29:47 2019

@author: dfeng
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pandas as pd
import utm

from segmentation import W2_Segmentation

import pdb


class W2_Particle_Tracking(object):
    """
    general class for reading and visualizing CE-QUAL-W2 modeled particle tracking
    """
    
    
    def __init__(self, workdir, **kwargs):
        self.__dict__.update(kwargs)
        
        self.workdir = workdir
        
    
    def ReadParticles(self, branchID='1'):
        """
        read Part1.dat, ...
        """
        filename = '%s\\Part%s.dat'%(self.workdir, branchID)
        
        
        ## Step 1: read the timeframe info: JDAY
        f = open(filename, 'r')
        JDAYs = []
        ## Step 1: Read the number of time steps
        for s in f:
            line = s.split()
            if 'ZONE' in line:
                #print 'Reading Model Run Time JDAY = %s ... \n'%line[2]
                #pdb.set_trace()
                JDAYs.append(line[2][:-2])
        f.close()
                
        ## Step 2: read the particle locations at each time step
        self.output_sections = []
        recording = False
        f = open(filename, 'r')
        for i in range(len(JDAYs)):
            ttem = JDAYs[i]
            new_section = []
            for s in f:
                line = s.split()
                if recording is False:
                    if ttem+'",' in line:
                        print 'Reading Model Run Time JDAY = %s ... \n'%ttem
                        recording = True
                        #pdb.set_trace()
                        new_section.append(line)
                elif recording is True:
                    new_section.append(line)
                    if 'TEXT' in line:
                        recording = False
                        break
            
            self.output_sections.append(new_section[1:-1])
            #pdb.set_trace()
            
        
        if len(JDAYs) != len(self.output_sections):
            raise Exception('Read data not successful, double check!!\n')
            
            
        ## Step 3: get time info     
        self.runtimes = self.JDAY_con(JDAYs)
        
        ## Step 4: get X and Z coordinates    
        self.X = []
        self.Z = []
        for i in range(len(self.output_sections)):
            section = self.output_sections[i]
            xtem = [float(l[0]) for l in section]
            ztem = [float(p[1]) for p in section]
            self.X.append(xtem)
            self.Z.append(ztem)
            #pdb.set_trace()
            
    def ReadFinalcsv(self):
        """
        Read the csv file, final particles location
        """
        filename = '%s\\finalparticle.csv'%(self.workdir)
        
        column_names = ['Part', 'Seg', 'Xlocation', 'Layer', 'VerticalDistfromTop', \
                        'LateralDistfromLeftBank', 'Branch', 'ParticleInModel', \
                        'JDAYleftsystem', 'DetentionTime', 'RemovalMechanism', \
                        'SedVelocity', 'DateStart']
        
        #df = pd.read_fwf(filename, skiprows=1,names = column_names)
        df = pd.read_csv(filename, skiprows=1,names = column_names, delimiter=',')
        
        self.Part = [int(p) for p in df['Part'].values]
        self.Seg = [int(s) for s in df['Seg'].values]
        self.Layer = [int(y) for y in df['Layer'].values]
        self.Branch = [int(b) for b in df['Branch'].values]
        self.Xlocation = [float(x) for x in df['Xlocation'].values]
        self.LaterDist = [float(l) for l in df['LateralDistfromLeftBank'].values]
        
        
        
        
        
        #pdb.set_trace()
    
    
            
    def VisFinal(self, PlotGrid=True, saveshp=True):
        """
        visualize the particle trajectories at a given time step
        set PlotTemp == False for now, need to figure out how to contour plot 1D array
        """
        
        self.ReadFinalcsv()
        
        ## call the segmentation class to get some information
        filename = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\particle_tracking_test\20191111_baseline4\Bth_WB1.npt'
        WB = W2_Segmentation(filename)
        WB.readBathymetry()
        
        
        x = []
        y = []
        for i in range(len(self.Part)):
        #for i in [0,1,2,3,4,5]:
            ## for a specific particle
            particleID = self.Part[i]
            segID = self.Seg[i]
            branchID = self.Branch[i]
            layerID = self.Layer[i]
            laterDist = self.LaterDist[i]
            print "Reading particle %s ... \n"%str(particleID)
            
            ## search the corresponding seg info
            segs = WB.lyrs_segID[WB.lyrs_branchID==branchID]
            seg_length = np.asarray(WB.seg_length)[WB.lyrs_branchID==branchID]
            seg_ori = np.asarray(WB.seg_ori)[WB.lyrs_branchID==branchID]
            lyrs = WB.lyrs[WB.lyrs_branchID==branchID]
            ## truncate boundary cells 
            segs = segs[1:-1]
            seg_length = seg_length[1:-1]
            seg_ori = seg_ori[1:-1]
            lyrs = lyrs[1:-1]
            
            ## read corresponding seg points info
            segs, Pnts= WB.BranchPnt(branchID=branchID)
            
            ## goes longitudinal direction
            xtem = Pnts[segID-2][0] - self.Xlocation[i]*np.sin(seg_ori[segID-2])
            ytem = Pnts[segID-2][1] - self.Xlocation[i]*np.cos(seg_ori[segID-2])
            
            ## goes Lateral direction starting from left bank
            #x.append( xtem +  lyrs[segID-2][layerID]/2.*np.cos(np.pi*2-seg_ori[segID-2]) )
            #y.append( ytem +  lyrs[segID-2][layerID]/2.*np.sin(np.pi*2-seg_ori[segID-2]) )
            x.append( xtem - lyrs[segID-2][layerID]/2.*np.cos(np.pi*2-seg_ori[segID-2]) +  laterDist*np.cos(np.pi*2-seg_ori[segID-2]) )
            y.append( ytem - lyrs[segID-2][layerID]/2.*np.sin(np.pi*2-seg_ori[segID-2]) +  laterDist*np.sin(np.pi*2-seg_ori[segID-2]) )
        
            #pdb.set_trace()
        
        if saveshp == True:
            self.writeShp(x, y)
        
        
        plt.rcParams.update({'font.size': 18})
        #fig = plt.figure(figsize=(11.5,10))
        fig = plt.figure(figsize=(11.5,8))
        ax = fig.add_subplot(111)
        ax.plot(x, y, '.', color = 'r', markersize=8)
        
        if PlotGrid==True:
            
            WB.VisSeg2()
            
            ax.plot(WB.Pnts1[:,0], WB.Pnts1[:,1], '.-k')
            ax.plot(WB.Pnts2[:,0], WB.Pnts2[:,1], '.-b')
            ax.plot(WB.Pnts3[:,0], WB.Pnts3[:,1], '.-b')
            ax.plot(WB.Pnts4[:,0], WB.Pnts4[:,1], '.-b')
            ax.plot(WB.Pnts5[:,0], WB.Pnts5[:,1], '.-b')
            
            for i in range(len(WB.westPnts1)):
                ax.plot([WB.westPnts1[i,0], WB.eastPnts1[i,0]], [WB.westPnts1[i,1], WB.eastPnts1[i,1]], '-k')
            for i in range(len(WB.westPnts2)):
                ax.plot([WB.westPnts2[i,0], WB.eastPnts2[i,0]], [WB.westPnts2[i,1], WB.eastPnts2[i,1]], '-k')
            for i in range(len(WB.westPnts3)):
                ax.plot([WB.westPnts3[i,0], WB.eastPnts3[i,0]], [WB.westPnts3[i,1], WB.eastPnts3[i,1]], '-k')
            for i in range(len(WB.westPnts4)):
                ax.plot([WB.westPnts4[i,0], WB.eastPnts4[i,0]], [WB.westPnts4[i,1], WB.eastPnts4[i,1]], '-k')
            for i in range(len(WB.westPnts5)):
                ax.plot([WB.westPnts5[i,0], WB.eastPnts5[i,0]], [WB.westPnts5[i,1], WB.eastPnts5[i,1]], '-k')
            
            for i in range(len(WB.segs1)):
                ax.annotate('%s'%str(WB.segs1[i]), (WB.Pnts1[i,0], WB.Pnts1[i,1]), color='b', fontsize=8)
            for i in range(len(WB.segs2)):
                ax.annotate('%s'%str(WB.segs2[i]), (WB.Pnts2[i,0], WB.Pnts2[i,1]), color='b', fontsize=8)
            for i in range(len(WB.segs3)):
                ax.annotate('%s'%str(WB.segs3[i]), (WB.Pnts3[i,0], WB.Pnts3[i,1]), color='b', fontsize=8)
            for i in range(len(WB.segs4)):
                ax.annotate('%s'%str(WB.segs4[i]), (WB.Pnts4[i,0], WB.Pnts4[i,1]), color='b', fontsize=8)
            for i in range(len(WB.segs5)):
                ax.annotate('%s'%str(WB.segs5[i]), (WB.Pnts5[i,0], WB.Pnts5[i,1]), color='b', fontsize=8)
        
        #ax.set_xlim([-20000,10000])
        ax.set_xlabel('Easting [m]')
        ax.set_ylabel('Northing [m]')
        ax.set_aspect(True)
        #plt.show()
        plt.savefig('particle_location_horizontal.png')
        plt.close()
        #pdb.set_trace()
        
    
    def writeShp(self, x, y):
        """
        write particle coordinates to GIS shapefile
        https://gis.stackexchange.com/questions/119160/using-pyshp-to-create-polygon-shapefiles
        https://salsa.debian.org/debian-gis-team/pyshp
        """
        import shapefile as shp
        
        ## step 1: convert UTM to latitude and longitude
        lat = np.zeros_like(x)
        lon = np.zeros_like(y)        
        for i in range(len(x)):
            lat[i], lon[i] = utm.to_latlon(x[i], y[i], 14, 'U')
    
        
        ## step 2: write to shapefile
        w = shp.Writer('final_particle')
        w.field('Lon','C')
        w.field('Lat','C') #float - needed for coordinates
        w.field('PID', 'C')
        
        for i in range(len(x)):
            w.point(lon[i], lat[i])
            w.record("{:.3f}".format(lon[i]), "{:.3f}".format(lat[i]), str(i))
            #w.record(str(i))
    
        w.close()
        
#        prj = open("final_particle.prj", "w") 
#        epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]' 
#        prj.write(epsg) 
#        prj.close()
        
        prj = open("final_particle.prj", "w") 
        epsg = 'GEOGCS["WGS 84",'
        epsg += 'DATUM["WGS_1984",'
        epsg += 'SPHEROID["WGS 84",6378137,298.257223563]]'
        epsg += ',PRIMEM["Greenwich",0],'
        epsg += 'UNIT["degree",0.0174532925199433]]'
        prj.write(epsg)
        prj.close()
        
        #pdb.set_trace()
    
    
    def JDAY_con(self,JDAY):
        """
        funtion adapted from Justin's code
        """
        basedate = datetime(2010,12,31)
        return_date = [basedate+timedelta(days = float(ii)) for ii in JDAY]
        return return_date    
        
        
        
#### For testing ####       
        
if __name__ == "__main__":
    
    #wdir = r'C:\Users\dfeng\Downloads\v42\Examples\Particle tracking_DF'
    #wdir = r'C:\Users\dfeng\Downloads\v42\Tests\20191106_baseline'
    #wdir = r'C:\Users\dfeng\Downloads\v42\Tests\20191108_baseline2'
    #wdir = r'C:\Users\dfeng\Downloads\v42\Tests\20191111_baseline3'
    #wdir = r'C:\Users\dfeng\Downloads\v42\Tests\20191111_baseline4'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\particle_tracking_test\20191115_1404_test0'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\particle_tracking_test\20191115_1640_test1'
    wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\particle_tracking_test\20191121_1112_test2'
    WPT = W2_Particle_Tracking(wdir)
    WPT.VisFinal(PlotGrid=True)