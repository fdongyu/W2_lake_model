# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 15:55:18 2019

@author: dfeng
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

import pdb


## file list
## initialparticle.csv                        the initial state of all particles added or to be added to the system
## finalparticle.csv                          the final state of all particles released in the model
## Envrprf_v_particle.csv, envrprf_t_particle.csv, envrprf_depth_particle.csv     histograms of each particle released in the model domain for velocity, temperature and depth.

##  This file has X, Z, Water(=1 if in the water and -1 if above water or in bottom), U(HorizVel), W(VertVel), Temperature, 
##  DissolvedOxygen, KLayer, and ISegment for each time interval of output.  background to overlay over the particle file
## Branch1.dat, Branch2.dat, ...

## This file has the location of all the particles as a function of time
## Part1.dat, Part2.dat, ...    



class W2_Particle_Tracking(object):
    """
    general class for reading and visualizing CE-QUAL-W2 modeled particle tracking
    """
    mask_value = -99
    
    def __init__(self, workdir, **kwargs):
        self.__dict__.update(kwargs)
        
        self.workdir = workdir
        
    def ReadBackground(self, branchID='1'):
        """
        read Branch1.dat, ... 
        The background flow
        """
        filename = '%s\\Branch%s.dat'%(self.workdir, branchID)
        
        ## Step 1: read the timeframe info: JDAY
        f = open(filename, 'r')
        JDAYs = []
        for s in f:
            line = s.split()
            if 'ZONE' in line:
                #print 'Reading Model Run Time JDAY = %s ... \n'%line[2]
                #pdb.set_trace()
                JDAYs.append(line[2][:-2])
        f.close()
        
        ## Step 2: read the background flow at each time step
        self.backflow_sections = []
        recording = False
        f = open(filename, 'r')
        for i in range(len(JDAYs)):
            ttem = JDAYs[i]
            new_section = []
            for s in f:
                line = s.split()
                if recording is False:
                    if ttem+'",' in line:
                        #pdb.set_trace()
                        print 'Reading Model Run Time JDAY = %s ... \n'%ttem
                        recording = True
                        #pdb.set_trace()
                        new_section.append(line)
                elif recording is True:
                    new_section.append(line)
                    if len(line)==0:
                        recording = False
                        break
            
            if i == len(JDAYs)-1:
                self.backflow_sections.append(new_section[1:])
            else:
                self.backflow_sections.append(new_section[1:-1])
        
        ## Step 3: get time info     
        self.runtimes_flow = self.JDAY_con(JDAYs)
        
        ##Step 4: get X and Z coordinates    
        self.X_flow = []
        self.Z_flow = []
        self.U = []
        self.W = []
        self.T = []
        for i in range(len(self.backflow_sections)):
            section = self.backflow_sections[i]
            xtem = [float(l[0]) for l in section]
            ztem = [float(p[1]) for p in section]
            utem = [float(r[3]) for r in section]
            wtem = [float(q[4]) for q in section]
            ttem = [float(t[5]) for t in section]
            self.X_flow.append(xtem)
            self.Z_flow.append(ztem)
            self.U.append(utem)
            self.W.append(wtem)
            self.T.append(ttem)
        
        #pdb.set_trace()
            
        
        
        
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
            
            
            
        #pdb.set_trace()
        
    def VisParticles(self, timestep=-1, PlotFlow=True, PlotTemp=True, PlotGrid=True):
        """
        visualize the particle trajectories at a given time step
        set PlotTemp == False for now, need to figure out how to contour plot 1D array
        """
        self.ReadParticles('1')
        
        X = self.X[timestep]
        Z = self.Z[timestep]
        
        if PlotFlow == True or PlotTemp == True:
            self.ReadBackground('1')
            X_flow = self.X_flow[timestep]
            Z_flow = self.Z_flow[timestep]
        
        plt.rcParams.update({'font.size': 18})
        #fig = plt.figure(figsize=(11.5,10))
        fig = plt.figure(figsize=(11.5,8))
        ax = fig.add_subplot(111)
        ax.plot(X, Z, '.', color = 'r')
        if PlotFlow == True:
            U = self.U[timestep]
            W = self.W[timestep]
            
            ## mask data
            X_flow = np.asarray(X_flow)
            Z_flow = np.asarray(Z_flow)
            U = np.asarray(U)
            W = np.asarray(W)
            maskuv = np.logical_or(U != 0,W != 0)
            
            scale = 1.
            scale = 100./scale
            Q = ax.quiver(X_flow[maskuv], Z_flow[maskuv], U[maskuv]*100., W[maskuv]*100.,
                               zorder=5, width=0.001, headwidth=4, headlength=4.5,
                              scale=scale, color='b')
            qk = ax.quiverkey(Q, 0.15, 0.15, 1, r'$1 \frac{cm}{s}$', labelpos='W',fontproperties={'weight': 'bold','size':20})
        
        if PlotTemp == True:
            import matplotlib.tri as tri
            #from scipy.interpolate import griddata
            T = self.T[timestep]
            U = np.asarray(self.U[timestep])
            W = np.asarray(self.W[timestep])
            
            ## mask data
            X_flow = np.asarray(X_flow)
            Z_flow = np.asarray(Z_flow)
            T = np.asarray(T)
            
            triang = tri.Triangulation(X_flow, Z_flow)
            isbad = np.less_equal(np.asarray(self.T[-1]), 0) 
            #isbad = np.equal(U, 0) & np.equal(W, 0)
            mask = np.any(np.where(isbad[triang.triangles], True, False), axis=1)
            triang.set_mask(mask)
            
            
            #T[T==0] = self.mask_value
            #T = np.ma.masked_array(T,mask=T==self.mask_value)
            
            #T_limits = [0,35]
            T_limits = [T.min(),T.max()]
            levels = np.linspace(T_limits[0], T_limits[1], 100)
            cs = ax.tricontourf(triang, T, cmap=plt.cm.bone, levels=levels)
            #cs = ax.tricontourf(X_flow, Z_flow, T, cmap=plt.cm.bone, levels=levels)
            
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cb = fig.colorbar(cs, cax=cax, orientation='vertical')
            cb.ax.tick_params(labelsize=18)
            cb.ax.yaxis.offsetText.set_fontsize(14)
            cb.set_label('Temperature', fontsize=18)
            
        if PlotGrid == True:
            from bathymetry import W2_Bathymetry
            filename = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
            WB = W2_Bathymetry(filename)
            pat = WB.VisBranch2(1)
            for sq in pat:
                ax.add_patch(sq)
            ax.autoscale_view()
        
          
        timestr = datetime.strftime(self.runtimes[timestep],'%Y-%m-%d')
        ax.title.set_text('Time: %s'%timestr)
        #ax.set_xlim([4000, 12000])
        #ax.set_ylim([135, 150])
        ax.set_ylim([135, 160])
        ax.set_xlabel('Distance from upstream (m)')
        ax.set_ylabel('Water Depth (m)')
        #ax.yaxis.grid(True)
        #ax.xaxis.grid(True)
        #plt.show()
        plt.savefig('particle_tracks_%s.png'%str(timestep))
        #plt.savefig('example_%s.png'%str(timestep))
        plt.close()
        
    
        
    def VisParticles_full(self):
        """
        visualize the particle trajectories, the particles at all time steps are plotted
        """
        self.ReadParticles('1')
    
        X = [item for sublist in self.X for item in sublist]
        Z = [item for sublist in self.Z for item in sublist]
        
        plt.rcParams.update({'font.size': 18})
        fig = plt.figure(figsize=(11.5,6))
        ax = fig.add_subplot(111)
        ax.plot(X, Z, '.', color = 'k')
        ax.set_xlim([4000, 12000])
        ax.set_ylim([140, 150])
        ax.yaxis.grid(True)
        ax.xaxis.grid(True)
        plt.savefig('simple_particle_tracks.png')
        plt.close()
        #plt.show()
        
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
    wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\particle_tracking_test\20191115_1640_test1'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\particle_tracking_test\20191121_1112_test2'
    WPT = W2_Particle_Tracking(wdir)
    WPT.VisParticles(11, PlotGrid=True)
    
    