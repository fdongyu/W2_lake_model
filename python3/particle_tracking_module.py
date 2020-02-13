# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:15:39 2020

@author: dfeng
"""

"""
Calculate particle tracking with surface and bottom velocities
"""

import numpy as np
import matplotlib.pyplot as plt
import random

from w2_contour import W2_Contour
from bathymetry import W2_Bathymetry
from segmentation import W2_Segmentation

import pdb

class Particle_Tracking_Module(W2_Contour):
    
    mask_value = -99
    
    #branch2:
    DHS2 = 33
    #branch3:
    DHS3 = 34
    #bracnh4
    DHS4 = 42
    #branch5
    DHS5 = 43
    
    Dx = 1e-5  ##  longitudinal dispersion coefficient
    
    def __init__(self, workdir, **kwargs):
        
        self.workdir = workdir
        
        #self.Nt = 300  ## number of days to calculate particle transport
      
        self.read_output()
    
    
    def read_output(self):
        """
        read full model output
        """
        #### read full output
        #### self.X_flow, self.Z_flow, self.U, self.W
        self.Readcpl()
        
        
    def particle_tracking_model_1D(self, Np, Nt, InitialSeg, branchID, dt=1, transportSurface=True, transportBottom=True, travelTime=True):
        """
        particle tracking with velocity 
        Np -- number of particles
        Nt -- particle tracking period (unit: day)
        InitialSeg -- initial spill release segment ID
        dt -- time interval, default 1 day, unit: day
        """
        
        dt *= 24*3600.   #### conversion from day to seconds
        
        
        #### read surface and bottom velocities
        if branchID == 1:
            self.X_surface, self.Z_surface, self.U_surface, \
            self.X_bottom, self.Z_bottom, self.U_bottom = self.read_velocity(Nt, branchID=1)
            
        elif branchID == 5:
            X_surface1, Z_surface1, U_surface1, \
            X_bottom1, Z_bottom1, U_bottom1 = self.read_velocity(Nt, branchID=1)
            X_surface5, Z_surface5, U_surface5, \
            X_bottom5, Z_bottom5, U_bottom5 = self.read_velocity(Nt, branchID=5)
            
            #### read bathymetry information
            Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
            WB = W2_Bathymetry(Bthfile)
            pat = WB.VisBranch2(branchID=1)
            #### adding branch 5 to main branch 
            self.X_surface = []
            self.Z_surface = []
            self.U_surface = []
            
            self.X_bottom = []
            self.Z_bottom = []
            self.U_bottom = []
            for t in range(Nt):
                
                ## surface
                xind_surface = self.findNearest(WB.X[self.DHS5-1], X_surface1[t][:])
                xtem_surface_branch1 = np.asarray(X_surface1[t][xind_surface:]) - X_surface1[t][xind_surface-1] \
                                + X_surface5[t][-1]
                self.X_surface.append( X_surface5[t] + xtem_surface_branch1.tolist() )
                self.Z_surface.append( Z_surface5[t] + Z_surface1[t][xind_surface:] )
                self.U_surface.append( U_surface5[t] + U_surface1[t][xind_surface:] )
                
                ## bottom
                xind_bottom = self.findNearest(WB.X[self.DHS5-1], X_bottom1[t][:])
                xtem_bottom_branch1 = np.asarray(X_bottom1[t][xind_bottom:]) - X_bottom1[t][xind_bottom-1] \
                                + X_bottom5[t][-1]
                self.X_bottom.append( X_bottom5[t] + xtem_bottom_branch1.tolist() )
                self.Z_bottom.append( Z_bottom5[t] + Z_bottom1[t][xind_bottom:] )
                self.U_bottom.append( U_bottom5[t] + U_bottom1[t][xind_bottom:] )            
            
            
        #pdb.set_trace()
        
        #### read bathymetry information
        Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
        WB = W2_Bathymetry(Bthfile)
        pat = WB.VisBranch2(branchID)
        
        
        #### particle tracking calculation
        if transportSurface:
            
            #### particle location array
            self.location_x_surface = np.zeros([Np, Nt])   ####[Number of particles, time period]
            self.grid_x_surface = np.zeros([Nt])   #### surface water level at each x grid
            
            #### initial particle location        
            self.location_x_surface[:,0] = WB.X[InitialSeg-1]
            
            #### first order Euler algorithm: x(t+1) = x(t) + U*dt + R*sqrt(6 * Dx *dt) 
            for i in range(Np):
                for t in range(Nt-1):
                    xtem = np.abs(self.X_surface[t] - self.location_x_surface[i, t])
                    #### check if 
                    if xtem.min() < 1000:
                        #### query index
                        ind = np.argwhere(xtem==xtem.min())[0][0]
                        utem = self.U_surface[t][ind]
                        R = random.uniform(0,2) - 1     ## random number between [-1,1]
                        self.location_x_surface[i,t+1] = self.location_x_surface[i, t] + utem *dt + R*np.sqrt(6*self.Dx*dt)
                    elif xtem.min() > 1000:   ## there is no close grid point, water dries at this location
                        utem = 0
                        self.location_x_surface[i,t+1] = self.location_x_surface[i, t] + utem *dt
                    #if t in range(236, 238):
                    ## at these steps, water at the first several cells dries, X_surface starts at 9659, while location_x_surface is 8440. 
                    ## so particles do not move at these time steps 
                    
            
            for t in range(Nt):
                self.grid_x_surface[t] = self.Z_surface[t][0]
        
        
        if transportBottom:
            
            #### particle location array
            self.location_x_bottom = np.zeros([Np, Nt])
            self.grid_x_bottom = np.zeros([Nt])    #### bottom water level at each x grid
            
            #### initial particle location
            self.location_x_bottom[:,0] = WB.X[InitialSeg-1]
        
            #### first order Euler algorithm
            for i in range(Np):
                for t in range(Nt-1):
                    xtem = np.abs(self.X_bottom[t] - self.location_x_bottom[i, t])
                    #### check if 
                    if xtem.min() < 1000:
                        #### query index
                        ind = np.argwhere(xtem==xtem.min())[0][0]
                        utem = self.U_bottom[t][ind]
                        R = random.uniform(0,2) - 1     ## random number between [-1,1]
                        self.location_x_bottom[i,t+1] = self.location_x_bottom[i, t] + utem *dt + R*np.sqrt(6*self.Dx*dt)
                    elif xtem.min() > 1000:   ## there is no close grid point, water dries at this location
                        utem = 0
                        self.location_x_bottom[i,t+1] = self.location_x_bottom[i, t] + utem *dt
            
            for t in range(Nt):
                self.grid_x_bottom[t] = self.Z_bottom[t][0]
                
        #pdb.set_trace()
#        #### For testing only: visualize particle locations
#        iy = 0
#        plt.rcParams.update({'font.size': 16})
#        fig = plt.figure(figsize=(14,10))
#        ax = fig.add_subplot(211)
#        for i in range(Np):
#            ax.plot(self.location_x_surface[i], self.grid_x_surface+iy, 'o')
#            iy+=5
#            
#        ax2 = fig.add_subplot(212)
#        for i in range(Np):
#            ax2.plot(self.location_x_bottom[i], self.grid_x_bottom-iy, 'o')
#            iy-=5
#        plt.show()
        
        if travelTime:
            self.travel_time(Np, Nt, InitialSeg, branchID, self.location_x_bottom, write2shp=True)
            
        
        
    def travel_time(self, Np, Nt, InitialSeg, branchID, location_x, write2shp):
        """
        calculate travel time based on the particle tracking
        """
        
        if branchID == 1:
            
            #### read segment information
            Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
            WB = W2_Bathymetry(Bthfile)
            pat = WB.VisBranch2(branchID)
        
            #### create empty array for travel time
            Ttime = np.zeros([Np, WB.X.shape[0]])
            
            #### calculate travel time
            for i in range(Np):
                for tstep in range(Nt):
                    location_x_tem = location_x[i, tstep]
                        
                    ind = self.Xsearch(location_x_tem, WB.X)        
                    
                    if tstep == 0:
                        Ttime[i, InitialSeg-1:ind+1] = tstep + 1
                        ind_tem = ind
                    else:
                        ind_nonzero = np.nonzero(Ttime[i,:])[0].max()  ## only add travel time to zero elements
                        if ind > max(ind_tem, ind_nonzero):
                            Ttime[i, max(ind_tem+1, ind_nonzero+1):ind+1] = tstep + 1
                            print (Ttime[i,:])
                        ind_tem = ind
            #pdb.set_trace()
                    
            
            if write2shp:
                from myshapefile import writeShpLines_one_branch
                
                #### call segment class for plotting
                Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
                WS = W2_Segmentation(Bthfile)
                WS.VisSeg2()
                
                #### calculate the average among particles
                #### be careful about this average among particles 
                #Ttime = np.mean(Ttime, axis=0, dtype=np.int)
                #### simply average may yield smaller travel time at the downstream end
                #### because some particle may not travel to the last segment, which has zero travel time there
                #### so only average among non-zero values. 
                Ttime_avg = np.zeros([WB.X.shape[0]])
                for i in range(WB.X.shape[0]):
                    if i >= InitialSeg-1:
                        Ttime_avg[i] = Ttime[:,i][np.nonzero(Ttime[:,i])].mean()
                
                
                Ttime_avg = Ttime_avg[1:-1]
                Ttime_avg[Ttime_avg!=0] = Ttime_avg[-1] - Ttime_avg[Ttime_avg!=0]
                
                writeShpLines_one_branch(WS, Ttime_avg, shpname='particle_surface_traveltime_branch1')
                
        
        if branchID == 5:
            """
            Under development
            """
            #### read segment information for branch 5
            Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
            WB = W2_Bathymetry(Bthfile)
            pat = WB.VisBranch2(branchID)
            
            x_branch5 = WB.X   #### segment x coordinates for branch 5
            
            #### read segment information for branch 1
            Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
            WB = W2_Bathymetry(Bthfile)
            pat = WB.VisBranch2(branchID=1)
            
            x_branch1 = WB.X
            
            #### create empty array for travel time
            #### should not include the inactive cell at the end of branch5 when combining
            x_combined = x_branch5.tolist()[0:-1] + \
                        (x_branch1[self.DHS5-1:] - x_branch1[self.DHS5-1] + x_branch5[-2]).tolist()
            x_combined = np.asarray(x_combined)
            Ttime = np.zeros([Np, x_combined.shape[0]])
            
            
            #### calculate travel time
            for i in range(Np):
                for tstep in range(Nt):
                    location_x_tem = location_x[i, tstep]
                        
                    ind = self.Xsearch(location_x_tem, x_combined)        
                    
                    if tstep == 0:
                        Ttime[i, InitialSeg-1:ind+1] = tstep + 1
                        ind_tem = ind
                    else:
                        ind_nonzero = np.nonzero(Ttime[i,:])[0].max()  ## only add travel time to zero elements
                        if ind > max(ind_tem, ind_nonzero):
                            Ttime[i, max(ind_tem+1, ind_nonzero+1):ind+1] = tstep + 1
                            print (Ttime[i,:])
                        ind_tem = ind
            
            #### separate the travel time to branch 1 segments and branch 5 segments 
            Ttime5 = Ttime[:, 0: len(x_branch5[0:-1])+1]
            
            Ttime1 = np.zeros([Np, len(x_branch1)])
            Ttime1[:, self.DHS5:] = Ttime[:, len(x_branch5[0:-1])+1:]
            
            if write2shp:
                from myshapefile import writeShpLines
                
                #### call segment class for plotting
                Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
                WS = W2_Segmentation(Bthfile)
                WS.VisSeg2()
                
                #### calculate the average among particles
                Ttime_avg1 = np.zeros([Ttime1.shape[1]])
                for i in range(Ttime1.shape[1]):
                    if i >= self.DHS5 and len(Ttime1[:,i].nonzero()[0]) != 0:
                        Ttime_avg1[i] = Ttime1[:,i][np.nonzero(Ttime1[:,i])].mean()
                        
                Ttime_avg5 = np.zeros([Ttime5.shape[1]])
                for i in range(Ttime5.shape[1]):
                    if i >= InitialSeg-1:
                        Ttime_avg5[i] = Ttime5[:,i][np.nonzero(Ttime5[:,i])].mean()
                
                
                
                Ttime_avg1 = Ttime_avg1[1:-1]
                Ttime_avg5 = Ttime_avg5[1:-1]
                
                Ttimes_avg = [Ttime_avg1, Ttime_avg5]
                MaxTime = Ttimes_avg[0][-1]
                for Ttime in Ttimes_avg:
                    Ttime[Ttime!=0] = MaxTime - Ttime[Ttime!=0]
                
                
                writeShpLines(WS, Ttimes_avg, shpname='particle_bottom_traveltime_branch5')
                
                
                
            
            
            
    
    def Xsearch(self, xx, xgrid):
        """
        search the index of xx in xgrid
        note xx has to be greater than xgrid[i]
        """
        xtem = xx - xgrid
        
        xtem = xtem.tolist()
        
        m = min(i for i in xtem if i>=0)
        
        return xtem.index(m)
        
    
    def read_velocity(self, Nt, branchID=1):
        """
        read surface and bottom velocities, specific for each branch
        """
        
        #### read bathymetry information
        Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
        WB = W2_Bathymetry(Bthfile)
        pat = WB.VisBranch2(branchID)
        
        X_surface = []     ## [time period, x grid points]
        Z_surface = []
        U_surface = []
        
        X_bottom = []
        Z_bottom = []
        U_bottom = []
        
        
        for tstep in range(Nt):
            
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
        
            ## read velocity
            U = self.U[tstep][ind0:ind1+1]
            W = self.W[tstep][ind0:ind1+1]
            U = np.asarray(U)
            W = np.asarray(W)
            
            #### find surface and bottom velocity at tstep
            X_surface_tem, Z_surface_tem, U_surface_tem,    \
            X_bottom_tem,  Z_bottom_tem,  U_bottom_tem  =   \
                            self.surface_bottom_velocity(X_flow, Z_flow, U, W)
            
            X_surface.append(X_surface_tem)
            Z_surface.append(Z_surface_tem)
            U_surface.append(U_surface_tem)
            
            X_bottom.append(X_bottom_tem)
            Z_bottom.append(Z_bottom_tem)
            U_bottom.append(U_bottom_tem)
        
            #pdb.set_trace()
        return X_surface, Z_surface, U_surface, X_bottom, Z_bottom, U_bottom
    
        
        
    def surface_bottom_velocity(self, Xin, Zin, Uin, Win):
        """
        find surface and bottom velocity
        find the index of maximum or second maximum number
        https://stackoverflow.com/questions/6910641/how-do-i-get-indices-of-n-maximum-values-in-a-numpy-array?rq=1
        """
        
        #### find the index at each x grid
        #### Xin values at each x grid are the same
        ind_xgrid = [np.argwhere(Xin==xtem).flatten() for xtem in np.unique(Xin)]
        
        ## find surface and bottom velocity at each x grid
        X_surface = []
        Z_surface = []
        U_surface = []
        
        X_bottom = []
        Z_bottom = []
        U_bottom = []
        
        
        Nx = len(ind_xgrid) -1   ### Number of x grid points, note the last point has no velocity
        
        for i in range(Nx):
        #for i in [35]:
            
            ## index at each x grid point
            ind_x_tem = ind_xgrid[i]   
               
            #### find surface and bottom index
            mask = np.logical_or(Uin[ind_x_tem]!=self.mask_value, Win[ind_x_tem]!=self.mask_value) ## find invalid value  
            ind_surface = Zin[ind_x_tem][mask].argsort()[-2]    ## maximum (or second maximum) Z
            ind_bottom = Zin[ind_x_tem][mask].argsort()[2]    ## minimum (or second minimum) Z, note velocity at the bottom index
                                                              ## = 1 (second minimum) may have zero values, so make it 2 (third minimum)
            
            #### Note, there might some bugs for bottom velocity. 
            #### If velocity is zero, the particle will not move a lot, simply move with the dispersion ~O(10 m)
            #### But at branch 5,  U_bottom[t][12] is always 0, so particles get stagnant at ~10447 m, (X_bottom[t][12] = 10676 m)
            #### so double check: if at this segment, the bottom velocity is zero, if 0, take one layer up
            #### because at segment 12 (branch 5), the velocities at the three bottom layers are always zero
            icount = 2
            if i != Nx-1:
                while Uin[ind_x_tem][ind_bottom] == 0:
                    ind_bottom = Zin[ind_x_tem][mask].argsort()[icount] 
                    icount += 1
            
            
            ## surface
            X_surface_tem = Xin[ind_x_tem][ind_surface]
            Z_surface_tem = Zin[ind_x_tem][ind_surface]
            U_surface_tem = Uin[ind_x_tem][ind_surface]
            
            X_surface.append(X_surface_tem)
            Z_surface.append(Z_surface_tem)
            U_surface.append(U_surface_tem)
            
            ## bottom
            X_bottom_tem = Xin[ind_x_tem][ind_bottom]
            Z_bottom_tem = Zin[ind_x_tem][ind_bottom]
            U_bottom_tem = Uin[ind_x_tem][ind_bottom]
            
            X_bottom.append(X_bottom_tem)
            Z_bottom.append(Z_bottom_tem)
            U_bottom.append(U_bottom_tem)
            
        #pdb.set_trace()
            
        return X_surface, Z_surface, U_surface, X_bottom, Z_bottom, U_bottom
    
    
    def findNearest(self, xx, Xall):
        """
        find the index of the nearest element
        """
        
        xtem = np.abs(Xall - xx)
        
        return np.argwhere(xtem==xtem.min())[0][0]
            
    
        
        
        
if __name__ == "__main__": 
    
    wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191213_1533_tracer_test'
    
    PTM = Particle_Tracking_Module(wdir)
    #PTM.particle_tracking_model_1D(10, 250, 15, branchID=1, transportSurface=True, transportBottom=True, travelTime=True)
    PTM.particle_tracking_model_1D(10, 250, 15, branchID=5, transportSurface=True, transportBottom=True, travelTime=True)