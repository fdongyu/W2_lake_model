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
        
        
    def particle_tracking_model_1D(self, Np, Nt, InitialSeg, branchID, dt=1, transportSurface='True', transportBottom='True'):
        """
        particle tracking with velocity 
        Np -- number of particles
        Nt -- particle tracking period (unit: day)
        InitialSeg -- initial spill release segment ID
        dt -- time interval, default 1 day, unit: day
        """
        
        dt *= 24*3600.
        
        if branchID == 1:
            X_surface, Z_surface, U_surface, \
            X_bottom, Z_bottom, U_bottom = self.read_velocity(Nt, branchID=1)
            
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
            X_surface = []
            Z_surface = []
            U_surface = []
            
            X_bottom = []
            Z_bottom = []
            U_bottom = []
            for t in range(Nt):
                
                ## surface
                xind_surface = self.findNearest(WB.X[self.DHS5-2], X_surface1[t][:])
                xtem_surface_branch1 = np.asarray(X_surface1[t][xind_surface:]) - X_surface1[t][xind_surface-1] \
                                + X_surface5[t][-1]
                X_surface.append( X_surface5[t] + xtem_surface_branch1.tolist() )
                Z_surface.append( Z_surface5[t] + Z_surface1[t][xind_surface:] )
                U_surface.append( U_surface5[t] + U_surface1[t][xind_surface:] )
                
                ## bottom
                xind_bottom = self.findNearest(WB.X[self.DHS5-2], X_bottom1[t][:])
                xtem_bottom_branch1 = np.asarray(X_bottom1[t][xind_bottom:]) - X_bottom1[t][xind_bottom-1] \
                                + X_bottom5[t][-1]
                X_bottom.append( X_bottom5[t] + xtem_bottom_branch1.tolist() )
                Z_bottom.append( Z_bottom5[t] + Z_bottom1[t][xind_bottom:] )
                U_bottom.append( U_bottom5[t] + U_bottom1[t][xind_bottom:] )            
            
            
        #pdb.set_trace()
        
        #### read bathymetry information
        Bthfile = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
        WB = W2_Bathymetry(Bthfile)
        pat = WB.VisBranch2(branchID)
        
        if transportSurface:
            
            #### particle location array
            location_x_surface = np.zeros([Np, Nt])   ####[Number of particles, time period]
            grid_x_surface = np.zeros([Nt])   #### surface water level at each x grid
            
            #### initial particle location        
            location_x_surface[:,0] = WB.X[InitialSeg-1]
            
            #### first order Euler algorithm: x(t+1) = x(t) + U*dt + R*sqrt(6 * Dx *dt) 
            for i in range(Np):
                for t in range(Nt-1):
                    xtem = np.abs(X_surface[t] - location_x_surface[i, t])
                    #### check if 
                    if xtem.min() < 1000:
                        #### query index
                        ind = np.argwhere(xtem==xtem.min())[0][0]
                        utem = U_surface[t][ind]
                        R = random.uniform(0,2) - 1     ## random number between [-1,1]
                        location_x_surface[i,t+1] = location_x_surface[i, t] + utem *dt + R*np.sqrt(6*self.Dx*dt)
                    elif xtem.min() > 1000:   ## there is no close grid point, water dries at this location
                        utem = 0
                        location_x_surface[i,t+1] = location_x_surface[i, t] + utem *dt
                    #if t in range(236, 238):
                    ## at these steps, water at the first several cells dries, X_surface starts at 9659, while location_x_surface is 8440. 
                    ## so particles do not move at these time steps 
                    #    pdb.set_trace()
            
            for t in range(Nt):
                grid_x_surface[t] = Z_surface[t][0]
        
        
        if transportBottom:
            
            #### particle location array
            location_x_bottom = np.zeros([Np, Nt])
            grid_x_bottom = np.zeros([Nt])    #### bottom water level at each x grid
            
            #### initial particle location
            location_x_bottom[:,0] = WB.X[InitialSeg-1]
        
            #### first order Euler algorithm
            for i in range(Np):
                for t in range(Nt-1):
                    xtem = np.abs(X_bottom[t] - location_x_bottom[i, t])
                    #### check if 
                    if xtem.min() < 1000:
                        #### query index
                        ind = np.argwhere(xtem==xtem.min())[0][0]
                        utem = U_bottom[t][ind]
                        R = random.uniform(0,2) - 1     ## random number between [-1,1]
                        location_x_bottom[i,t+1] = location_x_bottom[i, t] + utem *dt + R*np.sqrt(6*self.Dx*dt)
                    elif xtem.min() > 1000:   ## there is no close grid point, water dries at this location
                        utem = 0
                        location_x_bottom[i,t+1] = location_x_bottom[i, t] + utem *dt
            
            for t in range(Nt):
                grid_x_bottom[t] = Z_bottom[t][0]
                
        pdb.set_trace()
        
        
        #### visualize particle locations
        iy = 0
        plt.rcParams.update({'font.size': 16})
        fig = plt.figure(figsize=(14,10))
        ax = fig.add_subplot(211)
        for i in range(Np):
            ax.plot(location_x_surface[i], grid_x_surface+iy, 'o')
            iy+=5
            
        ax2 = fig.add_subplot(212)
        for i in range(Np):
            ax2.plot(location_x_bottom[i], grid_x_bottom-iy, 'o')
            iy-=5
        plt.show()
        
        
    
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
            icount = 3
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
    PTM.particle_tracking_model_1D(10, 250, 15, branchID=5, transportSurface='True', transportBottom='True')
    