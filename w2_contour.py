# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 16:15:42 2019

@author: dfeng
"""
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

import pdb


class W2_Contour(object):
    """
    general class for reading and visualizing CE-QUAL-W2 modeled concentration
    """
    
    mask_value = -99
    
    def __init__(self, workdir, **kwargs):
        self.__dict__.update(kwargs)
        
        self.workdir = workdir
        
    def Readcpl(self):
        """
        read cpl file ... 
        """
        filename = '%s\\cpl.opt'%self.workdir
        
        ## Step 1: read the timeframe info: JDAY
        f = open(filename, 'r')
        JDAYs = []
        for s in f:
            line = s.split()
            if 'ZONE' in line:
                #print 'Reading Model Run Time JDAY = %s ... \n'%line[2]
                #pdb.set_trace()
                JDAYs.append(line[2][:-1])
        f.close()
        
        
        ## Step 2: read the output variables at each time step
        backflow_sections = []
        recording = False
        f = open(filename, 'r')
        for i in range(len(JDAYs)):
            ttem = JDAYs[i]
            new_section = []
            for s in f:
                line = s.split()
                #pdb.set_trace()
                if recording is False:
                    if ttem+'"' in line:
                        #pdb.set_trace()
                        print 'Reading Model Run Time JDAY = %s ... \n'%ttem
                        recording = True
                        #pdb.set_trace()
                        new_section.append(line)
                elif recording is True:
                    new_section.append(line)
                    if 'TEXT' in line:
                        recording = False
                        break
            
            #if i == len(JDAYs)-1:
            backflow_sections.append(new_section[1:-1])
            #else:
            #    self.backflow_sections.append(new_section[1:-1])
        
        ## Step 3: get time info     
        self.runtimes = self.JDAY_con(JDAYs)
        
        #pdb.set_trace()
        ##Step 4: get X and Z coordinates    
        self.X_flow = []
        self.Z_flow = []
        self.U = []
        self.W = []
        T = []
        RHO = []
        TDS = []
        tracer = []
        
        for i in range(len(backflow_sections)):
            section = backflow_sections[i]
            xtem = [float(l[0]) for l in section]
            ztem = [float(p[1]) for p in section]
            utem = [float(r[2]) for r in section]
            wtem = [float(q[3]) for q in section]
            Ttem = [float(t[4]) for t in section]
            rhotem = [float(rho[5]) for rho in section]
            tdstem = [float(tds[6]) for tds in section]
            
            try:
                trtem = [float(tr[7]) for tr in section]
            except:
                ttrrtem = []
                for tr in section:
                    try:
                        ttrrtem.append(float(tr[7]))
                    except:
                        ttrrtem.append(0)
                trtem.append(ttrrtem)
                
            self.X_flow.append(xtem)
            self.Z_flow.append(ztem)
            self.U.append(utem)
            self.W.append(wtem)
            T.append(Ttem)
            RHO.append(rhotem)
            TDS.append(tdstem)
            tracer.append(trtem)
            
        ## create python dictionary for water quality variables
        self.var_output = {}
        self.var_output.update({'TDS':{'value': TDS,'limits':[140,200], 'long_name': 'Total dissolved solids'}})
        self.var_output.update({'Tracer':{'value': tracer,'limits':[0,2], 'long_name': 'Conservative tracer'}})
        self.var_output.update({'T':{'value': T,'limits':[0,35], 'long_name': 'Water temperature (C)'}})
        
        #pdb.set_trace()
        
        
    def VisContour(self, varname='Tracer', timestep=-1, branchID=1, Plotuv=False, PlotGrid=False):
        """
        Create the contour plot of a variable given the variable name
        variable names are provided in cpl.opt
        e.g. 'U', 'W', 'Tracer'
        The python dictionary for variables can be found in vardict.py
        """
        self.Readcpl()
        
        
        ## search for index for each branch 
        ## algorithm find the distance between two elements in self.X_flow that are large, 
        ## this is where the two branches separate
        #dist = np.diff(self.X_flow[0])
        dist = np.diff(self.X_flow[timestep])
        inds = np.where(dist>1200)[0]
        
        ## this way is not necessary but clear
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
            #ind1 = len(self.X_flow[0])
            ind1 = len(self.X_flow[timestep])
            
            
            
        X_flow = self.X_flow[timestep][ind0:ind1+1]
        Z_flow = self.Z_flow[timestep][ind0:ind1+1]
        
        X_flow = np.asarray(X_flow)
        Z_flow = np.asarray(Z_flow)


        plt.rcParams.update({'font.size': 18})
        fig = plt.figure(figsize=(11.5,8))
        ax = fig.add_subplot(111)
        
        
        if PlotGrid == True:
            from bathymetry import W2_Bathymetry
            filename = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
            WB = W2_Bathymetry(filename)
            pat = WB.VisBranch2(branchID)
            for sq in pat:
                ax.add_patch(sq)
            ax.autoscale_view()
            
            ## align each branch with grid 
            ## not sure how X is defined in the model 
            ## algorithm starting from the end of the branch  
            #pdb.set_trace()
            dx =  WB.X.max() - X_flow.max()
            X_flow += dx
            
            #pdb.set_trace()
        
        
        if Plotuv:
            #X_flow = np.asarray(X_flow)
            #Z_flow = np.asarray(Z_flow)
            
            U = self.U[timestep][ind0:ind1+1]
            W = self.W[timestep][ind0:ind1+1]
            U = np.asarray(U)
            W = np.asarray(W)
            mask = np.logical_or(U != self.mask_value,W != self.mask_value)
            #U = np.ma.masked_array(U,mask=U==self.mask_value)
            #W = np.ma.masked_array(W,mask=W==self.mask_value)
            
            scale = 1.
            scale = 100./scale
            Q = ax.quiver(X_flow[mask], Z_flow[mask], np.asarray(U[mask])*100., np.asarray(W[mask])*100.,
                              zorder=5, width=0.001, headwidth=4, headlength=4.5,
                              scale=scale, color='r')
            qk = ax.quiverkey(Q, 0.15, 0.15, 1, r'$1 \frac{cm}{s}$', labelpos='W',fontproperties={'weight': 'bold','size':20})
        
        
        var = self.var_output[varname]['value'][timestep][ind0:ind1+1]
        var = np.asarray(var)
        #tracer[tracer==self.mask_value] = np.nan
        var = np.ma.masked_array(var,mask=var==self.mask_value)
        
        #xgrid, zgrid = np.meshgrid(X_flow, Z_flow)
        #Tgrid = griddata((X_flow, Z_flow), T, (xgrid, zgrid))
        #ax.contourf(xgrid, zgrid, Tgrid, 10, cmap=plt.cm.bone)
        #levels = np.linspace(self.var_output[varname]['limits'][0], self.var_output[varname]['limits'][1], 100)
        levels = np.linspace(var.min(), var.max(), 100)
        #pdb.set_trace()
        cmap = plt.set_cmap('bone_r')
        cs = ax.tricontourf(X_flow, Z_flow, var, cmap=cmap, levels=levels)
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        cb = fig.colorbar(cs, cax=cax, orientation='vertical')
        cb.ax.tick_params(labelsize=12)
        cb.ax.yaxis.offsetText.set_fontsize(12)
        cb.set_label('%s'%self.var_output[varname]['long_name'], fontsize=14)
            
            
        timestr = datetime.strftime(self.runtimes[timestep],'%Y-%m-%d')
        ax.title.set_text('Time: %s'%timestr)
        #ax.set_xlim([0, 60000])
        ax.set_ylim([135, 160])
        ax.set_xlabel('Distance from upstream (m)')
        ax.set_ylabel('Water Depth (m)')
        #ax.yaxis.grid(True)
        #ax.xaxis.grid(True)
        plt.show()
        #plt.savefig('%s\\%s_%s_%s.png'%(self.workdir, varname, str(branchID), str(timestep)))
        #plt.savefig('example.png')
        #plt.close()
        
    
    def AnimateContour(self, varname='Tracer', branchID=1, PlotGrid=False):
        """
        create animation 
        """
        import matplotlib.animation as animation
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        
        
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)
        
        self.Readcpl()
        
        plt.rcParams.update({'font.size': 18})
        fig = plt.figure(figsize=(12.5,8))
        ax = fig.add_subplot(111)
        
        if PlotGrid == True:
            from bathymetry import W2_Bathymetry
            filename = '%s\\%s'%(self.workdir, 'Bth_WB1.npt')
            WB = W2_Bathymetry(filename)
            pat = WB.VisBranch2(branchID)
            
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        
        def animate(ii):
            ax.clear()
            ## search for index for each branch 
            ## algorithm find the distance between two elements in self.X_flow that are large, 
            ## this is where the two branches separate
            dist = np.diff(self.X_flow[ii])
            inds = np.where(dist>1200)[0]
            ## this way is not necessary but clear
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
                ind1 = len(self.X_flow[ii])
                
            X_flow = self.X_flow[ii][ind0:ind1+1]
            Z_flow = self.Z_flow[ii][ind0:ind1+1]
            X_flow = np.asarray(X_flow)
            Z_flow = np.asarray(Z_flow)
            
            if PlotGrid == True: 
                for sq in pat:
                    ax.add_patch(sq)
                ax.autoscale_view()  
                ## align each branch with grid 
                dx =  WB.X.max() - X_flow.max()
                X_flow += dx
            
            U = self.U[ii][ind0:ind1+1]
            W = self.W[ii][ind0:ind1+1]
            U = np.asarray(U)
            W = np.asarray(W)
            mask = np.logical_or(U != self.mask_value,W != self.mask_value)
            
            scale = 1.
            scale = 100./scale
            Q = ax.quiver(X_flow[mask], Z_flow[mask], np.asarray(U[mask])*100., np.asarray(W[mask])*100.,
                              zorder=5, width=0.001, headwidth=4, headlength=4.5,
                              scale=scale, color='r')
            qk = ax.quiverkey(Q, 0.15, 0.15, 1, r'$1 \frac{cm}{s}$', labelpos='W',fontproperties={'weight': 'bold','size':20})
            
            
            var = self.var_output[varname]['value'][ii][ind0:ind1+1]
            var = np.asarray(var)
            var = np.ma.masked_array(var,mask=var==self.mask_value)
            levels = np.linspace(var.min(), var.max(), 100)
            cmap = plt.set_cmap('bone_r')
            cs = ax.tricontourf(X_flow, Z_flow, var, cmap=cmap, levels=levels)
            
        
            cb = fig.colorbar(cs, cax=cax, orientation='vertical')
            cb.ax.tick_params(labelsize=12)
            cb.ax.yaxis.offsetText.set_fontsize(12)
            cb.set_label('%s'%self.var_output[varname]['long_name'], fontsize=14)
            
            
            timestr = datetime.strftime(self.runtimes[ii],'%Y-%m-%d')
            ax.title.set_text('Time: %s'%timestr)
            ax.set_ylim([135, 160])
            ax.set_xlabel('Distance from upstream (m)')
            ax.set_ylabel('Water Depth (m)')
            
            return cs, cax
            
            
        anim = animation.FuncAnimation(fig, animate, frames=100, interval=600, blit=False)
        anim.save('%s\\tracer.mp4'%self.workdir, writer=writer)

        
        
        
    def JDAY_con(self,JDAY):
        """
        funtion adapted from Justin's code
        """
        basedate = datetime(2010,12,31)
        return_date = [basedate+timedelta(days = float(ii)) for ii in JDAY]
        return return_date
        
#### For testing ####       
        
if __name__ == "__main__":
    #wdir = r'C:\Users\dfeng\Downloads\v42\Tests\20191112_tracer'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191113_tracer_test'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191127_1115_tracer_test'
    wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191127_1336_tracer_test'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191127_1631_tracer_test'
    #wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\tracer_test\20191202_1100_tracer_test'
    WC = W2_Contour(wdir)
    #WC.VisContour('Tracer', timestep=35, branchID=1, Plotuv=True, PlotGrid=True)
    WC.AnimateContour('Tracer', branchID=1, PlotGrid=True)