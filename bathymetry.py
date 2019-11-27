# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 09:30:56 2019

@author: dfeng
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.collections as coll

import pdb

class W2_Bathymetry(object):
    """
    general class for reading and visualizing CE-QUAL-W2 model bathymetry
    """
    Nlyr = 34
    Nseg = 122
    
    def __init__(self, bathfile, **kwargs):
        self.__dict__.update(kwargs)
        
        self.bathfile = bathfile
        
    
    def readVar(self): 
        
        self.seg_length = self.readSeg('length')
        self.seg_ori = self.readSeg('Orientation') 
        self.lyr_height = self.readSeg('Height')
        
        
        #pdb.set_trace()
        
        
    def readSeg(self, varname='length'):
        """
        funtion that reads W2 bathymetry basic info
        """
        output_section = []
        recording = False
        f = open(self.bathfile, 'r')
        for s in f:
            line = s.split()
            if recording is False:
                if varname in line:
                    #print "Reading Segment %s ...\n"%varname
                    recording = True
            elif recording is True:
                output_section.append(line)
                if len(line) == 0:
                    recording = False
                    break
        f.close()
        
        return [float(item) for sublist in output_section[:-1] for item in sublist]
    
    
    def readLyr(self):
        """
        funtion that reads W2 segment layer info
        """
        
        self.lyrs = np.zeros([self.Nseg, self.Nlyr])
        self.lyrs_branchID = []
        self.lyrs_segID = []
        
        ## Step 1: read segment ID    
        f = open(self.bathfile, 'r')
        for s in f:
            line = s.split()
            if 'Branch' in line:
                #print line
                #print "Reading Segment %s ...\n"%line 
                if len(line) == 6:
                    self.lyrs_branchID.append(line[3])
                    self.lyrs_segID.append(line[5])
                elif len(line) == 5:
                    self.lyrs_branchID.append(line[2])
                    self.lyrs_segID.append(line[4])
                elif len(line) == 4:
                    self.lyrs_branchID.append(line[1])
                    self.lyrs_segID.append(line[3])
                
        f.close()
        
        ## Step 2: read Layer from each segment
        for i in range(len(self.lyrs_segID)):
            segID = self.lyrs_segID[i]
            lyr_tem = self.readSegLyr(segID)  
            #print lyr_tem
            self.lyrs[i,:len(lyr_tem)] = lyr_tem
            
        #pdb.set_trace()
        
    
    def readSegLyr(self, segID):
        """
        function that reads the Segment layer given the Seg ID
        """
        output_section = []
        recording = False
        f = open(self.bathfile, 'r')
        for s in f:
            line = s.split()
            if recording is False:
                if 'Branch' in line and segID in line[3:]:
                    #print "Reading Segment %s ...\n"%segID
                    recording = True
            elif recording is True:
                output_section.append(line)
                if len(line) == 0:
                    recording = False
                    break
        f.close()
        
        #pdb.set_trace()
        return [float(item) for sublist in output_section[:-1] for item in sublist]
    
    
    def VisBranch(self, branchID=1):
        """
        Given the branch ID, plot the grid at that branch
        """
        
        ## read Segment info
        self.seg_length = self.readSeg('length')
        self.lyr_height = self.readSeg('Height')
        
        ## read layer info for all segments
        self.readLyr()
        
        lyrs_branchID = np.asarray([int(float(ID)) for ID in self.lyrs_branchID])
        lyrs_segID = np.asarray([int(float(ID)) for ID in self.lyrs_segID])
        
        ## segments in the specified branch
        segs = lyrs_segID[lyrs_branchID==branchID]
        seg_length = np.asarray(self.seg_length)[lyrs_branchID==branchID]
        
        lyrs = self.lyrs[lyrs_branchID==branchID,:]
        #pdb.set_trace()
        ############################################
        ## The bottom and left rectangle coordinates
        X = np.cumsum(seg_length) 
        Z = np.cumsum(np.asarray(self.lyr_height)) 
        
        X = X.tolist()
        Z = Z.tolist()
        
        X = [0] + X[:-1]
        Z = [0] + Z[:-1]
        
        Z = [156.68 - zz for zz in Z]  ## check W2ControlGUI
        ############################################
        
        
        plt.rcParams.update({'font.size': 18})
        fig = plt.figure(figsize=(11.5,8))
        ax = fig.add_subplot(111)
        
#        pat = []
        for i in range(len(X)):
            for j in range(len(Z)):
#                 ax.axvline(x=X[i], color='k')
#                 ax.axhline(y=Z[j], color='k')
                if lyrs[i,j] == 0:  ## masked grid
                    sq = patches.Rectangle((X[i], Z[j]), seg_length[i], self.lyr_height[j], fill=True)
                elif lyrs[i,j] != 0:  ## computational grid
                    sq = patches.Rectangle((X[i], Z[j]), seg_length[i], self.lyr_height[j], fill=False)
#                pat.append(sq)
                ax.add_patch(sq)

        #pc = coll.PatchCollection(pat, facecolor=None, edgecolor='k')
        #ax.add_collection(pc)
        #ax.relim()
        ax.set_ylim([135, 160])
        ax.title.set_text('Branch %s'%str(branchID))
        ax.set_xlabel('Distance from upstream (m)')
        ax.set_ylabel('Water Depth (m)')
        ax.autoscale_view()
        #plt.gca().invert_yaxis()
        #plt.show()
        plt.savefig('figures_grids\grid_branch_%s.png'%str(branchID))
        plt.close()
        #pdb.set_trace()
        
    def VisBranch2(self, branchID=1):
        """
        Given the branch ID, provide the grid to align with velocity output
        """
        
        ## read Segment info
        self.seg_length = self.readSeg('length')
        self.lyr_height = self.readSeg('Height')
        
        ## read layer info for all segments
        self.readLyr()
        
        lyrs_branchID = np.asarray([int(float(ID)) for ID in self.lyrs_branchID])
        lyrs_segID = np.asarray([int(float(ID)) for ID in self.lyrs_segID])
        
        ## segments in the specified branch
        segs = lyrs_segID[lyrs_branchID==branchID]
        seg_length = np.asarray(self.seg_length)[lyrs_branchID==branchID]
        
        lyrs = self.lyrs[lyrs_branchID==branchID,:]
        #pdb.set_trace()
        ############################################
        ## The bottom and left rectangle coordinates
        X = np.cumsum(seg_length) 
        Z = np.cumsum(np.asarray(self.lyr_height)) 
        
        X = X.tolist()
        Z = Z.tolist()
        
        X = [0] + X[:-1]
        Z = [0] + Z[:-1]
        
        #Z = [148.3 - zz for zz in Z]
        if branchID == 1:
            X = [xx - 862.5 for xx in X]
        elif branchID == 2:
            X = [xx - 652.7 for xx in X]   ## starting from segment 47
        elif branchID == 3:
            X = [xx - 788.6 for xx in X]   ## starting from segment 62
        elif branchID == 4:
            X = [xx - 1132 for xx in X]   ## starting from segment 76
        elif branchID == 5:
            X = [xx - 673.8 for xx in X]   ## starting from segment 86
            
            
        
        #Z = [156.07 - 0.61/2. - zz for zz in Z]  ## check W2ControlGUI
        Z = [156.07 - zz for zz in Z]
        
        ############################################
        pat = []
        for i in range(len(X)):
            for j in range(len(Z)):
#                 ax.axvline(x=X[i], color='k')
#                 ax.axhline(y=Z[j], color='k')
                if lyrs[i,j] == 0:  ## masked grid
                    sq = patches.Rectangle((X[i], Z[j]), seg_length[i], self.lyr_height[j], linewidth=0.1, fill=True)
                elif lyrs[i,j] != 0:  ## computational grid
                    sq = patches.Rectangle((X[i], Z[j]), seg_length[i], self.lyr_height[j], linewidth=0.1, fill=False)
                pat.append(sq)
#                ax.add_patch(sq)
        
        
        return pat
        

        
#### For testing ####       
        
if __name__ == "__main__":
    
    filename = r'C:\Users\dfeng\Downloads\v42\Tests\20191111_baseline3\Bth_WB1.npt'
    WB = W2_Bathymetry(filename)
    WB.VisBranch(1)
    