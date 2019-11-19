# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 11:50:15 2019

@author: dfeng
"""

"""
BRANCH G      US      DS     UHS     DHS     UQB     DQB   NLMIN   SLOPE  SLOPEC
BR1            2      45       0      00       0       0       1 0.00000     0.0
BR2           48      59       0      33       0       0       1 0.00000     0.0
BR3           62      73       0      34       0       0       1 0.00000     0.0
BR4           76      83       0      42       0       0       1 0.00000     0.0
BR5           86     121       0      43       0       0       1 0.00000     0.0
"""

import numpy as np
import matplotlib.pyplot as plt

from bathymetry import W2_Bathymetry

import pdb


class W2_Segmentation(object):
    """
    general class for reading and visualizing CE-QUAL-W2 model segmentation (horizontal view of grids)
    This is important for visualizing particles in the horizontal view
    """
    
    Nlyr = 34
    Nseg = 122
    
    def __init__(self, bathfile, **kwargs):
        self.__dict__.update(kwargs)
        
        self.bathfile = bathfile
        
    def readBathymetry(self):
        """
        read segment and layer info from the bathymetry file
        """
        WB = W2_Bathymetry(self.bathfile)
        WB.readVar()
        WB.readLyr()
        
        self.seg_length = WB.seg_length   #segment length
        self.seg_ori = WB.seg_ori    # segment orientation
        self.lyr_height = WB.lyr_height  # layer height or each segment
        
        self.lyrs_branchID = np.asarray([int(float(ID)) for ID in WB.lyrs_branchID])
        self.lyrs_segID = np.asarray([int(float(ID)) for ID in WB.lyrs_segID])
        self.lyrs = WB.lyrs
        
    
    def VisSeg(self):
        """
        visualize the alignment of segments
        """
        
        self.readBathymetry()
        
        segs, Pnts, segs2, Pnts2, segs3, Pnts3, segs4, Pnts4, segs5, Pnts5 = self.Seg2Pnt()
        
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111)
        
        ax.plot(Pnts[:,0], Pnts[:,1], '.-k')
        ax.plot(Pnts2[:,0], Pnts2[:,1], '.-b')
        ax.plot(Pnts3[:,0], Pnts3[:,1], '.-b')
        ax.plot(Pnts4[:,0], Pnts4[:,1], '.-b')
        ax.plot(Pnts5[:,0], Pnts5[:,1], '.-b')
        
        for i in range(len(segs)):
            ax.annotate('%s'%str(segs[i]), (Pnts[i,0], Pnts[i,1]), color='r', fontsize=10)
        for i in range(len(segs2)):
            ax.annotate('%s'%str(segs2[i]), (Pnts2[i,0], Pnts2[i,1]), color='r', fontsize=10)
        for i in range(len(segs3)):
            ax.annotate('%s'%str(segs3[i]), (Pnts3[i,0], Pnts3[i,1]), color='r', fontsize=10)
        for i in range(len(segs4)):
            ax.annotate('%s'%str(segs4[i]), (Pnts4[i,0], Pnts4[i,1]), color='r', fontsize=10)
        for i in range(len(segs5)):
            ax.annotate('%s'%str(segs5[i]), (Pnts5[i,0], Pnts5[i,1]), color='r', fontsize=10)
        
        ax.set_xlim([-20000,10000])
        ax.set_aspect(True)
        plt.show()
        
        #pdb.set_trace()
    
    def Seg2Pnt(self):
        """
        convert the segment information to points
        """
        ## tracing points 
        ## determine the starting point coordinates (0,0)
        ## then the next point is 
        
        segs, seg_length, seg_ori, lyrs = self.SegInfo(1)
        
        Pnts1 = [[0,0]]
            
        for i in range(len(segs)):
            #pdb.set_trace()
            Ptem = [Pnts1[i][0]-seg_length[i]*np.sin(seg_ori[i]), Pnts1[i][1]-seg_length[i]*np.cos(seg_ori[i])]
            Pnts1.append(Ptem)
        
        pdb.set_trace()
        segs2, Pnts2 = self.Seg2Pnt_branch(Pnts1, branchID=2)
        
        segs3, Pnts3 = self.Seg2Pnt_branch(Pnts1, branchID=3)
        
        segs4, Pnts4 = self.Seg2Pnt_branch(Pnts1, branchID=4)
        
        segs5, Pnts5 = self.Seg2Pnt_branch(Pnts1, branchID=5)
        
        
        return segs[1:-1], np.asarray(Pnts1)[1:-1], segs2, Pnts2, segs3, Pnts3, segs4, Pnts4, segs5, Pnts5
        
    
    def Seg2Pnt_branch(self, Pnts1, branchID):
        """
        for branches other than the main branch
        """
        
        segs1, seg_length1, seg_ori1, lyrs1 = self.SegInfo(1)
        
        segs, seg_length, seg_ori, lyrs = self.SegInfo(branchID)
        ## delete boundary cells first
        segs = segs[1:-1]
        seg_length = seg_length[1:-1]
        seg_ori = seg_ori[1:-1]
    
        ## Step 1: find the coordinates of the starting point for each branch
        if branchID == 2:
            DHS=33
            ## startting segment coordination
            ssc = Pnts1[DHS-1]
            Pnts = [ [ ssc[0]+lyrs1[DHS-1]*np.cos(np.pi*2-seg_ori1[DHS-1]), ssc[1]+lyrs1[DHS-1]*np.sin(np.pi*2-seg_ori1[DHS-1]) ] ]
#            ns = len(segs)  
#            for i in range(len(segs)):
#                Ptem = [Pnts[i][0]+seg_length[ns-1-i]*np.sin(seg_ori[ns-1-i]), Pnts[i][1]+seg_length[ns-1-i]*np.cos(seg_ori[ns-1-i])]
#                Pnts.append(Ptem)
            
        elif branchID == 3:
            DHS=34
            ssc = Pnts1[DHS-1]
            Pnts = [ [ ssc[0]-lyrs1[DHS-1]*np.cos(np.pi*2-seg_ori1[DHS-1]), ssc[1]-lyrs1[DHS-1]*np.sin(np.pi*2-seg_ori1[DHS-1]) ] ]
#            ns = len(segs)
#            for i in range(len(segs)):
#                Ptem = [Pnts[i][0]+seg_length[ns-1-i]*np.sin(seg_ori[ns-1-i]), Pnts[i][1]+seg_length[ns-1-i]*np.cos(seg_ori[ns-1-i])]
#                Pnts.append(Ptem)
        elif branchID == 4:
            DHS=42
            ssc = Pnts1[DHS-1]
            Pnts = [ [ ssc[0]+lyrs1[DHS-1]*np.cos(np.pi*2-seg_ori1[DHS-1]), ssc[1]+lyrs1[DHS-1]*np.sin(np.pi*2-seg_ori1[DHS-1]) ] ]
#            ns = len(segs)
#            for i in range(len(segs)):
#                Ptem = [Pnts[i][0]+seg_length[ns-1-i]*np.sin(seg_ori[ns-1-i]), Pnts[i][1]+seg_length[ns-1-i]*np.cos(seg_ori[ns-1-i])]
#                Pnts.append(Ptem)
        elif branchID == 5:
            DHS=43
            ssc = Pnts1[DHS-1]
            Pnts = [ [ ssc[0]-lyrs1[DHS-1]*np.cos(np.pi*2-seg_ori1[DHS-1]), ssc[1]-lyrs1[DHS-1]*np.sin(np.pi*2-seg_ori1[DHS-1]) ] ]
#            ns = len(segs)
#            for i in range(len(segs)):
#                Ptem = [Pnts[i][0]+seg_length[ns-1-i]*np.sin(seg_ori[ns-1-i]), Pnts[i][1]+seg_length[ns-1-i]*np.cos(seg_ori[ns-1-i])]
#                Pnts.append(Ptem)
    
        ns = len(segs)  
        for i in range(len(segs)):
            Ptem = [Pnts[i][0]+seg_length[ns-1-i]*np.sin(seg_ori[ns-1-i]), Pnts[i][1]+seg_length[ns-1-i]*np.cos(seg_ori[ns-1-i])]
            Pnts.append(Ptem)

        
        #pdb.set_trace()
        return segs[::-1], np.asarray(Pnts)
    
    
        
    def SegInfo(self, branchID=1):
        """
        Input branchID, output the segmentation information for that branch
        """
        segs = self.lyrs_segID[self.lyrs_branchID==branchID]
        seg_length = np.asarray(self.seg_length)[self.lyrs_branchID==branchID]
        seg_ori = np.asarray(self.seg_ori)[self.lyrs_branchID==branchID]
        
        lyrs = self.lyrs[self.lyrs_branchID==branchID]
#        seg_ori_degree = []
#        for ori in seg_ori:
#            seg_ori_degree.append(np.degrees(ori))
        
        return segs, seg_length, seg_ori, lyrs[:,1]
    
    
#### For testing ####       
        
if __name__ == "__main__":
    
    filename = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\particle_tracking_test\20191111_baseline4\Bth_WB1.npt'
    WB = W2_Segmentation(filename)
    WB.VisSeg()