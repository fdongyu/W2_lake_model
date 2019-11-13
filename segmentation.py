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
        
        branchID=1
        
        segs = self.lyrs_segID[self.lyrs_branchID==branchID]
        seg_length = np.asarray(self.seg_length)[self.lyrs_branchID==branchID]
        seg_ori = np.rad2deg(np.asarray(self.seg_ori))[self.lyrs_branchID==branchID]
        
        seg_ori_new = []
        for ori in seg_ori:
            if ori<90:
                seg_ori_new.append(np.radians(360-ori))
            else:
                seg_ori_new.append(np.radians(ori))
        
        pdb.set_trace()
        
    
    
#### For testing ####       
        
if __name__ == "__main__":
    
    filename = r'C:\Users\dfeng\Downloads\v42\Tests\20191111_baseline3\Bth_WB1.npt'
    WB = W2_Segmentation(filename)
    WB.VisSeg()