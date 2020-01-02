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

import os
import numpy as np
import matplotlib.pyplot as plt
import utm

from bathymetry import W2_Bathymetry

import pdb


class W2_Segmentation(object):
    """
    general class for reading and visualizing CE-QUAL-W2 model segmentation (horizontal view of grids)
    This is important for visualizing particles in the horizontal view
    """
    
    Nlyr = 34
    Nseg = 122
    
    #LatLon = [33.238, -96.41]   ## This is for modified smooth segments
    LatLon = [33.241, -96.415]
    
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
        
        ## Step One: plot the segment link
        segs, Pnts, segs2, Pnts2, segs3, Pnts3, segs4, Pnts4, segs5, Pnts5 = self.Seg2Pnt()
        ## Note the length of the segs = length of Pnts + 1. and the segment is not points, it is the line between two points
        
        ## Step Two: each segment
        eastPnts, westPnts = self.SegLyr(Pnts, branchID=1)
        eastPnts2, westPnts2 = self.SegLyr(Pnts2, branchID=2)
        eastPnts3, westPnts3 = self.SegLyr(Pnts3, branchID=3)
        eastPnts4, westPnts4 = self.SegLyr(Pnts4, branchID=4)
        eastPnts5, westPnts5 = self.SegLyr(Pnts5, branchID=5)
        
        
        plt.rcParams.update({'font.size': 16})
        fig = plt.figure(figsize=(10,12))
        ax = fig.add_subplot(111)
        
        ax.plot(Pnts[:,0], Pnts[:,1], '.-k')
        ax.plot(Pnts2[:,0], Pnts2[:,1], '.-b')
        ax.plot(Pnts3[:,0], Pnts3[:,1], '.-b')
        ax.plot(Pnts4[:,0], Pnts4[:,1], '.-b')
        ax.plot(Pnts5[:,0], Pnts5[:,1], '.-b')
        
        
        for i in range(len(westPnts)):
            ax.plot([westPnts[i,0], eastPnts[i,0]], [westPnts[i,1], eastPnts[i,1]], '-k')
        for i in range(len(westPnts2)):
            ax.plot([westPnts2[i,0], eastPnts2[i,0]], [westPnts2[i,1], eastPnts2[i,1]], '-k')
        for i in range(len(westPnts3)):
            ax.plot([westPnts3[i,0], eastPnts3[i,0]], [westPnts3[i,1], eastPnts3[i,1]], '-k')
        for i in range(len(westPnts4)):
            ax.plot([westPnts4[i,0], eastPnts4[i,0]], [westPnts4[i,1], eastPnts4[i,1]], '-k')
        for i in range(len(westPnts5)):
            ax.plot([westPnts5[i,0], eastPnts5[i,0]], [westPnts5[i,1], eastPnts5[i,1]], '-k')
        
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
        
        #ax.set_xlim([-20000,10000])
        ax.set_aspect(True)
        
        fig.tight_layout()
        plt.show()
        #plt.savefig('segmentations.png')
        #plt.close()
        
        
        #pdb.set_trace()
        
    def VisSeg2(self):
        
        self.readBathymetry()
        
        ## Step One: plot the segment link
        self.segs1,  self.Pnts1,  self.segs2,  self.Pnts2,  self.segs3,  self.Pnts3,  self.segs4,  self.Pnts4,  self.segs5,  self.Pnts5 = self.Seg2Pnt()
        ## Note the length of the segs = length of Pnts + 1. and the segment is not points, it is the line between two points
        
        ## Step Two: each segment
        self.eastPnts1, self.westPnts1 = self.SegLyr(self.Pnts1, branchID=1)
        self.eastPnts2, self.westPnts2 = self.SegLyr(self.Pnts2, branchID=2)
        self.eastPnts3, self.westPnts3 = self.SegLyr(self.Pnts3, branchID=3)
        self.eastPnts4, self.westPnts4 = self.SegLyr(self.Pnts4, branchID=4)
        self.eastPnts5, self.westPnts5 = self.SegLyr(self.Pnts5, branchID=5)
        
        
    def VisSeg_shp(self):
        """
        create shapefile to visualize the W2 segments
        """
        import shapefile as shp
        import utm
        
        self.readBathymetry()
        
        ## Step One: plot the segment link
        self.segs1,  self.Pnts1,  self.segs2,  self.Pnts2,  self.segs3,  self.Pnts3,  self.segs4,  self.Pnts4,  self.segs5,  self.Pnts5 = self.Seg2Pnt()
        ## Note the length of the segs = length of Pnts + 1. and the segment is not points, it is the line between two points
        
        ## Step Two: each segment
        self.eastPnts1, self.westPnts1 = self.SegLyr(self.Pnts1, branchID=1)
        self.eastPnts2, self.westPnts2 = self.SegLyr(self.Pnts2, branchID=2)
        self.eastPnts3, self.westPnts3 = self.SegLyr(self.Pnts3, branchID=3)
        self.eastPnts4, self.westPnts4 = self.SegLyr(self.Pnts4, branchID=4)
        self.eastPnts5, self.westPnts5 = self.SegLyr(self.Pnts5, branchID=5)
        
        w = shp.Writer('shp\W2_segments')
        w.field('branch_ID','C')
        w.field('segment_ID','C')
        w.field('Lon_west','C')  #float - needed for coordinates
        w.field('Lat_west','C') 
        w.field('Lon_east','C')
        w.field('Lat_east','C')
        

        #lines = []
        for i in range(len(self.westPnts1)):
            westlat, westlon = utm.to_latlon(self.westPnts1[i,0], self.westPnts1[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(self.eastPnts1[i,0], self.eastPnts1[i,1], 14, 'U')
            #lines.append([[westlat, westlon], [eastlat, eastlon]])
            w.line([[[westlon, westlat], [eastlon, eastlat]]])
            w.record('branch1', str(self.segs1[i]), "{:.3f}".format(westlon), "{:.3f}".format(westlat), \
                     "{:.3f}".format(eastlon), "{:.3f}".format(eastlat))
            
        for i in range(len(self.westPnts2)):
            westlat, westlon = utm.to_latlon(self.westPnts2[i,0], self.westPnts2[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(self.eastPnts2[i,0], self.eastPnts2[i,1], 14, 'U')
            w.line([[[westlon, westlat], [eastlon, eastlat]]])
            w.record('branch2', str(self.segs2[::-1][i]), "{:.3f}".format(westlon), "{:.3f}".format(westlat), \
                     "{:.3f}".format(eastlon), "{:.3f}".format(eastlat))
            
        for i in range(len(self.westPnts3)):
            westlat, westlon = utm.to_latlon(self.westPnts3[i,0], self.westPnts3[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(self.eastPnts3[i,0], self.eastPnts3[i,1], 14, 'U')
            w.line([[[westlon, westlat], [eastlon, eastlat]]])
            w.record('branch3', str(self.segs3[::-1][i]), "{:.3f}".format(westlon), "{:.3f}".format(westlat), \
                     "{:.3f}".format(eastlon), "{:.3f}".format(eastlat))
        
        for i in range(len(self.westPnts4)):
            westlat, westlon = utm.to_latlon(self.westPnts4[i,0], self.westPnts4[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(self.eastPnts4[i,0], self.eastPnts4[i,1], 14, 'U')
            w.line([[[westlon, westlat], [eastlon, eastlat]]])
            w.record('branch4', str(self.segs4[::-1][i]), "{:.3f}".format(westlon), "{:.3f}".format(westlat), \
                     "{:.3f}".format(eastlon), "{:.3f}".format(eastlat))
            
        for i in range(len(self.westPnts5)):
            westlat, westlon = utm.to_latlon(self.westPnts5[i,0], self.westPnts5[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(self.eastPnts5[i,0], self.eastPnts5[i,1], 14, 'U')
            w.line([[[westlon, westlat], [eastlon, eastlat]]])
            w.record('branch5', str(self.segs5[::-1][i]), "{:.3f}".format(westlon), "{:.3f}".format(westlat), \
                     "{:.3f}".format(eastlon), "{:.3f}".format(eastlat))
            
        w.close()
            
        prj = open("shp\W2_segments.prj", "w") 
        epsg = 'GEOGCS["WGS 84",'
        epsg += 'DATUM["WGS_1984",'
        epsg += 'SPHEROID["WGS 84",6378137,298.257223563]]'
        epsg += ',PRIMEM["Greenwich",0],'
        epsg += 'UNIT["degree",0.0174532925199433]]'
        prj.write(epsg)
        prj.close()
            
        
    
    def VisSeg_shp_GDAL(self):  
        """
        create segment shapefile using GDAL
        LineCollection(segs, linewidths=(0.5, 1, 1.5, 2),
                               colors=colors, linestyle='solid')
        """
        from matplotlib.collections import LineCollection
        import matplotlib
        import utm
        
        self.readBathymetry()
        
        ## Step One: plot the segment link
        self.segs1,  self.Pnts1,  self.segs2,  self.Pnts2,  self.segs3,  self.Pnts3,  self.segs4,  self.Pnts4,  self.segs5,  self.Pnts5 = self.Seg2Pnt()
        ## Note the length of the segs = length of Pnts + 1. and the segment is not points, it is the line between two points
        
        ## Step Two: each segment
        self.eastPnts1, self.westPnts1 = self.SegLyr(self.Pnts1, branchID=1)
        self.eastPnts2, self.westPnts2 = self.SegLyr(self.Pnts2, branchID=2)
        self.eastPnts3, self.westPnts3 = self.SegLyr(self.Pnts3, branchID=3)
        self.eastPnts4, self.westPnts4 = self.SegLyr(self.Pnts4, branchID=4)
        self.eastPnts5, self.westPnts5 = self.SegLyr(self.Pnts5, branchID=5)
        
        
        latlon = dict()
        latlon.setdefault('westlat1', [])
        latlon.setdefault('westlon1', [])
        latlon.setdefault('eastlat1', [])
        latlon.setdefault('eastlon1', [])
        
        latlon.setdefault('westlat2', [])
        latlon.setdefault('westlon2', [])
        latlon.setdefault('eastlat2', [])
        latlon.setdefault('eastlon2', [])
        
        latlon.setdefault('westlat3', [])
        latlon.setdefault('westlon3', [])
        latlon.setdefault('eastlat3', [])
        latlon.setdefault('eastlon3', [])
        
        latlon.setdefault('westlat4', [])
        latlon.setdefault('westlon4', [])
        latlon.setdefault('eastlat4', [])
        latlon.setdefault('eastlon4', [])
        
        latlon.setdefault('westlat5', [])
        latlon.setdefault('westlon5', [])
        latlon.setdefault('eastlat5', [])
        latlon.setdefault('eastlon5', [])
        
        lines1 = []
        for i in range(len(self.westPnts1)):
            westlat, westlon = utm.to_latlon(self.westPnts1[i,0], self.westPnts1[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(self.eastPnts1[i,0], self.eastPnts1[i,1], 14, 'U')
            lines1.append([[westlon, westlat], [eastlon, eastlat]])
            latlon['westlat1'].append(westlat)
            latlon['westlon1'].append(westlon)
            latlon['eastlat1'].append(eastlat)
            latlon['eastlon1'].append(eastlon)
            
        lines2 = []
        for i in range(len(self.westPnts2)):
            westlat, westlon = utm.to_latlon(self.westPnts2[i,0], self.westPnts2[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(self.eastPnts2[i,0], self.eastPnts2[i,1], 14, 'U')
            lines2.append([[westlon, westlat], [eastlon, eastlat]])
            latlon['westlat2'].append(westlat)
            latlon['westlon2'].append(westlon)
            latlon['eastlat2'].append(eastlat)
            latlon['eastlon2'].append(eastlon)
            
        lines3 = []
        for i in range(len(self.westPnts3)):
            westlat, westlon = utm.to_latlon(self.westPnts3[i,0], self.westPnts3[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(self.eastPnts3[i,0], self.eastPnts3[i,1], 14, 'U')
            lines3.append([[westlon, westlat], [eastlon, eastlat]])
            latlon['westlat3'].append(westlat)
            latlon['westlon3'].append(westlon)
            latlon['eastlat3'].append(eastlat)
            latlon['eastlon3'].append(eastlon)
            
        lines4 = []
        for i in range(len(self.westPnts4)):
            westlat, westlon = utm.to_latlon(self.westPnts4[i,0], self.westPnts4[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(self.eastPnts4[i,0], self.eastPnts4[i,1], 14, 'U')
            lines4.append([[westlon, westlat], [eastlon, eastlat]])
            latlon['westlat4'].append(westlat)
            latlon['westlon4'].append(westlon)
            latlon['eastlat4'].append(eastlat)
            latlon['eastlon4'].append(eastlon)
            
        lines5 = []    
        for i in range(len(self.westPnts5)):
            westlat, westlon = utm.to_latlon(self.westPnts5[i,0], self.westPnts5[i,1], 14, 'U')
            eastlat, eastlon = utm.to_latlon(self.eastPnts5[i,0], self.eastPnts5[i,1], 14, 'U')
            lines5.append([[westlon, westlat], [eastlon, eastlat]])
            latlon['westlat5'].append(westlat)
            latlon['westlon5'].append(westlon)
            latlon['eastlat5'].append(eastlat)
            latlon['eastlon5'].append(eastlon)
            
        
        lc = LineCollection(lines1)
        outfile = os.path.join(os.getcwd(), 'segments.shp') 
        
        lines = [lines1, lines2, lines3, lines4, lines5]
        
        ## write to shapefile
        self.Contour2Shp(lines, latlon, outfile)
        
        
        
    def BranchPnt(self, branchID=1):
        """
        visualize the alignment of segments
        """
        
        self.readBathymetry()
        
        ## Step One: plot the segment link
        segs1, Pnts1, segs2, Pnts2, segs3, Pnts3, segs4, Pnts4, segs5, Pnts5 = self.Seg2Pnt()
        ## Note the length of the segs = length of Pnts + 1. and the segment is not points, it is the line between two points
        
        ## Step Two: each segment
#        eastPnts1, westPnts2 = self.SegLyr(Pnts1, branchID=1)
#        eastPnts2, westPnts2 = self.SegLyr(Pnts2, branchID=2)
#        eastPnts3, westPnts3 = self.SegLyr(Pnts3, branchID=3)
#        eastPnts4, westPnts4 = self.SegLyr(Pnts4, branchID=4)
#        eastPnts5, westPnts5 = self.SegLyr(Pnts5, branchID=5)
        
        if branchID==1:
            return segs1, Pnts1
        elif branchID==2:
            return segs2, Pnts2
        elif branchID==3:
            return segs3, Pnts3
        elif branchID==4:
            return segs4, Pnts4
        elif branchID==5:
            return segs5, Pnts5
        
        
    
    def Seg2Pnt(self):
        """
        convert the segment information to points
        """
        ## tracing points 
        ## determine the starting point coordinates (0,0)
        ## then the next point is 
        
        segs, seg_length, seg_ori, lyrs = self.SegInfo(1)
        
        xy = utm.from_latlon(self.LatLon[0], self.LatLon[1])[:2]
        
        #Pnts1 = [[0,0]]
        Pnts1 = [xy]
            
        for i in range(len(segs)):
            #pdb.set_trace()
            Ptem = [Pnts1[i][0]-seg_length[i]*np.sin(seg_ori[i]), Pnts1[i][1]-seg_length[i]*np.cos(seg_ori[i])]
            Pnts1.append(Ptem)
        
        
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
        
        midPnts1 = []
        Pnts1 = np.asarray(Pnts1)
        for i in range(len(Pnts1)-1):
            tem = (Pnts1[i]+Pnts1[i+1])/2.
            midPnts1.append(tem.tolist())       
        midPnts1 = np.asarray(midPnts1)
        
        segs, seg_length, seg_ori, lyrs = self.SegInfo(branchID)
        ## delete boundary cells first
        segs = segs[1:-1]
        seg_length = seg_length[1:-1]
        seg_ori = seg_ori[1:-1]
    
        ## Step 1: find the coordinates of the starting point for each branch
        if branchID == 2:
            DHS=33
            ## startting segment coordination
            #ssc = Pnts1[DHS-1]
            ssc = midPnts1[DHS-1]
            Pnts = [ [ ssc[0]+lyrs1[DHS-1]/2.*np.cos(np.pi*2-seg_ori1[DHS-1]), ssc[1]+lyrs1[DHS-1]/2.*np.sin(np.pi*2-seg_ori1[DHS-1]) ] ]
#            ns = len(segs)  
#            for i in range(len(segs)):
#                Ptem = [Pnts[i][0]+seg_length[ns-1-i]*np.sin(seg_ori[ns-1-i]), Pnts[i][1]+seg_length[ns-1-i]*np.cos(seg_ori[ns-1-i])]
#                Pnts.append(Ptem)
            
        elif branchID == 3:
            DHS=34
            #ssc = Pnts1[DHS-1]
            ssc = midPnts1[DHS-1]
            Pnts = [ [ ssc[0]-lyrs1[DHS-1]/2.*np.cos(np.pi*2-seg_ori1[DHS-1]), ssc[1]-lyrs1[DHS-1]/2.*np.sin(np.pi*2-seg_ori1[DHS-1]) ] ]
#            ns = len(segs)
#            for i in range(len(segs)):
#                Ptem = [Pnts[i][0]+seg_length[ns-1-i]*np.sin(seg_ori[ns-1-i]), Pnts[i][1]+seg_length[ns-1-i]*np.cos(seg_ori[ns-1-i])]
#                Pnts.append(Ptem)
        elif branchID == 4:
            DHS=42
            #ssc = Pnts1[DHS-1]
            ssc = midPnts1[DHS-1]
            Pnts = [ [ ssc[0]+lyrs1[DHS-1]/2.*np.cos(np.pi*2-seg_ori1[DHS-1]), ssc[1]+lyrs1[DHS-1]/2.*np.sin(np.pi*2-seg_ori1[DHS-1]) ] ]
#            ns = len(segs)
#            for i in range(len(segs)):
#                Ptem = [Pnts[i][0]+seg_length[ns-1-i]*np.sin(seg_ori[ns-1-i]), Pnts[i][1]+seg_length[ns-1-i]*np.cos(seg_ori[ns-1-i])]
#                Pnts.append(Ptem)
        elif branchID == 5:
            DHS=43
            #ssc = Pnts1[DHS-1]
            ssc = midPnts1[DHS-1]
            Pnts = [ [ ssc[0]-lyrs1[DHS-1]/2.*np.cos(np.pi*2-seg_ori1[DHS-1]), ssc[1]-lyrs1[DHS-1]/2.*np.sin(np.pi*2-seg_ori1[DHS-1]) ] ]
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
    
    def SegLyr(self, Pnts, branchID=1):
        """
        Input a list of Points of the specified branch
        Ouput a list of Points for each segments
        """
        segs, seg_length, seg_ori, lyrs = self.SegInfo(branchID)
        ## delete boundary cells first
        segs = segs[1:-1]
        seg_length = seg_length[1:-1]
        seg_ori = seg_ori[1:-1]
        lyrs = lyrs[1:-1]
        
        if branchID in [2,3,4,5]:
            Pnts = Pnts[::-1]
        
        ## calculate the mid point at the center of a segment
        midPnts = []
        for i in range(len(Pnts)-1):
            tem = (Pnts[i]+Pnts[i+1])/2.
            midPnts.append(tem.tolist())
        
        midPnts = np.asarray(midPnts)
        eastPnts = np.zeros_like(midPnts)
        westPnts = np.zeros_like(midPnts)
        for i in range(len(midPnts)):
            eastPnts[i,0] = midPnts[i,0] + lyrs[i]/2.*np.cos(np.pi*2-seg_ori[i])
            eastPnts[i,1] = midPnts[i,1] + lyrs[i]/2.*np.sin(np.pi*2-seg_ori[i])
            westPnts[i,0] = midPnts[i,0] - lyrs[i]/2.*np.cos(np.pi*2-seg_ori[i])
            westPnts[i,1] = midPnts[i,1] - lyrs[i]/2.*np.sin(np.pi*2-seg_ori[i])
        
            #pdb.set_trace()
        return eastPnts, westPnts
    
    
    def Contour2Shp(self, lines, latlon, outfile):
        """
        Converts a matplotlib contour object to a shapefile
        
        The shapefile has the filed Contour corresponding to the contour levels set
        in matplotlib.pyplot.contour
        
        **The projection is hardwired to WGS84 for now. This needs updating.
        """
        # Convert the contour values to a shapefile
        # see this example:
        #    http://invisibleroads.com/tutorials/gdal-shapefile-points-save.html
        
        import osgeo.ogr
        
        # Create the projection
        spatialReference = osgeo.osr.SpatialReference()
        spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
        
        # Create the shape file
        driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
        shapeData = driver.CreateDataSource(outfile)
        
        # Create the layer
        layer = shapeData.CreateLayer('Contour', spatialReference, osgeo.ogr.wkbLineString)
        layerDefinition = layer.GetLayerDefn()
        
        
        # Create fields containing segment infos
        field_def = osgeo.ogr.FieldDefn('branchID', osgeo.ogr.OFTInteger)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('SegID', osgeo.ogr.OFTInteger)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('Lon_west', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('Lat_west', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('Lon_east', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('Lat_east', osgeo.ogr.OFTReal)
        layer.CreateField(field_def)
          
            
#        ctr=0
#        for i in range(len(lines[0])):
#            ctr+=1
#            line = osgeo.ogr.Geometry(osgeo.ogr.wkbLineString)
#            # Add points individually to the line
#            xy = lines[0][i]
#    
#            line.AddPoint_2D(xy[0][0],xy[0][1])
#            line.AddPoint_2D(xy[1][0],xy[1][1])
#            # Update the feature with the line data
#            featureIndex = ctr
#            feature = osgeo.ogr.Feature(layerDefinition)
#            feature.SetGeometry(line)
#            feature.SetFID(featureIndex)
#            feature.SetGeometryDirectly(line)
#            
#            # Set the attribute table
#            feature.SetField('SegID', int(self.segs1[i]))   # convert to int() is necessary, osgeo cannot recognize numpy int32 type
#            feature.SetField('Lon_west', "{:.3f}".format(latlon['westlon1'][i]))
#            feature.SetField('Lat_west', "{:.3f}".format(latlon['westlat1'][i]))
#            feature.SetField('Lon_east', "{:.3f}".format(latlon['eastlon1'][i]))
#            feature.SetField('Lat_east', "{:.3f}".format(latlon['eastlat1'][i]))
#            #feature.SetStyleString("PEN(c:#FF0000,w:5px)")   
#            layer.CreateFeature(feature)
#            
        
        def add_feature(branchID, layer, lines, segs, westlon, westlat, eastlon, eastlat):
            """
            function that adds feature to layer
            """
            
            ctr=0
            for i in range(len(lines)):
                ctr+=1
                line = osgeo.ogr.Geometry(osgeo.ogr.wkbLineString)
                # Add points individually to the line
                xy = lines[i]
        
                line.AddPoint_2D(xy[0][0],xy[0][1])
                line.AddPoint_2D(xy[1][0],xy[1][1])
                # Update the feature with the line data
                featureIndex = ctr
                feature = osgeo.ogr.Feature(layerDefinition)
                feature.SetGeometry(line)
                feature.SetFID(featureIndex)
                feature.SetGeometryDirectly(line)
                
                # Set the attribute table
                feature.SetField('branchID', int(branchID)) 
                feature.SetField('SegID', int(segs[i]))   # convert to int() is necessary, osgeo cannot recognize numpy int32 type
                feature.SetField('Lon_west', "{:.3f}".format(westlon[i]))
                feature.SetField('Lat_west', "{:.3f}".format(westlat[i]))
                feature.SetField('Lon_east', "{:.3f}".format(eastlon[i]))
                feature.SetField('Lat_east', "{:.3f}".format(eastlat[i]))
                #feature.SetStyleString("PEN(c:#FF0000,w:5px)")   
                layer.CreateFeature(feature)
            
            
        add_feature(1, layer, lines[0], self.segs1, latlon['westlon1'], latlon['westlat1'], latlon['eastlon1'], latlon['eastlat1'])
        add_feature(2, layer, lines[1], self.segs2[::-1], latlon['westlon2'], latlon['westlat2'], latlon['eastlon2'], latlon['eastlat2'])
        add_feature(3, layer, lines[2], self.segs3[::-1], latlon['westlon3'], latlon['westlat3'], latlon['eastlon3'], latlon['eastlat3'])
        add_feature(4, layer, lines[3], self.segs4[::-1], latlon['westlon4'], latlon['westlat4'], latlon['eastlon4'], latlon['eastlat4'])
        add_feature(5, layer, lines[4], self.segs5[::-1], latlon['westlon5'], latlon['westlat5'], latlon['eastlon5'], latlon['eastlat5'])
        
            
            
        # Close the shape file
        shapeData.Destroy()
        print ('Complete - shapefile written to:\n      %s'%outfile)
        
        
    """
    Todo:
        A shapefile of contour lines can only be created from rasters 
        either in ArcGIS or using function gdal.ContourGenerate()
        see: https://stackoverflow.com/questions/22100453/gdal-python-creating-contourlines
        and  https://svn.osgeo.org/gdal/trunk/autotest/alg/contour.py
        also http://geoexamples.blogspot.com/2013/08/creating-vectorial-isobands-with-python.html
        So, we either create simple line objects with the attribute table showing the travel time and other info
        or create a raster file from the segment lines and then create the contour based on the raster file
    """
        
        
    
    
#### For testing ####       
        
if __name__ == "__main__":
    
    #filename = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\particle_tracking_test\20191111_baseline4\Bth_WB1.npt'
    filename = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\20191007_SWAT_W2_practice\5. Running CE-QUAL-W2 model\CEQUALW2\Bth_WB1.npt'
    WB = W2_Segmentation(filename)
    #WB.VisSeg()
    #WB.VisSeg_shp()
    WB.VisSeg_shp_GDAL()