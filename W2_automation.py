# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 12:17:27 2019

@author: dfeng
"""


import os
import sys
import subprocess 
from shutil import copy2

import pdb


def copy_file(src_file, dst_dir):
    """
    This function is used to move file between folders
    """
    if os.path.isfile(src_file):
        if os.path.isdir(dst_dir):
            pass
        else:
            os.makedirs(dst_dir)
        #print src_file
        srcname = src_file 
        filename = os.path.basename(src_file)
        dstname = os.path.join(dst_dir, filename)
        
        if os.path.isfile(srcname):
            copy2(srcname, dstname) 
            #print srcname,dstname,'success'
        elif os.path.isdir(dstname):
            os.remove(dstname)
            #print 'remove %s' % dstname
            copy2(srcname, dstname)
    if os.path.isfile(src_file):
        os.remove(srcname)
#    if __name__ == '__main__':
#        if len(sys.argv) != 3:
#            print 'need srcFile and dstDir'
#            sys.exit(-1)
#        srcFile = sys.argv[1]
#        dstDir = sys.argv[2]
#        copy_file(srcFile, dstDir)

def move_file(src_file, dst_dir):
    """
    This function is used to move file between folders without deleting the old one
    """
    if os.path.isfile(src_file):
        if os.path.isdir(dst_dir):
            pass
        else:
            os.makedirs(dst_dir)
        #print src_file
        srcname = src_file 
        filename = os.path.basename(src_file)
        dstname = os.path.join(dst_dir, filename)
        
        if os.path.isfile(srcname):
            copy2(srcname, dstname) 
            #print srcname,dstname,'success'
        elif os.path.isdir(dstname):
            os.remove(dstname)
            #print 'remove %s' % dstname
            copy2(srcname, dstname)
    if __name__ == '__main__':
        if len(sys.argv) != 3:
            print 'need srcFile and dstDir'
            sys.exit(-1)
        srcFile = sys.argv[1]
        dstDir = sys.argv[2]
        copy_file(srcFile, dstDir)




class W2_Operation(object):
    """
    general class for system operation of the W2 model
    """
    
    def __init__(self, workdir, **kwargs):
        
        self.__dict__.update(kwargs)
        
        self.workdir = workdir
        
    
    def multi_run(self, timelist):
        """
        run multiple scenarios provided with a list of end times
        """
        
        
        for i in range(len(timelist)):
            endtime = timelist[i]
            print 'Running case with endtime JDAY = %s ...\n'%str(endtime)
            self.UpdateCon(endtime)
            self.model_run()
            self.SaveFinalParticle(endtime)
            
    
    
    def SaveFinalParticle(self, intime):
        """
        save the finalparticle.csv for each simulation
        """
        
        dst_dir = '%s\\finalparticle'%self.workdir  ## directory to save file
        if os.path.isdir(dst_dir):
            pass
        else:
            os.makedirs(dst_dir)
            
        newfile = '%s\\finalparticle_%s.csv'%(self.workdir, str(intime))
        if os.path.isfile(newfile):
            os.remove(newfile)
        os.rename('%s\\finalparticle.csv'%self.workdir, newfile)
        
        copy_file(newfile, dst_dir)
        
        
    
        
    def UpdateCon(self, intime):
        """
        update configuration file
        """
        f=open('%s\\W2_con_template.npt'%self.workdir,'r')
        content=f.readlines()
        open('%s\\W2_con.npt'%self.workdir,'w').write('')
        h=open('%s\\W2_con.npt'%self.workdir,'w')
        
        if intime<10:
            intime2 = '   %s'%str(intime)
        elif intime>=10 and intime<100:
            intime2 = '  %s'%str(intime)
        elif intime>=1000:
            intime2 = ' %s'%str(intime)
        
        for line in content:
            h.write(line.replace('AUTO',intime2))
        h.close
        
        
    def model_run(self):
        """
        run the executable 
        """
        
        os.chdir(self.workdir)
        subprocess.call('w2_v4_64.exe', shell=True)
        
        


#### For testing ####       
        
if __name__ == "__main__":
    
    wdir = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\particle_tracking_test\20191125_1229_test5'
    WC = W2_Operation(wdir)
    #timelist = [11,12,13,14,15,16,17]
    timelist = [18,19,20,21,22,23,24]
    WC.multi_run(timelist)
    