#module for storing the files/parameters/... for the tracking code.
#Is it necessary to reload the import to re-call all functions?

import glob
import os, errno
import numpy as np
import datetime as dt

rEarth = 6370.e3 #radius of spherical Earth (m)
dFilter = 300.e3 #radius for whether local extremum is regional extremum
areaOverlap = .1 #fraction of tpv area overlap for determining correspondence

latThresh = 30.*np.pi/180. #segment N of this latitude
trackMinMaxBoth = 0 #0-min, 1-max, 2-both
info = '30N_eraI'

fDirData = '/data02/cases/summer2006/eraI/pv/'
filesData = sorted(glob.glob(fDirData+'eraI_theta-u-v_2pvu_2006-06-01_09-30*'), key=os.path.getmtime)
print filesData
fileMap = fDirData+'wrfout_mapProj.nc' #for inputType=wrf_trop

#time information of input data
deltaT = 6.*60.*60. #timestep between file times (s)
timeStart = dt.datetime(2006,6,1,0) #time=timeStart+iTime*deltaT
timeDelta = dt.timedelta(seconds=deltaT)
#select time intervals within filesData[iFile]...end[-1] means use all times
iTimeStart_fData = [0]
iTimeEnd_fData = [-1]
if (True): #a quick check of specified times
  nFiles = len(filesData)
  if (len(iTimeStart_fData) != nFiles or len(iTimeEnd_fData) != nFiles):
    print "Uhoh, wrong iTime*_data settings in my_settings.py"
    import sys
    sys.exit()

fDirSave = '/data01/tracks/summer06/testOverlap/'
#fDirSave = '/data01/tracks/wrf/algo/'
if not os.path.exists(fDirSave):
    os.makedirs(fDirSave)

fMesh = filesData[0]  
fMetr = fDirSave+'fields_debug.nc'
fSeg = fDirSave+'seg_debug.nc'
fCorr = fDirSave+'correspond_debug.pkl'
fTrack = fDirSave+'tracks_debug.txt'
fMetrics = fDirSave+'metrics_debug.nc'

inputType = 'eraI'
doPreProc = False
doSeg = False
doMetrics = False
doCorr = True
doTracks = True

def silentremove(filename):
  #from http://stackoverflow.com/questions/10840533/most-pythonic-way-to-delete-a-file-which-may-not-exist
  print "Removing file (if it exists): ",filename
  try:
      os.remove(filename)
  except OSError as e: # this would be "except OSError, e:" before Python 2.6
      if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
          raise # re-raise exception if a different error occured
