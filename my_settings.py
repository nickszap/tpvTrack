#module for storing the files/parameters/... for the tracking code.
#Is it necessary to reload the import to re-call all functions?

import glob
import os, errno
import numpy as np

rEarth = 6370.e3 #radius of spherical Earth (m)
dFilter = 300.e3 #radius for whether nbr extremum is regional extremum
areaOverlap = .1 #fraction of tpv area overlap for determining correspondence

latThresh = 10.*np.pi/180. #segment N of this latitude
trackMinMaxBoth = 0 #0-min, 1-max, 2-both
info = '10N_wrf_trop'

#fDirData = '/data02/cases/summer2006/eraI/pv/'
fDirData = '/data01/tracks/wrf/data/'
filesData = sorted(glob.glob(fDirData+'wrfout_trop*'), key=os.path.getmtime)
print filesData

deltaT = 6.*60.*60. #timestep (s)
#select time intervals within filesData[iFile]...end[-1] means use all times
iTimeStart_fData = [0]
iTimeEnd_fData = [-1]
if (True): #a quick check of specified times
  nFiles = len(filesData)
  if (len(iTimeStart_fData) != nFiles or len(iTimeEnd_fData) != nFiles):
    print "Uhoh, wrong iTime*_data settings in my_settings.py"
    import sys
    sys.exit()

#fDirSave = '/data02/cases/test_segment/testUnified/summer2006/'
fDirSave = '/data01/tracks/wrf/algo/'
if not os.path.exists(fDirSave):
    os.makedirs(fDirSave)

fMesh = filesData[0]  
fMetr = fDirSave+'fields_debug.nc'
fSeg = fDirSave+'seg_debug.nc'
fCorr = fDirSave+'correspond_debug.pkl'
fTrack = fDirSave+'tracks_debug.txt'
fMetrics = fDirSave+'metrics_debug.nc'

inputType = 'wrf_trop'
doPreProc = True
doSeg = True
doMetrics = True
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
