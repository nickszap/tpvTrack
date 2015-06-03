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
trackMinMaxBoth = 0 #0-min, 1-max (2-both shouldn't be used w/o further development)
info = '30N'

fDirData = '/data01/tracks/cesmLE/mem.024.2026/'
filesData = sorted(glob.glob(fDirData+'trop_*.nc'), key=os.path.getmtime)
#fDirData = '/data01/tracks/summer07/eraI/'
#filesData = sorted(glob.glob(fDirData+'ERAI*.nc'), key=os.path.getmtime)
print filesData
fileMap = fDirData+'wrfout_mapProj.nc' #for inputType=wrf_trop

#time information of input data
deltaT = 6.*60.*60. #timestep between file times (s)
t0 = 2034; tRef = 2026; iTime = 365*4*(t0-tRef)
timeStart = dt.datetime(t0,1,1,0) #time=timeStart+iTime*deltaT
timeDelta = dt.timedelta(seconds=deltaT)

iTimeStart_fData = [iTime]
iTimeEnd_fData = [iTime+365*4] #[-1]
if (True): #a quick check of specified times
  nFiles = len(filesData)
  if (len(iTimeStart_fData) != nFiles or len(iTimeEnd_fData) != nFiles):
    print "Uhoh, wrong iTime*_data settings in my_settings.py"
    import sys
    sys.exit()

#fDirSave = '/data01/tracks/cesmLE/mem.024.2026/'
fDirSave = fDirData
if not os.path.exists(fDirSave):
    os.makedirs(fDirSave)

fMesh = filesData[0]  
fMetr = fDirSave+'fields_{0}.nc'.format(t0)
fSeg = fDirSave+'seg_{0}.nc'.format(t0)
fCorr = fDirSave+'correspond_horizPlusVert_{0}.nc'.format(t0)
fTrack = fDirSave+'tracks_low_horizPlusVert_{0}.nc'.format(t0)
fMetrics = fDirSave+'metrics_{0}.nc'.format(t0)

inputType = 'cesmLE'
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
