#module for storing the files/parameters/... for the tracking code.
#Is it necessary to reload the import to re-call all functions?

import glob
import os, errno
import numpy as np
import datetime as dt

if False:
  from mpi4py import MPI
  commWorld = MPI.COMM_WORLD
  myRank = commWorld.Get_rank()
  nRanks = commWorld.size
else:
  commWorld = []
  myRank = 0
  nRanks=1
  
def getLimits_startStop(iStartGlobal, iEndGlobal, iWork=myRank, nWork=nRanks):
  #assign contiguous chunks in a sequence to processors. just leave the leftovers to the last processor(s).
  #when the length isn't divisible by the number of workers, this isn't the best solution but we can optimize for that later.
  
  szChunk = int(np.ceil( (iEndGlobal-iStartGlobal)/float(nWork) )) #interval must cover length s.t. # elements for all but last worker
  if (szChunk<1):
    print 'Check logic in getLimits_startStop for your strange case w/ more workers than elements'
  
  iStart = iStartGlobal+iWork*szChunk
  iEnd = iStart+szChunk-1
  iEnd = min(iEndGlobal,iEnd)
  
  return (iStart,iEnd)

rEarth = 6370.e3 #radius of spherical Earth (m)
dFilter = 300.e3 #radius for whether local extremum is regional extremum
areaOverlap = .01 #fraction of tpv area overlap for candidate correspondence
segRestrictPerc = 10. #percentile of boundary amplitudes to restrict watershed basins [0,100]

latThresh = 30.*np.pi/180. #segment N of this latitude
trackMinMaxBoth = 0 #0-min, 1-max (2-both shouldn't be used w/o further development)
info = '30N'

fDirData = '/glade2/scratch2/szapiro/cesm_le/'
if False:
  filesData = [fDirData+'b.e11.B20TRC5CNBDRD.f09_g16.007.cam.h2.TROP_P.1990010100Z-2005123118Z.nc', fDirData+'b.e11.B20TRC5CNBDRD.f09_g16.007.cam.h2.TROP_T.1990010100Z-2005123118Z.nc', fDirData+'b.e11.B20TRC5CNBDRD.f09_g16.007.cam.h2.U.1990010100Z-2005123118Z.nc', fDirData+'b.e11.B20TRC5CNBDRD.f09_g16.007.cam.h2.V.1990010100Z-2005123118Z.nc']
elif False:
  filesData = [fDirData+'b.e11.BRCP85C5CNBDRD.f09_g16.007.cam.h2.TROP_P.2071010100Z-2080123118Z.nc', fDirData+'b.e11.BRCP85C5CNBDRD.f09_g16.007.cam.h2.TROP_T.2071010100Z-2080123118Z.nc', fDirData+'b.e11.BRCP85C5CNBDRD.f09_g16.007.cam.h2.U.2071010100Z-2080123118Z.nc', fDirData+'b.e11.BRCP85C5CNBDRD.f09_g16.007.cam.h2.V.2071010100Z-2080123118Z.nc']
elif True:
  filesData = [fDirData+'b.e11.B20TRC5CNBDRD.f09_g16.020.cam.h2.TROP_P.1990010100Z-2005123118Z.nc', fDirData+'b.e11.B20TRC5CNBDRD.f09_g16.020.cam.h2.TROP_T.1990010100Z-2005123118Z.nc', fDirData+'b.e11.B20TRC5CNBDRD.f09_g16.020.cam.h2.U.1990010100Z-2005123118Z.nc', fDirData+'b.e11.B20TRC5CNBDRD.f09_g16.020.cam.h2.V.1990010100Z-2005123118Z.nc']
elif False:
  filesData = [fDirData+'b.e11.BRCP85C5CNBDRD.f09_g16.020.cam.h2.TROP_P.2071010100Z-2080123118Z.nc', fDirData+'b.e11.BRCP85C5CNBDRD.f09_g16.020.cam.h2.TROP_T.2071010100Z-2080123118Z.nc', fDirData+'b.e11.BRCP85C5CNBDRD.f09_g16.020.cam.h2.U.2071010100Z-2080123118Z.nc', fDirData+'b.e11.BRCP85C5CNBDRD.f09_g16.020.cam.h2.V.2071010100Z-2080123118Z.nc']
print filesData
fileMap = fDirData+'wrfout_mapProj.nc' #for inputType=wrf_trop

#time information of input data
deltaT = 6.*60.*60. #timestep between file times (s)
timeStart = dt.datetime(1990,1,1,0) #time=timeStart+iTime*deltaT
timeDelta = dt.timedelta(seconds=deltaT)
#select time intervals within filesData[iFile]...end[-1] means use all times
if True:
  iTimeStart_fData = [0]*4
  iTimeEnd_fData = [-1]*4 #[4*7]*4
else:
  iWork=3; nWork=4; iStart, iStop = getLimits_startStop(0, 23360, iWork=iWork, nWork=nWork)
  iTimeStart_fData = [iStart]*4
  iTimeEnd_fData = [iStop]*4
  
if (True): #a quick check of specified times
  nFiles = len(filesData)
  if (len(iTimeStart_fData) != nFiles or len(iTimeEnd_fData) != nFiles):
    print "Uhoh, wrong iTime*_data settings in my_settings.py"
    import sys
    sys.exit()

#fDirSave = '/data01/tracks/summer07/tpvTrack/'
fDirSave = fDirData+'tpvTracks/020.1990-2005/'
#fDirSave = fDirData+'tpvTracks/007.2071-2080/'
#fDirSave = '/data01/tracks/wrf/algo/'
if not os.path.exists(fDirSave):
    os.makedirs(fDirSave)

fMesh = filesData[0]  
#fMetr = fDirSave+'fields_{0}.nc'.format(iWork) #fDirSave+'fields.nc'
fMetr = fDirSave+'fields.nc'
fSegFmt = fDirSave+'seg_{0}.nc'
fSeg = fSegFmt.format(myRank)
fSegFinal = fDirSave+'seg.nc'; #fSeg = fSegFinal #for after running seg in parallel...
fCorr = fDirSave+'correspond_horizPlusVert.nc'
fTrackFmt = fDirSave+'tracks_{0}.nc'
fTrack = fTrackFmt.format(myRank)
fTrackFinal = fDirSave+'tracks_low.nc'
fMetrics = fDirSave+'metrics.nc'

inputType = 'cesmLE'
doPreProc = True
doSeg = False
doMetrics = False
doCorr = False
doTracks = False

def silentremove(filename):
  #from http://stackoverflow.com/questions/10840533/most-pythonic-way-to-delete-a-file-which-may-not-exist
  print "Removing file (if it exists): ",filename
  try:
      os.remove(filename)
  except OSError as e: # this would be "except OSError, e:" before Python 2.6
      if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
          raise # re-raise exception if a different error occured
