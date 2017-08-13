#module for storing the files/parameters/... for the tracking code.
#Is it necessary to reload the import to re-call all functions?

import glob
import os, errno
import numpy as np
import datetime as dt
from mpi4py import MPI
import netCDF4
import netcdftime as nct

commWorld = MPI.COMM_WORLD
myRank = commWorld.Get_rank()
nRanks = commWorld.size
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
areaOverlap = .01 #fraction of tpv area overlap for determining correspondence

latThresh = 30.*np.pi/180. #segment N of this latitude
trackMinMaxBoth = 0 #0-min, 1-max (2-both shouldn't be used w/o further development)
info = '30N_erai'

year = '2015'
fDirData = '/erai/'+year+'/'
filesData = [fDirData+'/g.'+year+'.nc',fDirData+'/u.'+year+'.nc',fDirData+'/v.'+year+'.nc']

#time information of input data
timeStart = dt.datetime(2015,1,1,0,0)
timeEnd = dt.datetime(2015,1,5,18,0)
deltaT = 6.*60.*60. #timestep between file times (s)
timeDelta = dt.timedelta(seconds=deltaT)

nctime = nct.utime('hours since 1800-01-01 00:00:00')
tNumStart = nctime.date2num(timeStart)
tNumEnd = nctime.date2num(timeEnd)

data = netCDF4.Dataset(filesData[0],'r')
time = data.variables['time'][:] 
iTimeStart_fData = np.where(time==tNumStart)[0][0]
iTimeEnd_fData = np.where(time==tNumEnd)[0][0]


fDirSave = '/home/track/Jan_1-5_2015/'

if not os.path.exists(fDirSave):
    os.makedirs(fDirSave)

fMesh = filesData[0]  
fMetr = fDirSave+'fields.nc'
fSegFmt = fDirSave+'seg_{0}.nc'
fSeg = fSegFmt.format(myRank)
fSegFinal = fDirSave+'seg.nc'; #fSeg = fSegFinal #for after running seg in parallel...
fCorr = fDirSave+'correspond_horizPlusVert.nc'
fTrackFmt = fDirSave+'tracks_{0}.nc'
fTrack = fTrackFmt.format(myRank)
fTrackFinal = fDirSave+'tracks_low.nc'
fMetrics = fDirSave+'metrics.nc'

inputType = 'eraI'
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
