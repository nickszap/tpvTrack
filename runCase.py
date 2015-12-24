import os
import datetime

#mpiexec -n 5 /home/khalbert/anaconda/bin/python /raid1/nick/python_mpi.py

def demo_ensemble_singleCase():
  #track the 50 members in parallel.
  #edit my_settings by hand to run
  
  #anlStart = '2006-11-10-06'
  #tpvStart = '2006-11-15-00'
  
  cmd = 'mpiexec -n 50 -f ~/machines5 /home/khalbert/anaconda/bin/python /raid1/nick/tpvTracks/tigge/tpvTrack/driver.py'
  print cmd
  os.system(cmd)

def demo_ensemble():
  #write a new my_settings.py for this case
  
  fmySettings = '/raid1/nick/tpvTracks/tigge/tpvTrack/my_settings.py' #'test.py'
  
  #from a = sorted(glob.glob('/data02/cases/tigge/200*')); [i+'/' for i in a]
  #anlDirs = ['/raid1/nick/tpvTracks/tigge/2007-04-15-18/', '/raid1/nick/tpvTracks/tigge/2007-05-01-12/', '/raid1/nick/tpvTracks/tigge/2007-05-21-00/', '/raid1/nick/tpvTracks/tigge/2007-06-13-18/', '/raid1/nick/tpvTracks/tigge/2007-07-13-18/', '/raid1/nick/tpvTracks/tigge/2007-08-30-18/', '/raid1/nick/tpvTracks/tigge/2007-12-31-06/','/raid1/nick/tpvTracks/tigge/2008-04-27-12/', '/raid1/nick/tpvTracks/tigge/2008-05-20-00/', '/raid1/nick/tpvTracks/tigge/2008-06-08-18/', '/raid1/nick/tpvTracks/tigge/2008-07-22-06/', '/raid1/nick/tpvTracks/tigge/2008-07-23-00/', '/raid1/nick/tpvTracks/tigge/2008-08-17-18/', '/raid1/nick/tpvTracks/tigge/2008-12-02-18/']
  #anlDirs = ['/raid1/nick/tpvTracks/tigge/2006-11-13-12/', '/raid1/nick/tpvTracks/tigge/2006-11-20-00/', '/raid1/nick/tpvTracks/tigge/2006-12-22-06/', '/raid1/nick/tpvTracks/tigge/2006-12-27-18/']
  #anlDirs = ['/raid1/nick/tpvTracks/tigge/2007-03-14-12/']
  #anlDirs = ['/raid1/nick/tpvTracks/tigge/2009-01-29-06/', '/raid1/nick/tpvTracks/tigge/2009-03-24-12/', '/raid1/nick/tpvTracks/tigge/2009-04-28-06/', '/raid1/nick/tpvTracks/tigge/2009-05-09-18/', '/raid1/nick/tpvTracks/tigge/2009-06-17-18/', '/raid1/nick/tpvTracks/tigge/2009-08-05-18/', '/raid1/nick/tpvTracks/tigge/2009-08-29-12/', '/raid1/nick/tpvTracks/tigge/2009-09-25-00/', '/raid1/nick/tpvTracks/tigge/2009-10-19-06/', '/raid1/nick/tpvTracks/tigge/2009-11-18-18/', '/raid1/nick/tpvTracks/tigge/2010-01-31-18/', '/raid1/nick/tpvTracks/tigge/2010-04-01-06/', '/raid1/nick/tpvTracks/tigge/2010-04-26-06/', '/raid1/nick/tpvTracks/tigge/2010-05-07-06/', '/raid1/nick/tpvTracks/tigge/2010-05-26-18/', '/raid1/nick/tpvTracks/tigge/2010-06-25-18/', '/raid1/nick/tpvTracks/tigge/2010-08-03-18/', '/raid1/nick/tpvTracks/tigge/2010-08-19-18/', '/raid1/nick/tpvTracks/tigge/2011-01-13-06/', '/raid1/nick/tpvTracks/tigge/2011-03-21-18/', '/raid1/nick/tpvTracks/tigge/2011-03-27-18/', '/raid1/nick/tpvTracks/tigge/2011-05-06-18/', '/raid1/nick/tpvTracks/tigge/2011-06-04-18/', '/raid1/nick/tpvTracks/tigge/2011-09-19-18/', '/raid1/nick/tpvTracks/tigge/2011-09-20-06/', '/raid1/nick/tpvTracks/tigge/2011-12-20-18/', '/raid1/nick/tpvTracks/tigge/2012-01-20-18/', '/raid1/nick/tpvTracks/tigge/2012-03-06-00/', '/raid1/nick/tpvTracks/tigge/2012-05-05-18/', '/raid1/nick/tpvTracks/tigge/2012-06-11-06/', '/raid1/nick/tpvTracks/tigge/2012-06-18-06/', '/raid1/nick/tpvTracks/tigge/2012-08-06-18/', '/raid1/nick/tpvTracks/tigge/2013-02-01-00/', '/raid1/nick/tpvTracks/tigge/2013-05-18-00/', '/raid1/nick/tpvTracks/tigge/2013-06-30-00/', '/raid1/nick/tpvTracks/tigge/2013-09-14-06/', '/raid1/nick/tpvTracks/tigge/2013-12-24-18/', '/raid1/nick/tpvTracks/tigge/2014-04-20-18/', '/raid1/nick/tpvTracks/tigge/2014-07-12-00/', '/raid1/nick/tpvTracks/tigge/2014-08-05-12/', '/raid1/nick/tpvTracks/tigge/2014-09-22-18/', '/raid1/nick/tpvTracks/tigge/2014-11-09-06/', '/raid1/nick/tpvTracks/tigge/2015-03-18-06/', '/raid1/nick/tpvTracks/tigge/2015-06-17-06/', '/raid1/nick/tpvTracks/tigge/2015-06-19-12/', '/raid1/nick/tpvTracks/tigge/2015-06-27-18/']
  anlDirs = ['/raid1/nick/tpvTracks/tigge/2006-12-22-06/']
  
  # use from getData.py and then change Hours = 0. could just fix hours in time when writing my_settings.py
  #trackDates = [datetime.datetime(2007, 4, 20, 0, 0), datetime.datetime(2007, 5, 6, 0, 0), datetime.datetime(2007, 5, 26, 0, 0), datetime.datetime(2007, 6, 18, 0, 0), datetime.datetime(2007, 7, 18, 18, 0), datetime.datetime(2007, 9, 4, 18, 0), datetime.datetime(2008, 1, 5, 6, 0), datetime.datetime(2008, 5, 2, 12, 0), datetime.datetime(2008, 5, 25, 0, 0), datetime.datetime(2008, 6, 13, 0, 0), datetime.datetime(2008, 7, 27, 0, 0), datetime.datetime(2008, 7, 28, 0, 0), datetime.datetime(2008, 8, 22, 0, 0), datetime.datetime(2008, 12, 7, 0, 0)]
  #trackDates = [datetime.datetime(2006, 11, 18, 0, 0), datetime.datetime(2006, 11, 25, 0, 0), datetime.datetime(2006, 12, 27, 0, 0), datetime.datetime(2007, 1, 1, 0, 0)]
  #trackDates = [datetime.datetime(2007, 3, 19, 0, 0)]
  #trackDates = [datetime.datetime(2009, 2, 3, 0, 0), datetime.datetime(2009, 3, 29, 0, 0), datetime.datetime(2009, 5, 3, 0, 0), datetime.datetime(2009, 5, 14, 0, 0), datetime.datetime(2009, 6, 22, 0, 0), datetime.datetime(2009, 8, 10, 0, 0), datetime.datetime(2009, 9, 3, 0, 0), datetime.datetime(2009, 9, 30, 0, 0), datetime.datetime(2009, 10, 24, 0, 0), datetime.datetime(2009, 11, 23, 0, 0), datetime.datetime(2010, 2, 5, 0, 0), datetime.datetime(2010, 4, 6, 0, 0), datetime.datetime(2010, 5, 1, 0, 0), datetime.datetime(2010, 5, 12, 0, 0), datetime.datetime(2010, 5, 31, 0, 0), datetime.datetime(2010, 6, 30, 0, 0), datetime.datetime(2010, 8, 8, 0, 0), datetime.datetime(2010, 8, 24, 0, 0), datetime.datetime(2011, 1, 18, 0, 0), datetime.datetime(2011, 3, 26, 0, 0), datetime.datetime(2011, 4, 1, 0, 0), datetime.datetime(2011, 5, 11, 0, 0), datetime.datetime(2011, 6, 9, 0, 0), datetime.datetime(2011, 9, 24, 0, 0), datetime.datetime(2011, 9, 25, 0, 0), datetime.datetime(2011, 12, 25, 0, 0), datetime.datetime(2012, 1, 25, 0, 0), datetime.datetime(2012, 3, 11, 0, 0), datetime.datetime(2012, 5, 10, 0, 0), datetime.datetime(2012, 6, 16, 0, 0), datetime.datetime(2012, 6, 23, 0, 0), datetime.datetime(2012, 8, 11, 0, 0), datetime.datetime(2013, 2, 6, 0, 0), datetime.datetime(2013, 5, 23, 0, 0), datetime.datetime(2013, 7, 5, 0, 0), datetime.datetime(2013, 9, 19, 0, 0), datetime.datetime(2013, 12, 29, 0, 0), datetime.datetime(2014, 4, 25, 0, 0), datetime.datetime(2014, 7, 17, 0, 0), datetime.datetime(2014, 8, 10, 0, 0), datetime.datetime(2014, 9, 27, 0, 0), datetime.datetime(2014, 11, 14, 0, 0), datetime.datetime(2015, 3, 23, 0, 0), datetime.datetime(2015, 6, 22, 0, 0), datetime.datetime(2015, 6, 24, 0, 0), datetime.datetime(2015, 7, 2, 0, 0)]
  trackDates = [datetime.datetime(2006, 12, 27, 0, 0)]
  
  nCases = len(anlDirs)
  for iCase in xrange(nCases):
    fDir = anlDirs[iCase]
    t0 = trackDates[iCase]
    
    #write my_settings
    write_mySettings(fmySettings, fDir, t0)
    
    #cmd = 'mpiexec -n 50 -f ~/machines5 /home/khalbert/anaconda/bin/python /raid1/nick/tpvTracks/tigge/tpvTrack/driver.py {0} {1}'.format(fDir, timeString) #command line argument pass
    cmd = 'mpiexec -n 50 -f ~/machines5 /home/khalbert/anaconda/bin/python /raid1/nick/tpvTracks/tigge/tpvTrack/driver.py'
    os.system(cmd)

def write_mySettings(fNameOut, fDir, t0):

  s = '''
#module for storing the files/parameters/... for the tracking code.
#Is it necessary to reload the import to re-call all functions?

import glob
import os, errno
import numpy as np
import datetime as dt

from mpi4py import MPI
commWorld = MPI.COMM_WORLD
myRank = commWorld.Get_rank()
nRanks = commWorld.size

rEarth = 6370.e3 #radius of spherical Earth (m)
dFilter = 300.e3 #radius for whether local extremum is regional extremum
areaOverlap = .05 #fraction of tpv area overlap for candidate correspondence

latThresh = 30.*np.pi/180. #segment N of this latitude
trackMinMaxBoth = 0 #0-min, 1-max (2-both shouldn't be used w/o further development)
info = '30N_tigge'

fDirData = '{0}'
filesData = sorted(glob.glob(fDirData+'*.nc'), key=os.path.getmtime)
#fDirData = '/data01/tracks/summer07/eraI/'
#filesData = sorted(glob.glob(fDirData+'ERAI*.nc'), key=os.path.getmtime)
print filesData
fileMap = fDirData+'wrfout_mapProj.nc' #for inputType=wrf_trop

#time information of input data
deltaT = 6.*60.*60. #timestep between file times (s)
timeStart = dt.datetime({1},{2},{3},{4}) #time=timeStart+iTime*deltaT
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

#fDirSave = '/data01/tracks/summer07/tpvTrack/'
fDirSaveTemplate = fDirData+'{{0}}/'
fDirSave = fDirSaveTemplate
#fDirSave = '/data01/tracks/wrf/algo/'

fMesh = filesData[0]  
fMetr = fDirSave+'fields.nc'
fSeg = fDirSave+'seg.nc'
fCorr = fDirSave+'correspond_horizPlusVert.nc'
fTrack = fDirSave+'tracks_low_horizPlusVert.nc'
fMetrics = fDirSave+'metrics.nc'

def setFilenames(my_settings, iMem):
  my_settings.fDirSave = my_settings.fDirSaveTemplate.format(iMem)
  #fDirSave = '/data01/tracks/wrf/algo/'
  if not os.path.exists(fDirSave):
      os.makedirs(fDirSave)
      
  my_settings.fMetr = my_settings.fDirSave+'fields.nc'
  my_settings.fSeg = my_settings.fDirSave+'seg.nc'
  my_settings.fCorr = my_settings.fDirSave+'correspond_horizPlusVert.nc'
  my_settings.fTrack = my_settings.fDirSave+'tracks_low_horizPlusVert.nc'
  my_settings.fMetrics = my_settings.fDirSave+'metrics.nc'
  
  return (my_settings.fDirSave, my_settings.fMetr, my_settings.fSeg, my_settings.fCorr, my_settings.fTrack, my_settings.fMetrics)

inputType = 'eraI'
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

'''.format(fDir,t0.year,t0.month,t0.day,t0.hour)
  
  print 'Writing my_settings.py to: ', fNameOut
  f = open(fNameOut,'w')
  f.write(s)
  f.close()

if __name__=='__main__':
  demo_ensemble()


