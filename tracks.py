#The track program creates correspondences between basins at neighboring times.
#We want tracks that extend over the lifetime of each/all basins

import os
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import datetime as dt

import basinMetrics
import correspond

r2d = 180./np.pi

def form_majorTracks_iTimeStart(dataCorr, iTimeStart, iTimeEnd):
  #stitch together tracks of major 1-1 correspondences for tracks starting at iTimeStart, lasting to at most iTimeEnd
  trackList = []
  
  #start tracking only the basin's that weren't the tail of a major correspondence at the previous time
  sitesNotStarting = []
  if (iTimeStart>0): #all tracks start at the initial time
    sites0, sites1 = correspond.get_majorCorrespondences_iTime(dataCorr, iTimeStart-1)
    sitesNotStarting = sites1
  
  trackList_building = []
  sites0, sites1 = correspond.get_majorCorrespondences_iTime(dataCorr, iTimeStart)
  nSites0 = len(sites0)
  for iSite0 in xrange(nSites0):
    if (sites0[iSite0] not in sitesNotStarting):
      trackList_building.append([ sites0[iSite0],sites1[iSite0] ])
  
  #now continue the major tracks from iTimeStart+1 for as long as they last.
  #move from building to final trackList when track stops
  for iTime in xrange(iTimeStart+1,iTimeEnd+1):
    if ( len(trackList_building) == 0 ): #all tracks finished
      break
    sites0, sites1 = correspond.get_majorCorrespondences_iTime(dataCorr, iTime)
    #sites0 = sites0.tolist();
    
    newTrackList_building = []
    for trackSeq in trackList_building:
      siteInd = trackSeq[-1]
      if (siteInd not in sites0): #track finished
        trackList.append(trackSeq)
      else:
        #add continuing time's site to sequence
        iSite = sites0.index(siteInd)
        trackSeq.append(sites1[iSite])
        
        #add sequence to continue on for next time
        newTrackList_building.append(trackSeq)
    #use extended sequences for next time iteration
    trackList_building = newTrackList_building
    
  #copy tracks that lasted entire period
  for t in trackList_building: #does nothing if empty
    trackList.append(t)
    
  return trackList

def run_majorTracks_timeInterval(fNameTracks, fCorr, iTimeStart, iTimeEnd, nTimesGlobal, timeStartGlobal, deltaTGlobal, fMetrics=''):
  #written s.t. major tracks can be stitched together in parallel over time.
  #For inputs,
  # iTimeStart: first time index tracks can start
  # iTimeEnd: last time index tracks can start
  # nTimesGlobal: # times possible in longest track (1 less than # of segmentation times)
  
  
  dataTracks = write_tracks_metrics_netcdf_header(fNameTracks, 'test', nTimesGlobal, nTimesGlobal, timeStartGlobal, deltaTGlobal); iTrackGlobal=0
  dataCorr = netCDF4.Dataset(fCorr, 'r')
  
  for timeInd in xrange(iTimeStart,iTimeEnd+1):
    trackList = form_majorTracks_iTimeStart(dataCorr, timeInd, nTimesGlobal-1); print "Formed tracks for iTimeGlobal: ", timeInd
    #print trackList
    
    #write to file
    if (fMetrics==''):
      write_tracks_cells(fNameTracks, trackList)
    else:
      dataMetrics = netCDF4.Dataset(fMetrics,'r')
      
      #write_tracks_metrics_iTime(fNameTracks, timeInd, trackList, dataMetrics, timeStartGlobal, deltaTGlobal)
      iTrackGlobal = write_tracks_metrics_iTime_netcdf(dataTracks, timeInd, iTrackGlobal, trackList, dataMetrics)
      
      dataMetrics.close()
  
  dataCorr.close()
  dataTracks.close()
  
def write_tracks_cells(fNameTracks, trackList):
  print "Appending to file: "+fNameTracks
  f = open(fNameTracks,'a')
  for track in trackList:
    s = ' '.join(str(i) for i in track)
    f.write(s+'\n')
  f.close()

def read_tracks_cells(fNameTracks):
  #return trackList
  f = open(fNameTracks,'r')
  trackList = []
  for line in f:    
    cellStr = line.strip().split()
    trackSeq = [int(i) for i in cellStr]
    trackList.append(trackSeq)
  return trackList

timeStringFormat = "%Y-%m-%d-%H"
def write_tracks_metrics_netcdf_header(fName, info, nTimesInTrackMax, nTimes, timeStartGlobal, deltaTGlobal):
  '''
  For inputs,
  nTimesInTrackMax: maximum length of track
  nTimes: number of times in tracking interval
  nTimesInTrackMax <= nTimes
  
  Maybe dumping out the tracks to text file as:
  iTimeStartTrack1 nTimesTrack1 timeStart timeEnd
  (time1) metric1 metric2 ... for track1
  (time2) metric1 metric2 ... for track1
  -1
  ends up taking alot of time since we read many small chunks from basinMetrics.
  
  Writing tracks out as metric1[iTime,iTrack] will let us load basinMetrics[iTime], but we'll
  see how costly it is to write to scattered locations within the file. There are probably intermediate buffers and such.
  
  "ragged" or "vlen" arrays still don't make sense to me, so we'll use extra padding instead.
  '''
  data = netCDF4.Dataset(fName, 'w', format='NETCDF4')
  data.description = info
  
  # dimensions
  data.createDimension('nTracks', None)
  #data.createDimension('nTimesTrack', nTimesInTrackMax+1)
  data.createDimension('nTimesTrack', None) #unlimited dimension
  data.createDimension('nTimes', nTimes)
  
  tNow = dt.datetime.now().strftime(timeStringFormat)
  lenTime = len(tNow)
  data.createDimension('lenTimeString', lenTime)
  
  # variables
  data.createVariable('timeStamp', str, ('nTimes',))
  data.createVariable('iTimeStart', 'i4', ('nTracks',))
  data.createVariable('lenTrack', 'i4', ('nTracks',))
  data.createVariable('siteExtr', 'i4', ('nTracks','nTimesTrack',))
  
  for key in basinMetrics.metricKeys:
    data.createVariable(key, 'f8', ('nTracks','nTimesTrack',))
  
  for iTime0 in xrange(nTimes):
    tStart = timeStartGlobal+deltaTGlobal*iTime0; tStart = tStart.strftime(timeStringFormat)
    data.variables['timeStamp'][iTime0] = tStart
  
  return data

def write_tracks_metrics_iTime_netcdf(data, iTime0, iTrackGlobal0, trackList, dataMetrics):
  
  nTracks = len(trackList)
  trackLengths = np.array([len(track) for track in trackList], dtype=int)
  maxLength = np.max(trackLengths)
  #print "Maximum track length={0} at time {1}".format(maxLength, iTrackGlobal)
  print "Maximum track length={0} at time {1}".format(maxLength, data.variables['timeStamp'][iTime0])
  
  data.variables['lenTrack'][iTrackGlobal0:iTrackGlobal0+nTracks] = trackLengths[:]
  data.variables['iTimeStart'][iTrackGlobal0:iTrackGlobal0+nTracks] = iTime0
  
  for iTime in xrange(maxLength):
    iTimeGlobal = iTime0+iTime
    sites = dataMetrics.variables['sites'][iTimeGlobal,:]
    
    for key in basinMetrics.metricKeys:
      vals = dataMetrics.variables[key][iTimeGlobal,:]
    
      for iTrack in xrange(nTracks):
        if (trackLengths[iTrack]-1<iTime): #can only index an array of length 4 with [3]
          continue
        iTrackGlobal = iTrackGlobal0+iTrack
        #print iTrackGlobal, iTime, trackList[iTrack]
        site = trackList[iTrack][iTime]
        iSite = np.where(sites==site)[0][0]
        
        data.variables['siteExtr'][iTrackGlobal,iTime] = site
        data.variables[key][iTrackGlobal, iTime] = vals[iSite]
  
  iTrackGlobal = iTrackGlobal0+nTracks
  return iTrackGlobal
  
def write_tracks_metrics_iTime(fSave, iTime0, trackList, dataMetrics, timeStartGlobal, deltaTGlobal):
  '''
  format is:
  (header) metric1 metric2 ...
  iTimeStartTrack1 nTimesTrack1 timeStart timeEnd
  (time1) metric1 metric2 ... for track1
  (time2) metric1 metric2 ... for track1
  -1
  iTimeStartTrack2 nTimesTrack2 timeStart timeEnd
  (time1) metric1 metric2 ... for track2
  (time2) metric1 metric2 ... for track2
  -1
  .
  .
  .
  '''
  
  f = open(fSave,'a')
  
  #add the header if "appending" actually created the file
  if (f.tell() == 0): #file pointer at 0th byte
    varNames = basinMetrics.metricKeys
    s = ' '.join(varNames);
    f.write(s+'\n'); s = ''
  
  for track in trackList:
    s = ''
    nTimes = len(track);
    tStart = timeStartGlobal+deltaTGlobal*iTime0; tEnd = timeStartGlobal+deltaTGlobal*(iTime0+nTimes-1);
    tStart = tStart.strftime(timeStringFormat)
    tEnd = tEnd.strftime(timeStringFormat)
    s += '{0} {1} {2} {3}\n'.format(iTime0, nTimes, tStart, tEnd)
    for iTime in xrange(nTimes):
      site = track[iTime]
      vals = basinMetrics.get_metrics_basin(dataMetrics, iTime+iTime0, site)
      
      valsStr = '';
      for val in vals:
        valsStr += '{0:g} '.format(val)
      #valsStr = str(vals)[1:-1]
      s += valsStr+'\n'
    #end iTime  
    s += '-1\n'
    f.write(s); s = ''
  #end iTrack
  
  f.close()

def read_tracks_metrics(fNameTracks, metricNames):
  #Input the name of the track file and list of metricNames strings.
  #return list of numpy arrays list[iTrack][iTime,iMetric] with the metric properties of the tracked TPVs.
  
  nMetrics = len(metricNames)
  data = netCDF4.Dataset(fNameTracks,'r')
  nTracks = len(data.dimensions['nTracks'])
  
  trackList = []
  timeStartList = []
  
  for iTrack in xrange(nTracks):
    nTimes = data.variables['lenTrack'][iTrack]
    iTimeStart = data.variables['iTimeStart'][iTrack]
    timeStartList.append(data.variables['timeStamp'][iTimeStart])
    
    trackVals = np.empty((nTimes,nMetrics),dtype=float)
    for iMetric in xrange(nMetrics):
        key = metricNames[iMetric]
        trackVals[:,iMetric] = data.variables[key][iTrack,0:nTimes]
    trackList.append(trackVals)
        
  data.close()
  return trackList, timeStartList
  
def plot_tracks_cells(fTracks, mesh, fDirSave):
  f = open(fTracks,'r')
  
  m = Basemap(projection='ortho',lon_0=0,lat_0=89.5, resolution='l')
  
  plt.figure()
  m.drawcoastlines()
  
  for line in f:    
    cellStr = line.strip().split()
    trackList = [int(i) for i in cellStr]
    if (len(trackList)<3):
      continue
    
    lat, lon = mesh.get_latLon_inds(np.array(trackList,dtype=int))
    lat *= r2d; lon *= r2d
    x,y = m(lon,lat)
    #print lat; print lon
    
    #mark beginning and ending of track
    m.scatter(x[0],y[0], marker='+', color='g', s=45)
    m.scatter(x[-1],y[-1], marker='o', color='r', s=10)
    
    #plot track
    m.plot(x,y, 'b-')
    
  if (False):
    plt.show()
  else:
    fName = 'tracks_debug.png'
    fSave = fDirSave+fName
    print "Picture of tracks from {0}: {1}".format(fTracks,fSave)
    plt.savefig(fSave); plt.close()
      
  f.close()

def plot_tracks_metrics(fTracks, fSave):
  
  metricNames = ['thetaExtr', 'latExtr', 'lonExtr']
  latInd = metricNames.index('latExtr')
  lonInd = metricNames.index('lonExtr')
  varKey = 'thetaExtr'
  varInd = metricNames.index(varKey); varMin = 270.; varMax = 310.; #varMin= 320.; varMax = 380.;

  trackList, timeList = read_tracks_metrics(fTracks, metricNames)
  
  m = Basemap(projection='ortho',lon_0=0,lat_0=89.5, resolution='l')
  
  #ax = plt.figure()
  ax = plt.gca()
  m.drawcoastlines()
  
  for iTrack,track in enumerate(trackList):
    nTimes = track.shape[0]
    if (True):
      if (nTimes<56):
        continue
    
    lat = track[:,latInd]
    lon = track[:,lonInd]
    x,y = m(lon,lat)
    print timeList[iTrack], nTimes; print lat; print lon
    
    #mark beginning and ending of track
    m.scatter(x[0],y[0], marker='+', color='g', s=45)
    m.scatter(x[-1],y[-1], marker='o', color='r', s=10)
    
    #plot track, with color representing value
    #m.plot(x,y, 'b-')
    vals = track[:,varInd]
    colorline(x, y, z=vals, cmap=plt.get_cmap('RdBu_r'), norm=plt.Normalize(varMin, varMax), linewidth=3, alpha=1.0, ax=ax)
  
  #plt.colorbar()
  s = 'TPV {0}, [{1},{2}]'.format(varKey, varMin, varMax)
  plt.title(s)
  if (False):
    plt.show()
  else:
    print "Saving image of tracks from {0}: {1}".format(fTracks,fSave)
    plt.savefig(fSave); plt.close()
    
def demo_plotMetrics(fTracks):

  metricNames = ['thetaExtr', 'latExtr']
  #latInd = metricNames.index('latExtr')

  trackList = read_tracks_metrics(fTracks, metricNames)
  
  for iMetric,metricName in enumerate(metricNames):
  
    plt.figure()
  
    for track in trackList:
      nTimes = track.shape[0]
      if (True):
        if (nTimes<4):
          continue
      #lat = track[:,latInd]
      plt.plot(track[:,iMetric])
    
    plt.title(metricName)
    plt.show()

def demo_compareMetrics(fTracks):
  metricNames = ['rEquiv', 'vortMean']
  #latInd = metricNames.index('latExtr')

  trackList = read_tracks_metrics(fTracks, metricNames)
  
  plt.figure()

  for track in trackList:
    nTimes = track.shape[0]; #print track, track.shape
    if (True):
      if (nTimes<4):
        continue
    #lat = track[:,latInd]
    plt.scatter(track[:,0], track[:,1])
  
  s = '{0} vs {1}'.format(metricNames[0], metricNames[1])
  plt.title(s)
  plt.tight_layout()
  plt.ylim([1.e-6, 2.e-4]); plt.semilogy()
  plt.show()

def demo_plotLifetimes(fTracks):

  metricNames = ['latExtr']
  trackList = read_tracks_metrics(fTracks, metricNames)
  
  plt.figure()
  
  vals = []
  for track in trackList:
    nTimes = track.shape[0]
    vals.append(nTimes)
    if (np.sum(track[:,0]>70.)>7*4):
      print nTimes, track
  
  if (True):
    vals = [i for i in vals if i>6]
  
  plt.hist(vals, cumulative=True, bins=20)
  plt.title('Lifetime (timesteps)')
  plt.show()

#The following 2 fcts are taken from:
# http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    return segments


# Interface to LineCollection:
def colorline(x, y, z=None, cmap=plt.get_cmap('Blues_r'), norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0, ax=plt.gca()):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''
    
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))
           
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
        
    z = np.asarray(z)
    
    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)
    
    #ax = plt.gca()
    ax.add_collection(lc)
    
    return lc
    
def combineParallelFiles(fOut, filesIn, 
                          keys1d = ['iTimeStart','lenTrack'], 
                          keys2d=['siteExtr', 'circ', 'vortMean', 'ampMaxMin', 'rEquiv', 'thetaVol', 'ampMean', 'thetaExtr', 'latExtr', 'lonExtr']):
  #put all the tracks starting in contiguous chunks of times from the separate workers together into 1 file
  #'timeStamp' is correct for the global file in all workers' files.
  
  #rather than remake the header of a netcdf file, append to a copy of one of the existing
  cmd = 'cp {0} {1}'.format(filesIn[0], fOut)
  print cmd; os.system(cmd)
  
  dataOut = netCDF4.Dataset(fOut,'a')
  nFilesIn = len(filesIn)
  iTrackGlobal = 0
  for iFile in xrange(nFilesIn):
    fIn = filesIn[iFile];
    dataIn = netCDF4.Dataset(fIn,'r')
    nTracks = len(dataIn.dimensions['nTracks'])
    for key in keys1d:
      vals = dataIn.variables[key][:]
      dataOut.variables[key][iTrackGlobal:iTrackGlobal+nTracks] = vals[:]
    for key in keys2d:
      vals = dataIn.variables[key][:,:]
      dataOut.variables[key][iTrackGlobal:iTrackGlobal+nTracks,:] = vals[:,:]
    
    dataIn.close()
    iTrackGlobal = iTrackGlobal + nTracks
    
  dataOut.close()



