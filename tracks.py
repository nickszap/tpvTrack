#The track program creates correspondences between basins at neighboring times.
#We want tracks that extend over the lifetime of each/all basins

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

def form_track_site(dataCorr, iTimeStart, iTimeEnd, site0, trackOnlyMajor):
  # follow a given site throughout the correspondences "tree" and split tree into individual tracks
  
  tracks_checkContinue = [[site0]]
  trackList = []
  while (len(tracks_checkContinue)>0):
    #can continue a track if: have more times AND site corresponds to future site
    basinTrack = tracks_checkContinue.pop()
    nTimes = len(basinTrack); site0 = basinTrack[-1]
    
    iTime = iTimeStart+nTimes-1 #0-indexing for time (time0 will have 1 site, time1 will have 2 sites)
    if (iTime>=iTimeEnd): #no more times left
      trackList.append(basinTrack)
      continue
    
    corrSites, corrTypes = correspond.get_correspondingSites(dataCorr, iTime, site0)
    if (len(corrSites)<1): #don't connect to future site
      trackList.append(basinTrack)
      continue
    
    if (trackOnlyMajor):
      if (2 not in corrTypes):
        #don't connect to future site
        trackList.append(basinTrack)
        continue
      
    #if here, site0 connects to >=1 future sites so can continue track
    nSites1 = len(corrSites)
    for iSite1 in xrange(nSites1):
      if (trackOnlyMajor):
        #only consider branches that are major correspondences
        if (corrTypes[iSite1]<2):
          continue
      
      #continue by duplicating the entire previous track
      site1 = corrSites[iSite1]
      tracks_checkContinue.append(basinTrack+[site1])
    
  return trackList
  
def form_tracks_iTime(dataCorr, iTimeStart, iTimeEnd, sites0, trackOnlyMajor):
  trackList = []
  for site in sites0:
    #print "Forming tracks from correspondences for initial site {0} at time {1}".format(site, iTimeStart)
    siteTracks = form_track_site(dataCorr, iTimeStart, iTimeEnd, site, trackOnlyMajor)
    for t in siteTracks:
      if (len(t)>1): #drop the noise/artifacts
        trackList.append(t)
    
  return trackList
  
def run_tracks_timeInterval(fNameTracks, fCorr, iTimeStart, iTimeEnd, timeStartGlobal, deltaTGlobal, fMetrics='', trackOnlyMajor=False):
  #find tracks for sites at [iTimeStart,iTimeEnd) out to iTimeEnd (at most).
  #besides stitching the correspondences together, the tricky part is that we don't want to start a track
  #for a basin at each timestep...only start a track at "genesis".
  #define genesis as:
  #-basin not part of an existing track at that time
  
  #store the basins that are part of tracks at each time
  nTimes = iTimeEnd-iTimeStart
  sitesInTrack = [[] for i in xrange(nTimes+1)] #need +1 since storing iTimeStart basins as well
  
  dataTracks = write_tracks_metrics_netcdf_header(fNameTracks, 'test', nTimes, nTimes); iTrackGlobal=0
  dataCorr = netCDF4.Dataset(fCorr, 'r')
  
  for iTime in xrange(nTimes):
    timeInd = iTime+iTimeStart
    sites0, corrSites, typeCorr = correspond.read_corr_iTime(dataCorr, timeInd)
    
    #ignore sites already trajectoried in an existing track
    nSites0 = len(sites0)
    notInPrev = np.ones(nSites0,dtype=int)
    
    for iSite in xrange(nSites0):
      site0 = sites0[iSite]
      if (site0 in sitesInTrack[iTime]):
        notInPrev[iSite] = 0

    #print "{0}/{1} sites started at time {2}".format(np.sum(notInPrev>0), nSites0, timeInd) #this doesn't account for trackOnlyMajor
    sites0 = sites0[notInPrev>0]
  
    trackList = form_tracks_iTime(dataCorr, timeInd, iTimeEnd, sites0,trackOnlyMajor); print "Formed tracks for for iTimeGlobal: ", timeInd
    
    #update sitesInTrack
    for trackSeq in trackList:
      for i in xrange(len(trackSeq)):
        #basin in trackSeq[i] is at index i [iTimeStart,iTimeEnd]
        sitesInTrack[iTime+i].append(trackSeq[i])
    #print sitesInTrack
        
    #write to file
    if (fMetrics==''):
      write_tracks_cells(fNameTracks, trackList)
    else:
      dataMetrics = netCDF4.Dataset(fMetrics,'r')
      
      #write_tracks_metrics_iTime(fNameTracks, timeInd, trackList, dataMetrics, timeStartGlobal, deltaTGlobal)
      iTrackGlobal = write_tracks_metrics_iTime_netcdf(dataTracks, timeInd, iTrackGlobal, trackList, dataMetrics, timeStartGlobal, deltaTGlobal)
      
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
def write_tracks_metrics_netcdf_header(fName, info, nTimesInTrackMax, nTimes):
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
  data.createDimension('nTimesTrack', nTimesInTrackMax+1)
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
  return data

def write_tracks_metrics_iTime_netcdf(data, iTime0, iTrackGlobal0, trackList, dataMetrics, timeStartGlobal, deltaTGlobal):
  
  tStart = timeStartGlobal+deltaTGlobal*iTime0; tStart = tStart.strftime(timeStringFormat)
  data.variables['timeStamp'][iTime0] = tStart
  
  nTracks = len(trackList)
  trackLengths = np.array([len(track) for track in trackList], dtype=int)
  maxLength = np.max(trackLengths)
  print "Maximum track length={0} at time {1}".format(maxLength, tStart)
  
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
  varInd = metricNames.index(varKey); varMin = 280.; varMax = 325.

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
