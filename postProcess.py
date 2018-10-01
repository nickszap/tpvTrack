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

import tracks
import seaborn as sns

r2d = 180./np.pi

def plot_sfc(lat, lon, vals, fNameSave=None, title = ''):
  """
  Color plot on a map
  
  Arguments:
  lat - latitudes (in radians)
  lon - longitudes (in radians)
  vals - values to color
  fNameSave - if not None, name of figure to save
  title - title to add to plot
  """
  m = Basemap(projection='ortho',lon_0=0,lat_0=89.95, resolution='l')
  x,y = m(lon*r2d, lat*r2d)
  #print x.shape, y.shape

  m.drawcoastlines(linewidth=.5)
  #m.drawmapboundary()
  
  pPlot = m.pcolor(x,y,vals,tri=True, shading='flat',edgecolors='none',cmap=plt.cm.RdBu_r) #, vmin=280, vmax=360)

  plt.colorbar(pPlot)
  plt.title(title)
  if (fNameSave == None):
    plt.show()
  else:
    plt.savefig(fNameSave, bbox_inches='tight'); plt.close()

def calc_counts(dataSeg, dataTrack,nTimesMin=4):
  """Returns the counts of number of times a tracked basin (with lifetime>=nTimesMin) occurs in each cell"""
  
  #for composite-ing field using the basin as a filter
  nCells = len(dataSeg.dimensions['nCells'])
  countIn = np.zeros(nCells, dtype=int)
  
  nTracks = len(dataTrack.dimensions['nTracks'])
  for iTrack in xrange(nTracks):
    #tracking the cells in the basin to define the horizontal filter
    
    #progress
    if (iTrack%500 == 0):
      print "On track {0}/{1}".format(iTrack,nTracks)
    
    iTime0 = dataTrack.variables['iTimeStart'][iTrack]
    nTimes = dataTrack.variables['lenTrack'][iTrack]
    if (True):
      #nMin = 4; #print "Not counting tracks with nTimes<",nMin
      if (nTimes<nTimesMin):
        continue
    
    sitesInTrack = dataTrack.variables['siteExtr'][iTrack,0:nTimes]
    
    #for iTime in xrange(nTimes): #all times in track
    for iTime in xrange(1): #genesis times in track
    #for iTime in xrange(nTimes-1, nTimes): #lysis times in track
      iTimeGlobal = iTime0+iTime
      site = sitesInTrack[iTime]
      
      cell2Site = dataSeg.variables['cell2Site'][iTimeGlobal,:]
      inBasin = cell2Site==site
        
      countIn += inBasin #adds True==1
      
  return countIn

def demo_counts():
  """Example of plotting basin counts"""
  fDir = '/data02/tracks/summer07/tpvTrack/' #'/data02/tracks/summer06/jun1-sep30/'
  fMetr = fDir+'fields.nc'
  fSeg = fDir+'seg.nc'
  fTrack = fDir+'tracks_high.nc'
  
  dataSeg = netCDF4.Dataset(fSeg, 'r')
  dataTrack = netCDF4.Dataset(fTrack,'r')
  
  counts = calc_counts(dataSeg, dataTrack)
  
  dataMetr = netCDF4.Dataset(fMetr,'r')
  lat = dataMetr.variables['latCell'][:]
  lon = dataMetr.variables['lonCell'][:]
  
  plot_sfc(lat, lon, counts, fNameSave=None, title = 'Genesis count of 2007 highs')
    
def demo_plotMetrics(fTracks):
  """Example of plotting tracks' metrics"""
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
  """Example of joint plotting 2 tracks' metrics"""
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
  """Example histogram of tracked lifetimes"""
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

def demo_plot_minMax_life():
  """Example of plotting when metric extrema occur within each track"""
  filesTracks = ['/data02/tracks/summer06/jun1-sep30/tracks_test_tanh_.66.nc', '/data02/tracks/summer06/jun1-sep30/tracks_high.nc',
                 '/data02/tracks/summer07/tpvTrack/tracks_low.nc', '/data02/tracks/summer07/tpvTrack/tracks_high.nc']
  info = ['lows2006', 'highs2006', 'lows2007', 'highs2007']
  key = 'rEquiv'
  for iFile,fTrack in enumerate(filesTracks):
    mins = []; iTimeMins = [];
    maxs = []; iTimeMaxs = [];
  
    dataTrack = netCDF4.Dataset(fTrack,'r')
    
    nTracks = len(dataTrack.dimensions['nTracks'])
    for iTrack in xrange(nTracks):
      #progress
      if (iTrack%500 == 0):
        print "On track {0}/{1}".format(iTrack,nTracks)
      
      iTime0 = dataTrack.variables['iTimeStart'][iTrack]
      nTimes = dataTrack.variables['lenTrack'][iTrack]
      if (True):
        nMin = 16; #print "Not counting tracks with nTimes<",nMin
        if (nTimes<nMin):
          continue
      
      vals = dataTrack.variables[key][iTrack,0:nTimes]
      iMin = np.argmin(vals); iTimeMins.append(iMin/float(nTimes)); mins.append(vals[iMin])
      iMax = np.argmax(vals); iTimeMaxs.append(iMax/float(nTimes)); maxs.append(vals[iMax])
    
    iTimeMins = np.array(iTimeMins); mins = np.array(mins)
    iTimeMaxs = np.array(iTimeMaxs); maxs = np.array(maxs)
    print fTrack
    if (False):  
      plt.figure()
      plt.plot(iTimeMins, mins, 'bv')
      plt.plot(iTimeMaxs, maxs, 'r^')
      plt.show()
    else:
      #plt.figure()
      with sns.axes_style("white"):
        f = sns.jointplot(iTimeMins, mins, kind="hex", color='blue');
        #f.set_axis_labels(['Normalized lifetime', key])
        sns.axlabel('Normalized lifetime', key)
        saveName = 'min_{0}_{1}_{2}.png'.format(key, nMin, info[iFile])
        plt.savefig(saveName); plt.close()
        
        f = sns.jointplot(iTimeMaxs, maxs, kind="hex", color='red');
        #f.set_axis_labels(['Normalized lifetime', key])
        sns.axlabel('Normalized lifetime', key)
        saveName = 'max_{0}_{1}_{2}.png'.format(key, nMin, info[iFile])
        plt.savefig(saveName); plt.close()

def demo_compareMetrics_pretty():
  """Example of joint-plotting 2 metrics """
  filesTracks = ['/data02/tracks/summer06/jun1-sep30/tracks_test_tanh_.66.nc', '/data02/tracks/summer06/jun1-sep30/tracks_high.nc',
                 '/data02/tracks/summer07/tpvTrack/tracks_low.nc', '/data02/tracks/summer07/tpvTrack/tracks_high.nc']
  info = ['lows2006', 'highs2006', 'lows2007', 'highs2007']
  keys = ['rEquiv', 'ampMaxMin']
  for iFile,fTrack in enumerate(filesTracks):
    mins1 = []; mins2 = []
    maxs1 = []; maxs2 = []
  
    dataTrack = netCDF4.Dataset(fTrack,'r')
    
    nTracks = len(dataTrack.dimensions['nTracks'])
    for iTrack in xrange(nTracks):
      #progress
      if (iTrack%500 == 0):
        print "On track {0}/{1}".format(iTrack,nTracks)
      
      iTime0 = dataTrack.variables['iTimeStart'][iTrack]
      nTimes = dataTrack.variables['lenTrack'][iTrack]
      if (True):
        nMin = 16; #print "Not counting tracks with nTimes<",nMin
        if (nTimes<nMin):
          continue
      
      vals1 = dataTrack.variables[keys[0]][iTrack,0:nTimes]
      vals2 = dataTrack.variables[keys[1]][iTrack,0:nTimes]
      iMin = np.argmin(vals1)
      mins1.append(vals1[iMin]); mins2.append(vals2[iMin])
      iMax = np.argmin(vals1)
      maxs1.append(vals1[iMax]); maxs2.append(vals2[iMax])
    
    mins1 = np.array(mins1); mins2 = np.array(mins2);
    maxs1 = np.array(maxs1); maxs2 = np.array(maxs2);
    print fTrack
    if (False):  
      plt.figure()
      plt.plot(iTimeMins, mins, 'bv')
      plt.plot(iTimeMaxs, maxs, 'r^')
      plt.show()
    else:
      #plt.figure()
      with sns.axes_style("white"):
        f = sns.jointplot(mins1, mins2, kind="hex", color='blue');
        f.set_axis_labels(keys)
        #sns.axlabel(keys)
        saveName = 'min_{0}_{1}_{2}.png'.format(':'.join(keys), nMin, info[iFile])
        plt.savefig(saveName); plt.close()
        
        f = sns.jointplot(maxs1, maxs2, kind="hex", color='red');
        f.set_axis_labels(keys)
        #sns.axlabel(keys)
        saveName = 'max_{0}_{1}_{2}.png'.format(':'.join(keys), nMin, info[iFile])
        plt.savefig(saveName); plt.close()
        
def demo_calendarLife():
  """Example of plotting longest track that starts at each time"""
  filesTracks = ['/data02/tracks/summer06/jun1-sep30/tracks_test_tanh_.66.nc', '/data02/tracks/summer06/jun1-sep30/tracks_high.nc',
                 '/data02/tracks/summer07/tpvTrack/tracks_low.nc', '/data02/tracks/summer07/tpvTrack/tracks_high.nc']
  info = ['lows2006', 'highs2006','lows2007', 'highs2007']
  
  f = filesTracks[0]
  data = netCDF4.Dataset(f,'r')
  tStrings = data.variables['timeStamp'][:]
  data.close()
  times = [dt.datetime.strptime(i, tracks.timeStringFormat) for i in tStrings]
  times = np.array(times)
  deltaT = times[1]-times[0]; deltaT = deltaT.total_seconds()/(60*60*24); units = 'days';
  
  nTimes = len(times);
  nFiles = len(filesTracks);
  vals = np.zeros((nFiles,nTimes))
  
  for iFile in xrange(nFiles):
    data = netCDF4.Dataset(filesTracks[iFile],'r')
    
    nTracks = len(data.dimensions['nTracks'])
    iTimeStart = data.variables['iTimeStart'][:]
    lenTrack = data.variables['lenTrack'][:]
    for iTrack in xrange(nTracks):
      vals[iFile,iTimeStart[iTrack]] = max(vals[iFile,iTimeStart[iTrack]], lenTrack[iTrack])
      
    data.close()
  vals *= deltaT
  
  plt.figure()
  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d')) #('%Y/%m/%d/%H'))
  plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=7))
  plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
  #a.autofmt_xdate()
  
  for iFile in xrange(nFiles):
    plt.plot(times, vals[iFile,:], label=info[iFile])
  plt.legend()
  plt.xlabel("Start date")
  plt.ylabel("Lifetime ({0})".format(units))
  plt.show()

def demo_calendarMetric():
  """Example of plotting maximum metric of track at each time"""
  filesTracks = ['/data02/tracks/summer06/jun1-sep30/tracks_test_tanh_.66.nc', '/data02/tracks/summer06/jun1-sep30/tracks_high.nc',
                 '/data02/tracks/summer07/tpvTrack/tracks_low.nc', '/data02/tracks/summer07/tpvTrack/tracks_high.nc']
  info = ['lows2006', 'highs2006','lows2007', 'highs2007']
  
  f = filesTracks[0]
  data = netCDF4.Dataset(f,'r')
  tStrings = data.variables['timeStamp'][:]
  data.close()
  times = [dt.datetime.strptime(i, tracks.timeStringFormat) for i in tStrings]
  times = np.array(times)
  
  key = 'rEquiv'; units = 'km'
  
  nTimes = len(times);
  nFiles = len(filesTracks);
  vals = np.zeros((nFiles,nTimes), dtype=float)
  
  for iFile in xrange(nFiles):
    data = netCDF4.Dataset(filesTracks[iFile],'r')
    
    nTracks = len(data.dimensions['nTracks'])
    iTimeStart = data.variables['iTimeStart'][:]
    lenTrack = data.variables['lenTrack'][:]
    for iTrack in xrange(nTracks):
      nTrackTimes = lenTrack[iTrack]
      valsTrack = data.variables[key][iTrack,0:nTrackTimes]
      
      for iTime in xrange(nTrackTimes):
        iTimeGlobal = iTimeStart[iTrack]+iTime
        if (iTimeGlobal>=nTimes):
          print "Uhoh. longer time than possible: ", iTimeGlobal, nTimes
          continue
      
        vals[iFile,iTimeGlobal] = max(vals[iFile,iTimeGlobal], valsTrack[iTime])
      
    data.close()
  
  plt.figure()
  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d')) #('%Y/%m/%d/%H'))
  plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=7))
  plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
  #a.autofmt_xdate()
  
  for iFile in xrange(nFiles):
    plt.plot(times, vals[iFile,:], label=info[iFile])
  plt.legend()
  plt.xlabel("Date")
  plt.ylabel("{0} ({1})".format(key, units))
  plt.show()

if __name__=='__main__':
  demo_counts()
  #demo_plot_minMax_life()
  #demo_compareMetrics_pretty()
  #demo_calendarLife()
  #demo_calendarMetric()


