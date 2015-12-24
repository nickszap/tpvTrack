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
import subprocess
import cPickle as pickle

import tracks
import helpers
import seaborn as sns

r2d = 180./np.pi

def plot_sfc(lat, lon, vals, fNameSave=None, title = ''):
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

def calc_counts(dataSeg, dataTrack):
  #counts of number of basins in each cell
  
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
      nMin = 4; #print "Not counting tracks with nTimes<",nMin
      if (nTimes<nMin):
        continue
    
    sitesInTrack = dataTrack.variables['siteExtr'][iTrack,0:nTimes]
    
    #for iTime in xrange(nTimes):
    for iTime in xrange(1):
    #for iTime in xrange(nTimes-1, nTimes):
      iTimeGlobal = iTime0+iTime
      site = sitesInTrack[iTime]
      
      cell2Site = dataSeg.variables['cell2Site'][iTimeGlobal,:]
      inBasin = cell2Site==site
        
      countIn += inBasin #adds True==1
      
  return countIn

def demo_counts():
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

def demo_plot_minMax_life():
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

def demo_casesTracks():
  #
  
  t0 = '20070720'; tStamp = '2007-07-20-00'
  fERA = '/data01/tracks/summer07/tpvTrack/tracks_low.nc'
  data = netCDF4.Dataset(fERA,'r')
  lenTracks = data.variables['lenTrack'][:]; iTrack = np.argmax(lenTracks);
  timeStamps = data.variables['timeStamp'][:]; iTimeFile = timeStamps.tolist().index(tStamp)
  iTime0 = data.variables['iTimeStart'][iTrack];
  iTimeTrack = iTimeFile-iTime0
  latERA = data.variables['latExtr'][iTrack,iTimeTrack]; lonERA = data.variables['lonExtr'][iTrack,iTimeTrack]; print 'lat,lon for date: ', timeStamps[iTime0], latERA, lonERA, t0
  iTrackERA = iTrack; iTimeTrackERA = iTimeTrack
  data.close()

  fDir = '/data01/tracks/tigge/2007-07/tracks/'
  cmd = 'find {0} -type f | grep tracks_low_horizPlusVert.nc'.format(fDir) #find filenames with 'log.algo' under fDir
  print cmd;
  result = subprocess.check_output(cmd, shell=True); fNames = result.strip().split()
  fNames = [i for i in fNames];
  print fNames

  plt.figure()
  m = Basemap(projection='npstere', boundinglat=60, lon_0=0, resolution='l')#, round=True)
  m.drawcoastlines()
  
  maxTimesTrack = 0
  for f in fNames:
    data = netCDF4.Dataset(f,'r')

    #follow the tpv with longest life near the starting location (from ERA-I)...is there a better way?
    candTracks = np.where(data.variables['iTimeStart'][:]==0)[0]
    latCands = data.variables['latExtr'][candTracks,0]; lonCands = data.variables['lonExtr'][candTracks,0]
    d = helpers.calc_distSphere_multiple(6370., latERA*np.pi/180, lonERA*np.pi/180, latCands*np.pi/180, lonCands*np.pi/180)
    candTracks = candTracks[d<200]

    lenTracks = data.variables['lenTrack'][candTracks]
    iTrack = np.argmax(lenTracks); iTrack = candTracks[iTrack]; print 'iTrack: ', iTrack, f

    nTimesTrack = data.variables['lenTrack'][iTrack]
    latTrack = data.variables['latExtr'][iTrack,0:nTimesTrack]
    lonTrack = data.variables['lonExtr'][iTrack,0:nTimesTrack]
    data.close()
    maxTimesTrack = max(nTimesTrack,maxTimesTrack)

    x,y = m(lonTrack,latTrack)
    m.scatter(x[0],y[0], marker='+', color='g', s=100)
    m.scatter(x[-1],y[-1], marker='o', color='r', s=50)
    #lineStyle, wdth = get_ensembleStyle(f)
    m.plot(x,y) #, lineStyle)
  
  data = netCDF4.Dataset(fERA,'r')
  latERA = data.variables['latExtr'][iTrackERA,iTimeTrackERA:iTimeTrackERA+maxTimesTrack]; lonERA = data.variables['lonExtr'][iTrackERA,iTimeTrackERA:iTimeTrackERA+maxTimesTrack];
  data.close()
  x,y = m(lonERA,latERA)
  m.scatter(x[0],y[0], marker='+', color='g', s=100)
  m.scatter(x[-1],y[-1], marker='o', color='r', s=50)
  #lineStyle, wdth = get_ensembleStyle(f)
  m.plot(x,y, 'k-', linewidth=5.0)
    
  plt.show()

def demo_error_ensTPVs():

  key = 'thetaExtr'; units = 'K';
  tStampAnl = '2006-11-10-06' #time the TPV started
  tStampTIGGE = '2006-11-15-00' # #we start tracking w/in TIGGE after the start date of anl TPV (allows TPV to be 'well formed'). no constant dtIncrement since ~5 days later at 00Z
  fERA = '/data01/tracks/1979-2015/tracks/tracks_low.nc'
  
  data = netCDF4.Dataset(fERA,'r')
  lenTracks = data.variables['lenTrack'][:]; #iTrack = np.argmax(lenTracks)
  timeStamps = data.variables['timeStamp'][:]; iTimeFile = timeStamps.tolist().index(tStampAnl)
  iTimeStart = data.variables['iTimeStart'][:]
  
  iTrack = np.argmax(lenTracks* (iTimeStart==iTimeFile) ) #longest lived TPV that started at that startTime
  
  iTimeFile = timeStamps.tolist().index(tStampTIGGE) #align analysis tpv time to tigge tpv time
  iTime0 = iTimeStart[iTrack];
  iTimeTrack = iTimeFile-iTime0
  latERA = data.variables['latExtr'][iTrack,iTimeTrack]; lonERA = data.variables['lonExtr'][iTrack,iTimeTrack]; print 'lat,lon for date: ', tStampAnl, latERA, lonERA
  valsRef = data.variables[key][iTrack,iTimeTrack:]; print 'valsRef: ', valsRef
  data.close()

  fDir = '/raid1/nick/tpvTracks/tigge/{0}/'.format(tStampAnl)
  cmd = 'find {0} -type f | grep tracks_low_horizPlusVert.nc'.format(fDir) #find filenames with 'log.algo' under fDir
  print cmd;
  result = subprocess.check_output(cmd, shell=True); fNames = result.strip().split()
  fNames = [i for i in fNames] # if t0 in i];
  print fNames

  plt.figure()

  for f in fNames:
    data = netCDF4.Dataset(f,'r')

    #follow the tpv with longest life near the starting location (from ERA-I)...is there a better way?
    candTracks = np.where(data.variables['iTimeStart'][:]==0)[0]
    latCands = data.variables['latExtr'][candTracks,0]; lonCands = data.variables['lonExtr'][candTracks,0]
    d = helpers.calc_distSphere_multiple(6370., latERA*np.pi/180, lonERA*np.pi/180, latCands*np.pi/180, lonCands*np.pi/180)
    candTracks = candTracks[d<200]

    lenTracks = data.variables['lenTrack'][candTracks]
    iTrack = np.argmax(lenTracks); iTrack = candTracks[iTrack]; print 'iTrack: ', iTrack

    nTimesTrack = data.variables['lenTrack'][iTrack]
    valsTrack = data.variables[key][iTrack,0:nTimesTrack];
    print f, key, valsTrack.tolist(), data.variables['latExtr'][iTrack,0:nTimesTrack].tolist(), data.variables['lonExtr'][iTrack,0:nTimesTrack].tolist()
    #dVals = valsTrack
    #dVals = valsTrack-valsRef[0:nTimesTrack]
    dVals = valsTrack-valsTrack[0]-(valsRef[0:nTimesTrack]-valsRef[0])

    #lineStyle, wdth = get_ensembleStyle(f)
    plt.plot(dVals) #, lineStyle)
    
    data.close()
  plt.xlabel('Time')
  plt.ylabel(key+', '+units)
  plt.show()

def demo_groupEnsemble_variable():
  #put metric for analysis, perturbed ensemble, control ensemble, and high-res TPV into the same file for easier analysis.
  #Ensemble is not ordered 0,1,...,49
  
  #specific to each case/variable ------------------------------
  key = 'latExtr';
  
  #anlTimes = ['2007-03-14-12', '2007-04-15-18', '2007-05-01-12', '2007-05-21-00', '2007-06-13-18', '2007-07-13-18', '2007-08-30-18', '2007-12-31-06', '2008-04-27-12', '2008-05-20-00', '2008-06-08-18', '2008-07-22-06', '2008-07-23-00', '2008-08-17-18', '2008-12-02-18']
  #anlTimes = ['2007-03-14-12']
  anlTimes  = ['2006-11-10-06', '2006-11-13-12', '2006-11-20-00', '2006-12-22-06', '2006-12-27-18', '2007-03-14-12', '2007-04-15-18', '2007-05-01-12', '2007-05-21-00', '2007-06-13-18', '2007-07-13-18', '2007-08-30-18', '2007-12-31-06', '2008-04-27-12', '2008-05-20-00', '2008-06-08-18', '2008-07-22-06', '2008-07-23-00', '2008-08-17-18', '2008-12-02-18', '2009-01-29-06', '2009-03-24-12', '2009-04-28-06', '2009-05-09-18', '2009-06-17-18', '2009-08-05-18', '2009-08-29-12', '2009-09-25-00', '2009-10-19-06', '2009-11-18-18', '2010-01-31-18', '2010-04-01-06', '2010-04-26-06', '2010-05-07-06', '2010-05-26-18', '2010-06-25-18', '2010-08-03-18', '2010-08-19-18', '2011-01-13-06', '2011-03-21-18', '2011-03-27-18', '2011-05-06-18', '2011-06-04-18', '2011-09-19-18', '2011-09-20-06', '2011-12-20-18', '2012-01-20-18', '2012-03-06-00', '2012-05-05-18', '2012-06-11-06', '2012-06-18-06', '2012-08-06-18', '2013-02-01-00', '2013-05-18-00', '2013-06-30-00', '2013-09-14-06', '2013-12-24-18', '2014-04-20-18', '2014-07-12-00', '2014-08-05-12', '2014-09-22-18', '2014-11-09-06', '2015-03-18-06', '2015-06-17-06', '2015-06-19-12', '2015-06-27-18']
  
  #tiggeTimes = ['2007-03-19-00', '2007-04-20-00', '2007-05-06-00', '2007-05-26-00', '2007-06-18-00', '2007-07-18-00', '2007-09-04-00', '2008-01-05-00', '2008-05-02-00', '2008-05-25-00', '2008-06-13-00', '2008-07-27-00', '2008-07-28-00', '2008-08-22-00', '2008-12-07-00']
  #tiggeTimes = ['2007-03-19-00']
  tiggeTimes = ['2006-11-15-00', '2006-11-18-00', '2006-11-25-00', '2006-12-27-00', '2007-01-01-00', '2007-03-19-00', '2007-04-20-00', '2007-05-06-00', '2007-05-26-00', '2007-06-18-00', '2007-07-18-00', '2007-09-04-00', '2008-01-05-00', '2008-05-02-00', '2008-05-25-00', '2008-06-13-00', '2008-07-27-00', '2008-07-28-00', '2008-08-22-00', '2008-12-07-00', '2009-02-03-00', '2009-03-29-00', '2009-05-03-00', '2009-05-14-00', '2009-06-22-00', '2009-08-10-00', '2009-09-03-00', '2009-09-30-00', '2009-10-24-00', '2009-11-23-00', '2010-02-05-00', '2010-04-06-00', '2010-05-01-00', '2010-05-12-00', '2010-05-31-00', '2010-06-30-00', '2010-08-08-00', '2010-08-24-00', '2011-01-18-00', '2011-03-26-00', '2011-04-01-00', '2011-05-11-00', '2011-06-09-00', '2011-09-24-00', '2011-09-25-00', '2011-12-25-00', '2012-01-25-00', '2012-03-11-00', '2012-05-10-00', '2012-06-16-00', '2012-06-23-00', '2012-08-11-00', '2013-02-06-00', '2013-05-23-00', '2013-07-05-00', '2013-09-19-00', '2013-12-29-00', '2014-04-25-00', '2014-07-17-00', '2014-08-10-00', '2014-09-27-00', '2014-11-14-00', '2015-03-23-00', '2015-06-22-00', '2015-06-24-00', '2015-07-02-00']
  
  nCases = len(anlTimes)
  if (len(tiggeTimes) != nCases):
    print 'Uhoh. # times dont match between analysis and tigge'
    exit()
  
  for iCase in xrange(nCases):
    #tStampAnl = '2006-11-10-06' #time the analysis TPV started
    #tStampTIGGE = '2006-11-15-00' #we start tracking w/in TIGGE after the start date of anl TPV (allows TPV to be 'well formed'). no constant dtIncrement since ~5 days later at 00Z
    tStampAnl = anlTimes[iCase]
    tStampTIGGE = tiggeTimes[iCase]
    
    fERA = '/data01/tracks/1979-2015/tracks/tracks_low.nc'
    
    fDir = '/data02/cases/tigge/{0}/'.format(tStampAnl)
    cmd = 'find {0} -type f | grep tracks_low_horizPlusVert.nc'.format(fDir) #find filenames with 'log.algo' under fDir
    print cmd;
    result = subprocess.check_output(cmd, shell=True); fNames = result.strip().split()
    fNamesEns = [i for i in fNames] # if t0 in i];
    #fNames = sorted(fNames)
    print fNames
    
    fNameOut = fDir+key+'.pkl' #'.npz'
    #------------------------------
    
    #analysis TPV
    data = netCDF4.Dataset(fERA,'r')
    lenTracks = data.variables['lenTrack'][:]; #iTrack = np.argmax(lenTracks)
    timeStamps = data.variables['timeStamp'][:]; iTimeFile = timeStamps.tolist().index(tStampAnl)
    iTimeStart = data.variables['iTimeStart'][:]
    
    iTrack = np.argmax(lenTracks* (iTimeStart==iTimeFile) ) #longest lived TPV that started at that startTime
    
    iTimeFile = timeStamps.tolist().index(tStampTIGGE) #align analysis tpv time to tigge tpv time
    iTime0 = iTimeStart[iTrack];
    iTimeTrack = iTimeFile-iTime0
    latERA = data.variables['latExtr'][iTrack,iTimeTrack]; lonERA = data.variables['lonExtr'][iTrack,iTimeTrack]; print 'lat,lon for date: ', tStampAnl, latERA, lonERA
    nTimesTrack = lenTracks[iTrack]
    valsAnl = data.variables[key][iTrack,iTimeTrack:nTimesTrack]; #print 'valsRef: ', valsRef
    data.close()
    #------------------------------
    
    #ensemble TPV ------------------------------
    nMembers = len(fNamesEns)
    valsEns = [None]*nMembers #since variable number of times in each member. can np.asarray(valsEns) later if you'd like
    #for f in fNamesEns:
    for iMember in xrange(nMembers):
      f = fNamesEns[iMember]; print f
      data = netCDF4.Dataset(f,'r')

      #follow the tpv with longest life near the starting location (from ERA-I)...is there a better way?
      candTracks = np.where(data.variables['iTimeStart'][:]==0)[0]
      latCands = data.variables['latExtr'][candTracks,0]; lonCands = data.variables['lonExtr'][candTracks,0]
      d = helpers.calc_distSphere_multiple(6370., latERA*np.pi/180, lonERA*np.pi/180, latCands*np.pi/180, lonCands*np.pi/180)
      candTracks = candTracks[d<300]
      
      if (len(candTracks)==0):
        print 'No candidate matching tracks for this case'
        continue
      lenTracks = data.variables['lenTrack'][candTracks]
      iTrack = np.argmax(lenTracks); iTrack = candTracks[iTrack]; print 'iTrack: ', iTrack

      nTimesTrack = data.variables['lenTrack'][iTrack]
      valsTrack = data.variables[key][iTrack,0:nTimesTrack];
      
      valsEns[iMember] = valsTrack
      data.close()
    #------------------------------
    
    print 'Saving to: ', fNameOut
    pickle.dump([valsAnl, valsEns], open(fNameOut,'wb') ) #read back in with pickle.load( open(fNameOut,'rb') )
    #np.savez(fNameOut,anl=valsAnl, ens=valsEns) #savez not implemented

def demo_valsFromPickle():
  
  fDir = '/data02/cases/tigge/'
  #varInfo = 'Relative $\Delta\Theta_{core}$ (K)'
  #cmd = 'find {0} -type f | grep thetaExtr.pkl'.format(fDir) #find filenames with 'log.algo' under fDir
  varInfo = 'Radius (km)'
  cmd = 'find {0} -type f | grep rEquiv.pkl'.format(fDir) #find filenames with 'log.algo' under fDir
  print cmd;
  result = subprocess.check_output(cmd, shell=True); fNames = result.strip().split()
  fNamesCases = [i for i in fNames] # if t0 in i];
  print fNamesCases
  
  plt.figure()
  
  for f in fNamesCases:
    valsAnl, valsEns = pickle.load(open(f,'rb'))
    
    for vals in valsEns:
      if (vals != None): #couldn't find matching track
        nTimesTrack = len(vals)
        dVals = vals-vals[0]-(valsAnl[0:nTimesTrack]-valsAnl[0])

        #lineStyle, wdth = get_ensembleStyle(f)
        plt.plot(dVals) #, lineStyle)
        
  plt.xlabel('Time')
  plt.ylabel(varInfo)
  plt.show()

def demo_valsFromPickle_reliability():
  
  fDir = '/data02/cases/tigge/'
  #cmd = 'find {0} -type f | grep thetaExtr.pkl'.format(fDir) #find filenames with 'log.algo' under fDir
  #cmd = 'find {0} -type f | grep rEquiv.pkl'.format(fDir) #find filenames with 'log.algo' under fDir
  cmd = 'find {0} -type f | grep latExtr.pkl'.format(fDir) #find filenames with 'log.algo' under fDir
  print cmd;
  result = subprocess.check_output(cmd, shell=True); fNames = result.strip().split()
  fNamesCases = [i for i in fNames] # if t0 in i];
  print fNamesCases
  
  #iTime0 = 0; iTime1 = 4
  iTime0 = 16; iTime1 = max(8,iTime0)
  xVals = []
  yVals = []
  
  for f in fNamesCases:
    valsAnl, valsEns = pickle.load(open(f,'rb'))
    
    #valAnlCase = valsAnl[iTime1]-valsAnl[iTime0]
    valAnlCase = valsAnl[iTime0]
    for vals in valsEns:
      if (vals != None): #couldn't find matching track
        nTimesTrack = len(vals)
        if (nTimesTrack<=iTime1): #matching track not long enough
          continue
        #valMember = vals[iTime1]-vals[iTime0]
        valMember = vals[iTime0]

        xVals.append(valAnlCase)
        yVals.append(valMember)
  
  plt.figure()
  plt.plot(xVals,yVals,'bo')
  plt.plot(xVals,xVals,'k-')
  plt.xlabel('Analysis')
  plt.ylabel('Ensemble')
  #plt.ylabel(varInfo)
  plt.show()

def demo_valsFromPickle_rankHistogram():
  
  fDir = '/data02/cases/tigge/'
  cmd = 'find {0} -type f | grep thetaExtr.pkl'.format(fDir) #find filenames with 'log.algo' under fDir
  print cmd;
  result = subprocess.check_output(cmd, shell=True); fNames = result.strip().split()
  fNamesCases = [i for i in fNames] # if t0 in i];
  print fNamesCases
  
  iTime0 = 0; iTime1 = 4
  valsRank = []
  
  for f in fNamesCases:
    xVals = []
    valsEnsCase = []
    
    valsAnl, valsEns = pickle.load(open(f,'rb'))
    
    valAnlCase = valsAnl[iTime1]-valsAnl[iTime0]
    if (True):
      if (valAnlCase>0): #only look at intensification
        continue
    
    for vals in valsEns:
      if (vals != None): #couldn't find matching track
        nTimesTrack = len(vals)
        if (nTimesTrack<=iTime1): #matching track not long enough
          continue
        valMember = vals[iTime1]-vals[iTime0]
        valsEnsCase.append(valMember)
    
    valsEnsCase = sorted(valsEnsCase); nVals = len(valsEnsCase) #sorted goes [low,high] so low indices are more negative (e.g., more intense for theta-theta0)
    iAnl = np.searchsorted(valsEnsCase, valAnlCase)
    #valsRank.append(iAnl/float(nVals-1) )
    valsRank.append(iAnl/float(nVals) ) #since really nVals+1 slots with ensembleMembers+analysisMember
    
    print valsEnsCase, valAnlCase, iAnl
    
  plt.hist(valsRank)
  plt.show()
  
if __name__=='__main__':
  #demo_counts()
  #demo_plot_minMax_life()
  #demo_compareMetrics_pretty()
  #demo_calendarLife()
  #demo_calendarMetric()
  #demo_casesTracks()
  #demo_error_ensTPVs()
  #demo_groupEnsemble_variable()
  #demo_valsFromPickle()
  demo_valsFromPickle_reliability()
  #demo_valsFromPickle_rankHistogram()



