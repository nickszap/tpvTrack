#The track program creates correspondences between basins at neighboring times.
#We want tracks that extend over the lifetime of each/all basins

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
import glob
import datetime as dt
import os

import segment_ll_flat
import basinMetrics

r = 6370.e3
r2d = 180./np.pi

def get_correspondingSites(site0, sites0, sites1, isMatch):
  #isMatch[site0, site1]>0 if basins correspond between frames.
  #return list of sites1 corresponding to site0
  
  ind0 = sites0==site0
  if (not ind0.any()):
    return []
  isMatch01 = isMatch[ind0,:].squeeze()>0
  #print isMatch01.shape
  basins1 = sites1[isMatch01]
  return list(basins1)

def gather_timeCorrespondencesBasin(filesTrack, site0):
  #using the basin correspondences between time frames, accumulate over many timesteps.
  #return list of lists with corresponding basins at each list[iTime]
  
  nTimes = len(filesTrack);
  basinTrack = [[] for i in xrange(nTimes)] #list of nTimes empty lists
  basinTrack[0] = [site0]
  
  for iTime in xrange(nTimes-1):
    #file 0 links times 0-1
    #basinTrack[1] is basins at time 1 linked to basin[0]
    
    fTrack = filesTrack[iTime]
    data = np.load(fTrack)
    sites0 = data['sites0']; sites1 = data['sites1']; isMatch = data['isMatch']
    data.close()
    
    #have to follow all basins at previous time
    basins1 = []
    for iBasin in basinTrack[iTime]:
      basins1.extend(get_correspondingSites(iBasin, sites0, sites1, isMatch))
    basins1 = list(set(basins1)); print basins1, '\n'
    
    if (len(basins1)<1):
      print "Track for site {0} ends at timeIndex {1}".format(site0, iTime)
      break
    basinTrack[iTime+1] = basins1
  #print basinTrack
  
  return basinTrack

def gather_1To1Track(filesTrack, site0):
  #duplicate/split tracks of the timeCorrespondenceTree so (1) each tpv goes to <=1 other and (2) each tpv comes from <=1 other
  #in a sequence ordered over time.
  
  nTimes = len(filesTrack);
  basinTrack = []
  basinTrack.append([site0])
  
  for iTime in xrange(nTimes-1):
    #file 0 links times 0-1
    
    fTrack = filesTrack[iTime]
    data = np.load(fTrack)
    sites0 = data['sites0']; sites1 = data['sites1']; isMatch = data['isMatch']
    data.close()
    
    nTracks = len(basinTrack)
    for iTrack in xrange(nTracks):
      #track for timeIndex=0 (correspondences 0->1) should have 1 basin in sequence, the seed basin.
      #if less, the tpv lysis'ed
      trackSeq = basinTrack[iTrack]
      if (len(trackSeq)<iTime-1):
        continue
      #see if can continue on from basin at last time in track
      iBasin = trackSeq[-1]
      basins1 = get_correspondingSites(iBasin, sites0, sites1, isMatch)
      nBasins = len(basins1)
      if (nBasins<1):
        #track ended
        continue
      else:
        #can't pop this track off since mucks with looping over basinTracks, so replace instead
        iBasin = basins1[0]
        newSeq = []; newSeq.extend(trackSeq+[iBasin]); #copy trackSeq and add next basin onto the end
        basinTrack[iTrack] = newSeq
        #append the rest to after nTracks
        for iBasin in basins1[1:]:
          newSeq = []; newSeq.extend(trackSeq+[iBasin]); #copy trackSeq and add next basin onto the end
          basinTrack.append(newSeq)
    #end iTrack
  #end iTime
  if (False):
    for trackSeq in basinTrack:
      print trackSeq
  
  return basinTrack

def trackBasin_mostSimilar(fMesh, filesTrack, site0):
  #return array of the most similar basin track over time
  
  #gather possible basin correspondences over time ---------
  basinCorr = gather_timeCorrespondencesBasin(filesTrack, site0)
  
  #to track an individual tpv, want more than just correspondences.
  #so, we filter choices a la similarity cost function.
  dataMesh = np.load(fMesh)
  lat = dataMesh['lat']; lon = dataMesh['lon']; #nLat = len(lat); nLon = len(lon)
  dataMesh.close()
  
  basinTrack = filterTimeCorrespondence(basinCorr, lat, lon, r)
  return np.array(basinTrack)

def filterTimeCorrespondence(basinCorr, lat, lon, r):
  #suppose we want to have a track with 1 basin at each time frame.
  #given candidate basin correspondences over time, we can select the most similar
  #basin as the appropriate one.
  
  nTimes = len(basinCorr)
  basinTrack = []
  basinTrack.append(basinCorr[0][0])
  for iTime in xrange(1,nTimes):
    candBasins = basinCorr[iTime]; nCandidates = len(candBasins)
    if (nCandidates<1):
      print "Track ends before timeIndex {0}".format(iTime)
      break
    elif (nCandidates==1):
      basinTrack.append(candBasins[0])
    else:
      #choose the "most similar" basin to continue the track
      simVals = np.empty(nCandidates, dtype=float)
      site0 = basinTrack[iTime-1]
      for iCand in xrange(nCandidates):
        site1 = candBasins[iCand]
        simVals[iCand] = basinSimilarity(site0, site1, lat, lon, r)
      iCand = np.argmin(simVals)
      basinTrack.append(candBasins[iCand])
  print basinTrack
  
  return basinTrack
  
def basinSimilarity(site0, site1, lat, lon, r):
  #cost function approach where min value over basins is most similar.
  #many options for what to include in cost when comparing to basin at previous time:
  #-distance, distance from advected location, min(theta), ...
  
  nLon = len(lon)
  iLat0, iLon0 = segment_ll_flat.index_1dTo2d(site0, nLon)
  iLat1, iLon1 = segment_ll_flat.index_1dTo2d(site1, nLon)
  
  d = segment_ll_flat.calc_distSphere_multiple(r, lat[iLat0], lon[iLon0], lat[iLat1], lon[iLon1])
  
  return d

def writeFile_tracks(fSave, basinTrack, filesTrack, filesMetrics):
  '''
  format is:
  (header) nTracks metric1 metric2 ...
  nTimesTrack1
  (time1) metric1 metric2 ... for track1
  (time2) metric1 metric2 ... for track1
  -1
  nTimesTrack2
  (time1) metric1 metric2 ... for track2
  (time2) metric1 metric2 ... for track2
  -1
  .
  .
  .
  '''
  
  f = open(fSave,'w')
  nTracks = len(basinTrack)
  
  #write header
  fMetric = filesMetrics[0]
  data = np.load(fMetric)
  metricNames = [i for i in data]
  data.close()
  headerString = str(nTracks)+' '+' '.join(metricNames)
  f.write(headerString+'\n')
  
  s = ''
  for iTrack in xrange(nTracks):
    trackSeq = basinTrack[iTrack]
    nTimes = len(trackSeq)
    s += str(nTimes)+'\n'
    for iTime in xrange(nTimes):
      fTrack = filesTrack[iTime]
      fMetric = filesMetrics[iTime]
      
      #figure out index of this basin at this time in the metrics data
      site = trackSeq[iTime]
      if (site<0): #track doesn't start till later times
        continue
      data = np.load(fTrack)
      sites0 = list(data['sites0']);
      data.close()
      iBasin = sites0.index(site)
      
      #grab the metrics at this time for this basin
      data = np.load(fMetric)
      vals = [data[key][iBasin] for key in metricNames]
      data.close()
      #valsStr = np.array_str(vals, max_line_width=1000000)
      valsStr = str(vals)[1:-1] #leave out [ and ]
      s += valsStr+'\n'
    #end iTime  
    s += '-1\n'
    f.write(s); s = ''
  #end iTrack
    
  #f.write(s)
  
  f.close()

def demo_formTracks():
  fDir = '/data02/cases/summer2006/eraI/'
  filesMetr = sorted(glob.glob(fDir+'/seg/fields_*.npz'), key=os.path.getmtime)
  filesSeg = sorted(glob.glob(fDir+'/seg/seg_*.npz'), key=os.path.getmtime)
  filesMetrics = sorted(glob.glob(fDir+'/seg/metrics_*.npz'), key=os.path.getmtime)
  filesTrack = sorted(glob.glob(fDir+'/track/track_seg*.npz'), key=os.path.getmtime)
  filesTrack = filesTrack[0:20]
  
  #site0 = 13771 #for cfsr 2006/07/20 12Z
  #site0 = 5339 #for ERA-I 2006/07/20 0Z
  iTime = 0
  fTrack = filesTrack[iTime]
  data = np.load(fTrack)
  sites0 = data['sites0'];
  data.close()
  
  #gather/split individual tracks over time ---------
  basinTrack = []
  for site0 in sites0:
    basinTrack += gather_1To1Track(filesTrack, site0)
  print "nTracks: ", len(basinTrack)
  
  if (False):
    #print some tpv metric
    for trackSeq in basinTrack:
      for iTime in xrange(len(trackSeq)):
        site = trackSeq[iTime]
        val = calc_tpvMetric(iTime, site, filesMetr, filesSeg)
        print val
  
  if (True):
    #write tracks to file
    fSave = fDir+'tracks.txt'
    print "Writing tracks to: "+fSave
    writeFile_tracks(fSave, basinTrack, filesTrack, filesMetrics)
  
  if (True):
    #plot some output
    fMesh = filesMetr[0]
    dataMesh = np.load(fMesh)
    lat = dataMesh['lat']; lon = dataMesh['lon']; #nLat = len(lat); nLon = len(lon)
    dataMesh.close()
    
    plt.figure()
    m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
    m.drawcoastlines()
    m.drawmapboundary()
    for trackSeq in basinTrack:
      plot_1Track(np.array(trackSeq,dtype=int), lat, lon, m)
    plt.show()

def demo_formTracks_time():
  #form and write tracks for basins found at each time
  
  fDir = '/data02/cases/summer2006/eraI/'
  filesMetr = sorted(glob.glob(fDir+'/seg/fields_*.npz'), key=os.path.getmtime)
  filesSeg = sorted(glob.glob(fDir+'/seg/seg_*.npz'), key=os.path.getmtime)
  filesMetrics = sorted(glob.glob(fDir+'/seg/metrics_*.npz'), key=os.path.getmtime)
  filesTrack = sorted(glob.glob(fDir+'/track/track_seg*.npz'), key=os.path.getmtime)
  filesTrack = filesTrack[0:20]
  nTimes = len(filesTrack); #print nTimes
    
  basinTrack = []
  for iTime in xrange(0,nTimes-1):
    fTrack = filesTrack[iTime]
    data = np.load(fTrack)
    sites0 = data['sites0'];
    data.close()
    
    if (iTime>0): #for checking if basin(iTime) was from some basin(iTime-1)
      fTrackPrev = filesTrack[iTime-1]
      data = np.load(fTrackPrev)
      sites0Prev = data['sites0']; sites1Prev = data['sites1']; fromBasin = data['isMatch'].transpose()
      data.close()
    
    #gather/split individual tracks over time ---------
    nSkip = 0
    for site0 in sites0:
      basinsPrev = []
      if (iTime>0):
        basinsPrev = get_correspondingSites(site0, sites1Prev, sites0Prev, fromBasin)
      if (len(basinsPrev)>0): #skip duplicating tracks of basins already found
        nSkip = nSkip+1
        continue
      newTracks = gather_1To1Track(filesTrack[iTime:], site0)
      for trackSeq in newTracks:
        if (len(newTracks)<=1):
          continue
        basinTrack.append([-1]*iTime+trackSeq) #pad previous times with <0's
    print "{0}/{2} basins skipped at timeInd {1} since correspondence at previous time".format(nSkip, iTime, len(sites0))
  #end iTime    
  print "nTracks: ", len(basinTrack)
  print basinTrack[-1]
    
  if (False):
    #write tracks to file
    fSave = fDir+'tracks.txt'
    print "Writing tracks to: "+fSave
    writeFile_tracks(fSave, basinTrack, filesTrack, filesMetrics)
    print "Done writing tracks"
    
  if (True):
    #plot some output
    fMesh = filesMetr[0]
    dataMesh = np.load(fMesh)
    lat = dataMesh['lat']; lon = dataMesh['lon']; #nLat = len(lat); nLon = len(lon)
    dataMesh.close()
    
    plt.figure()
    m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
    m.drawcoastlines()
    m.drawmapboundary()
    for trackSeq in basinTrack:
    #for trackSeq in basinTrack[0:1000]:
      #print trackSeq
      seq = [i for i in trackSeq if i>=0]
      seq = np.array(seq,dtype=int)
      plot_1Track(seq, lat, lon, m)
    plt.show()    
  
def calc_tpvMetric(iTime, site, filesMetr, filesSeg):
  #wrapper for calculating some tpv value
  
  fMetr = filesMetr[iTime]
  fSeg = filesSeg[iTime]
  
  data = np.load(fMetr)
  lat = data['lat']; lon = data['lon']; nLat = len(lat); nLon = len(lon)
  u = data['u']; v = data['v']; thetaFlat = data['theta']
  inRegion = data['inRegion']
  #theta = segment_ll_flat.unflatten_1dTo2d(thetaFlat, nLat, nLon)
  data.close()
  
  data_seg = np.load(fSeg)
  cell2Site = data_seg['cell2Site'][:]
  cellIsMin = data_seg['cellIsMin'][:]
  cellIsMax = data_seg['cellIsMax'][:]
  data_seg.close()
  
  amp = basinMetrics.calc_amplitude(site, cell2Site, thetaFlat)
  
  return amp
  
def plot_1Track(flatInds, lat, lon, m):
  
  nLon = len(lon)
  iLat, iLon = segment_ll_flat.index_1dTo2d(flatInds, nLon)
  r2d = 180./np.pi;
  lats = lat[iLat]*r2d; lons = lon[iLon]*r2d
  
  x, y = m(lons, lats)
  m.plot(x, y,'b-',linewidth=.5)
  #m.scatter(x[::8], y[::8], marker='o')
  m.scatter(x[0], y[0], s=10, c='g', marker='^')
  m.scatter(x[-1], y[-1], s=10, c='r', marker='v')
  #plt.show()
  
if __name__=='__main__':
  #demo_formTracks()
  demo_formTracks_time()
  
  
