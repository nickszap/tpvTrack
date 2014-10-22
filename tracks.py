#The track program creates correspondences between basins at neighboring times.
#We want tracks that extend over the lifetime of each/all basins

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap

import basinMetrics
import correspond

r2d = 180./np.pi

def form_track_site(fNameCorr, iTimeStart, iTimeEnd, site0, trackOnlyMajor=False):
  # follow a given site throughout the correspondences "tree" and split tree into individual tracks
  
  tracks_checkContinue = [[site0]]
  trackList = []
  while (len(tracks_checkContinue)>0):
    #can continue a track if: have more times AND site corresponds to future site
    basinTrack = tracks_checkContinue.pop()
    nTimes = len(basinTrack); site0 = basinTrack[-1]
    
    iTime = iTimeStart+nTimes-1 #0-indexing for time
    if (iTime>iTimeEnd): #no more times left
      trackList.append(basinTrack)
      continue
    
    corrSites, corrTypes = correspond.get_correspondingSites(fNameCorr, iTime, site0)
    if (len(corrSites)<1): #don't connect to future site
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
  
def form_tracks_iTime(fNameCorr, iTimeStart, iTimeEnd, sites0):
  trackList = []
  for site in sites0:
    print "Forming tracks from correspondences for initial site {0} at time {1}".format(site, iTimeStart)
    siteTracks = form_track_site(fNameCorr, iTimeStart, iTimeEnd, site, trackOnlyMajor=True)
    trackList.extend(siteTracks)
    
  return trackList
  
def run_tracks(fNameTracks, fCorr, iTimeStart, iTimeEnd, fMetrics=''):
  #find tracks for sites at iTimeStart going to at most iTimeEnd
  
  iTime = iTimeStart
  sites0, corrSites, typeCorr = correspond.read_corr_iTime(fCorr, iTime)
  nSites0 = len(sites0)
  
  #sub-select possible tracks: --------------------
  #-major tracks
  #-tracks that start at this time
  #-tracks that started before this time
  
  notInPrev = np.ones(nSites0,dtype=int)
  
  if (iTime>0):
    sites0Prev, corrSites, typeCorr = correspond.read_corr_iTime(fCorr, iTime-1)
    for iSite in xrange(nSites0):
      for basins in corrSites:
        if sites0[iSite] in basins:
          notInPrev[iSite] = 0
          break
    #
  print "{0}/{1} sites started at time {2}".format(np.sum(notInPrev>0), nSites0, iTime) 
  sites0 = sites0[notInPrev>0]
  
  trackList = form_tracks_iTime(fCorr, iTimeStart, iTimeEnd, sites0)
  
  if (fMetrics==''):
    write_tracks_cells(fNameTracks, trackList)
  else:
    dataMetrics = netCDF4.Dataset(fMetrics,'r')
    write_tracks_metrics_iTime(fNameTracks, iTimeStart, trackList, dataMetrics)
    dataMetrics.close()

def write_tracks_cells(fNameTracks, trackList):
  print "Appending to file: "+fNameTracks
  f = open(fNameTracks,'a')
  for track in trackList:
    s = ' '.join(str(i) for i in track)
    f.write(s+'\n')
  f.close()
  
def write_tracks_metrics_iTime(fSave, iTime0, trackList, dataMetrics):
  '''
  format is:
  (header) nTracks metric1 metric2 ...
  timeStartTrack1 nTimesTrack1
  (time1) metric1 metric2 ... for track1
  (time2) metric1 metric2 ... for track1
  -1
  timeStartTrack2 nTimesTrack2
  (time1) metric1 metric2 ... for track2
  (time2) metric1 metric2 ... for track2
  -1
  .
  .
  .
  '''
  
  f = open(fSave,'a')
  
  for track in trackList:
    s = ''
    nTimes = len(track);
    s += '{0} {1}\n'.format(iTime0, nTimes)
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

def plot_tracks_cells(fTracks, mesh, fDirSave):
  f = open(fTracks,'r')
  
  m = Basemap(projection='ortho',lon_0=0,lat_0=89.5, resolution='l')
  r2d = 180./np.pi
  
  plt.figure()
  m.drawcoastlines()
  
  for line in f:    
    cellStr = line.strip().split()
    trackList = [int(i) for i in cellStr]
    if (len(trackList)<1):
      continue
    
    lat, lon = mesh.get_latLon_inds(np.array(trackList,dtype=int))
    lat *= r2d; lon *= r2d
    x,y = m(lon,lat)
    
    print trackList
    print lat; print lon
    
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
  
