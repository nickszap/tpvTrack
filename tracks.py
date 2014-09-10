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

def form_1To1Track(basinTrack, site0):
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

def form_track_site(fNameCorr, iTimeStart, iTimeEnd, site0):
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
    
    corrSites = correspond.get_correspondingSites(fNameCorr, iTime, site0)
    if (len(corrSites)<1): #don't connect to future site
      trackList.append(basinTrack)
      continue
      
    #if here, connect to >=1 future sites so can continue track
    for site1 in corrSites:
      #continue by duplicating the entire previous track
      tracks_checkContinue.append(basinTrack+[site1])
    
  return trackList
  
def form_tracks_iTime(fNameCorr, iTimeStart, iTimeEnd, sites0):
  trackList = []
  for site in sites0:
    print "Forming tracks from correspondences for initial site {0} at time {1}".format(site, iTimeStart)
    siteTracks = form_track_site(fNameCorr, iTimeStart, iTimeEnd, site)
    trackList.extend(siteTracks)
    
  return trackList
  
def run_tracks(fNameTracks, fCorr, iTimeStart, iTimeEnd, fMetrics=''):
  
  iTime = iTimeStart
  sites0, corrSites, typeCorr = correspond.read_corr_iTime(fCorr, iTime)
  nSites0 = len(sites0)
  notInPrev = np.ones(nSites0,dtype=int)
  
  #sub-select possible tracks: --------------------
  #-major tracks
  #-tracks that start at this time
  #-tracks that started before this time
  
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


