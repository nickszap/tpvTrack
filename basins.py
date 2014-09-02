import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
import glob
import datetime as dt
import os
#from scipy import stats

import segment_ll_flat
import track_ll as track
import sys; sys.path.append("/home/nickszap/Dropbox/pythonScripts/minAreaBoundingRectangle/")
import bbox_test

'''
#We can characterize the shapes of TPVs in various ways.
Given discrete basins in cells, we have >= info needed.
#For measures of asymmetry:
-eccentricity (fit ellipse, "spindles" to every boundary point,...)
-skewness about extremum
-moment of inertia about extremum
-distance of extremum to boundary 
  (can do spherical distance to boundary or 
  time=distance/speed with sign as vel.distance)

The boundary is cell inside w/ neighbor outside.
-arc length of boundary
-great circle distance vs distance w/in object

"Volume" measures:
-1/A \int(theta-theta0 dA)
'''

#circulation, amplitude, area, moment of inertia
#solidity: area(basin)/area(convexHull) but apparently convex hull of region on a sphere
#isn't straightforward
r = 6370.e3
r2d = 180./np.pi

def calc_aspectRatio(site, cell2Site, lat, lon):
  #project all points onto a plane. Then use properties of a bounding box fit around cells.
  
  if (np.max(lat)>20):
    print "Uhoh. Check that lat, lon input in radians"
  
  nLat = len(lat); nLon = len(lon)
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  latBasin = lat[iLat]; lonBasin = lon[iLon]; nPts = len(latBasin)
  
  lat0 = np.median(latBasin)*r2d; lon0 = np.median(lonBasin)*r2d
  m = Basemap(projection='ortho',lon_0=lon0,lat_0=lat0, resolution='l')
  x, y = m(lonBasin*r2d, latBasin*r2d);
  
  if (False):
    plt.figure()
    m.drawcoastlines()
    m.drawmapboundary()
    m.scatter(x, y, s=1, marker='o')
    plt.show()
  
  xyPoints = np.empty((nPts,2), dtype=float)
  xyPoints[:,0] = x[:]; xyPoints[:,1] = y[:]
  
  maxL, minL, rotAngle = bbox_test.runCase(xyPoints, False)
  print 'min/max radii of bounding box: ', minL/2.e3, maxL/2.e3
  return maxL/minL
  
def calc_majorAsymmetry(site, cell2Site, lat, lon, areaLatCell):
  #let the longest dimension of the bounding box be the major axis.
  #then have areas (+ and -) to both sides of extremum.
  #return ratio
  
  if (np.max(lat)>20):
    print "Uhoh. Check that lat, lon input in radians"
  
  nLat = len(lat); nLon = len(lon)
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  latBasin = lat[iLat]; lonBasin = lon[iLon]; nPts = len(latBasin)
  
  lat0 = np.median(latBasin)*r2d; lon0 = np.median(lonBasin)*r2d
  m = Basemap(projection='ortho',lon_0=lon0,lat_0=lat0, resolution='l')
  x, y = m(lonBasin*r2d, latBasin*r2d);
  
  xyPoints = np.empty((nPts,2), dtype=float)
  xyPoints[:,0] = x[:]; xyPoints[:,1] = y[:]
  maxL, minL, rotAngle = bbox_test.runCase(xyPoints, False)
  
  #project areas onto line centered at basin extremum
  iLat0, iLon0 = segment_ll_flat.index_1dTo2d(site, nLon)
  x0, y0 = m(lon[iLon0]*r2d, lat[iLat0]*r2d);
  dxLine = np.cos(rotAngle); dyLine = np.sin(rotAngle)
  #just need sign of dot product
  d = (x-x0)*dxLine+(y-y0)*dyLine
  sgnD = np.sign(d)
  
  if (True):
    plt.figure()
    plt.scatter(x,y)
    meanD = np.max(x-x0); plt.plot([x0,x0+dxLine*meanD], [y0,y0+dyLine*meanD], 'r')
    #m.drawcoastlines()
    #m.drawmapboundary()
    #m.scatter(x, y, s=1, marker='o')
    plt.show()
  
  areaCells = areaLatCell[iLat]/1.e6
  areaNeg = np.sum(areaCells[sgnD<0]); areaPos = np.sum(areaCells[sgnD>0])
  return max(areaNeg,areaPos)/min(areaNeg, areaPos)
  
def calc_circulation(site, cell2Site, vort, areaLatCell, nLat, nLon):
  #\int gradxu . dA
  
  #lat/lon cells in each basin
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  
  vortCells = vort[iLat, iLon]
  cellAreas = areaLatCell[iLat]
  return np.dot(vortCells, cellAreas)

def calc_secondMomentArea(site, cell2Site, areaLatCell, lat, lon):  
  #moment of inertia about extremum: \int r^2 dA
  #if want rotational moment, need some form of mass for sum(m r^2)
  #return in units of km^4
  
  nLat = len(lat); nLon = len(lon)
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  
  #moment about extremum
  iLat0, iLon0 = segment_ll_flat.index_1dTo2d(site, nLon)
  lat0 = lat[iLat0]; lon0 = lon[iLon0]
  
  d = segment_ll_flat.calc_distSphere_multiple(r/1.e3, lat0, lon0, lat[iLat], lon[iLon])
  #d /= 1.e3
  cellAreas = areaLatCell[iLat]/1.e6
  J = np.dot(d*d, cellAreas)
  return J
  
def calc_amplitude(site, cell2Site, thetaFlat):
  #max-min value
  
  inBasin = cell2Site==site;
  valsBasin = thetaFlat[inBasin]
  minVal = np.min(valsBasin); maxVal = np.max(valsBasin)
  return maxVal-minVal

def calc_fieldVolume(site, cell2Site, areaLatCell, thetaFlat, valRef, nLat, nLon):  
  #"volume": \int (theta-thetaRef) dA
  #if want rotational moment, need some form of mass for sum(m r^2)
  #return in units of km^4
  
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  
  #moment about extremum
  vals = thetaFlat[inBasin]
  #valRef = np.max(vals)
  vals = valRef-vals
  cellAreas = areaLatCell[iLat]/1.e6
  J = np.dot(vals, cellAreas)
  return J

def calc_amplitudeBoundary(site, cell2Site, nBoundaryFlat, thetaFlat):
  #<boundary>-extremum value
  
  vals = calc_boundaryValues(site, cell2Site, nBoundaryFlat, thetaFlat)
  valb = np.median(vals) #could do area weighted mean, mode, min,...
  #valb = stats.mode(vals)
  val0 = thetaFlat[site]
  return valb-val0
  
def calc_area(site, cell2Site, areaLatCell, nLat, nLon):
  
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  
  return np.sum(areaLatCell[iLat])

def calc_boundaryLength(site, cell2Site, nBoundaryFlat, dLatCell, nLat, nLon):
  
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  inBasin *= nBoundaryFlat>0
  
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  #nBoundaryBasinCell = nBoundaryFlat[inBasin]
  return np.sum(dLatCell[iLat])
  #return np.dot(nBoundaryBasinCell, dLatCell[iLat])

def calc_minMaxDistToBoundary(site, cell2Site, nBoundaryFlat, lat, lon):
  
  nLat = len(lat); nLon = len(lon)
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  inBasin *= nBoundaryFlat>0
  
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  iLat0, iLon0 = segment_ll_flat.index_1dTo2d(site, nLon)
  
  latb = lat[iLat]; lonb = lon[iLon]
  lat0 = lat[iLat0]; lon0 = lon[iLon0]
  
  d = segment_ll_flat.calc_distSphere_multiple(r, lat0, lon0, latb, lonb)
  minD = np.min(d); maxD = np.max(d)
  
  return (minD, maxD)
  
def calc_boundaryValues(site, cell2Site, nBoundaryFlat, thetaFlat):
  
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  inBasin *= nBoundaryFlat>0
  return thetaFlat[inBasin]

def find_boundaryCells(cell2Site, nLat, nLon, inRegion):
  #return isBoundary[lat,lon] >0  if cell2Site[nbr] != cell2Site[self]
  
  isBoundary = np.zeros((nLat, nLon), dtype=int)
  for iLat0 in xrange(nLat):
    for iLon0 in xrange(nLon):        
      ind0 = segment_ll_flat.index_2dTo1d(iLat0, iLon0, nLon)
      if (inRegion[ind0]<=0):
        continue;
      val0 = cell2Site[ind0]
      
      indNbrs = segment_ll_flat.nbrInds_ll_flat(iLat0, iLon0, nLat, nLon)
      valNbrs = cell2Site[indNbrs]
      
      nDiff = np.sum(valNbrs!=val0)
      isBoundary[iLat0, iLon0] = nDiff
      
  return isBoundary

def find_boundaryCells_basin(site, cell2Site, nLat, nLon):
  #return lat/lon indices of cells in basin with cell2Site[nbr] != cell2Site[self]
  
  nCells = nLat*nLon
  inBasin = cell2Site==site
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  isBoundary = np.zeros((nLat, nLon), dtype=int)
  for iCell in indsBasin:
    iLat0, iLon0 = segment_ll_flat.index_1dTo2d(iCell, nLon)
    indNbrs = segment_ll_flat.nbrInds_ll_flat(iLat0, iLon0, nLat, nLon)
    valNbrs = cell2Site[indNbrs]
    
    nDiff = np.sum(valNbrs!=site)
    isBoundary[iLat0, iLon0] = nDiff
  
  return isBoundary

def calc_distToSites(site0, allSites, lat, lon):
  
  nLat = len(lat); nLon = len(lon)
  
  iLat0, iLon0 = segment_ll_flat.index_1dTo2d(site0, nLon)
  iLat, iLon = segment_ll_flat.index_1dTo2d(allSites[allSites!=site0], nLon)
  
  latb = lat[iLat]; lonb = lon[iLon]
  lat0 = lat[iLat0]; lon0 = lon[iLon0]
  
  d = segment_ll_flat.calc_distSphere_multiple(r, lat0, lon0, latb, lonb)
  return d

def demo(fMetr, fSeg):
  
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
  
  sitesMin = cell2Site[cellIsMin>0]; sitesMax = cell2Site[cellIsMax>0];
  sites = sitesMin; #sites = np.append(sitesMin, sitesMax)
  nSites = len(sites)
  
  areaLatCell = track.calc_areaLatStrips(lat, r)/nLon
  vort = segment_ll_flat.calc_vertVorticity_ll(u, v, nLat, nLon, lat, r)
  
  if (False):
    #basin (areal) measures ---------------------------------------
    vals = np.empty(nSites, dtype=float)
    for iSite, site in enumerate(sites):
      circ = calc_circulation(site, cell2Site, vort, areaLatCell, nLat, nLon);
      print '\ncirculation for basin {0}: {1} km^2/s'.format(site, circ/1.e6)

      moment2 = calc_secondMomentArea(site, cell2Site, areaLatCell, lat, lon)
      print '2nd moment of area for basin {0}: {1:e} km^4'.format(site, moment2)

      amp = calc_amplitude(site, cell2Site, thetaFlat)
      print 'theta amplitude for basin {0}: {1} K'.format(site, amp)

      area = calc_area(site, cell2Site, areaLatCell, nLat, nLon)
      #print 'Area for basin {0}: {1} km^2\n'.format(site, area/1.e6)
      rEquiv = np.sqrt(area/np.pi) #pi r^2 = A
      print 'rEquiv for basin {0}: {1} km\n'.format(site, rEquiv/1.e3)

      vals[iSite] = rEquiv/1.e3

    if (True):
      plt.hist(vals, bins=10)
      plt.show()
  #end areal measures -----------------------------
  
  if (True):
    #basin boundary measures -----------------
    nBoundary = find_boundaryCells(cell2Site, nLat, nLon, inRegion)
    #segment_ll_flat.plot_var_ll(lat, lon, isBoundary, [0, 1])
    nBoundaryFlat = segment_ll_flat.flatten_2dTo1d(nBoundary, nLat, nLon)
    #isBoundary = nBoundary>0; 
    #isBoundaryFlat = segment_ll_flat.flatten_2dTo1d(isBoundary, nLat, nLon)
    
    #the length of the boundary is really the polygon border.
    #as a shortcut, we'll just do a distance for a boundary cell
    dLatCell = np.sqrt(areaLatCell)
    
    vals = np.empty(nSites, dtype=float)
    for iSite, site in enumerate(sites):
      dBoundary = calc_boundaryLength(site, cell2Site, nBoundaryFlat, dLatCell, nLat, nLon)
      area = calc_area(site, cell2Site, areaLatCell, nLat, nLon)
      rEquiv = np.sqrt(area/np.pi) #pi r^2 = A
      ratio = dBoundary/(2.*np.pi*rEquiv)
      vals[iSite] = ratio
      
      amp = calc_amplitudeBoundary(site, cell2Site, nBoundaryFlat, thetaFlat)   
      print 'Amplitude from boundary for basin {0}: {1} K'.format(site, amp)
      
    print vals
    plt.hist(vals,bins=10); plt.show()
    
  #end basin boundary measures -----------------

def driver_basinMetrics(filesSeg, filesMetr, filesTrack, sites):
  #pass in 0 or nSites filesSeg and filesMetr.
  #calculate basin properties. return array of 1 property
  
  nSeg = len(filesSeg); nMetr = len(filesMetr); nSites = len(sites)
  nTrack = len(filesTrack)
  vals = []
  for iSite in xrange(nSites):
    #id basin, mesh, and metr info
    ifSeg = iSite%nSeg; ifMetr = iSite%nMetr; ifTrack = iSite%nTrack
    fSeg = filesSeg[ifSeg]; fMetr = filesMetr[ifMetr]; site0 = sites[iSite]
    fTrack = filesTrack[ifTrack]
    
    #gather mesh+metr info
    data = np.load(fMetr)
    lat = data['lat']; lon = data['lon']; nLat = len(lat); nLon = len(lon)
    u = data['u']; v = data['v']; thetaFlat = data['theta']
    data.close()
    
    #gather seg info
    data_seg = np.load(fSeg)
    cell2Site = data_seg['cell2Site'][:]
    data_seg.close()
    
    #gather all sites at this time
    data = np.load(fTrack)
    sitesAtTime = data['sites0'][:]
    data.close()
    
    #start calculating fields
    areaLatCell = track.calc_areaLatStrips(lat, r)/nLon
    vort = segment_ll_flat.calc_vertVorticity_ll(u, v, nLat, nLon, lat, r)
    
    print '\n'+fMetr +'\n'+ fSeg +'\n'+fTrack
    
    #areal measures -------------------------------
    circ = calc_circulation(site0, cell2Site, vort, areaLatCell, nLat, nLon);
    print 'circulation for basin {0}: {1} km^2/s'.format(site0, circ/1.e6)

    moment2 = calc_secondMomentArea(site0, cell2Site, areaLatCell, lat, lon)
    print '2nd moment of area for basin {0}: {1:e} km^4'.format(site0, moment2)

    amp = calc_amplitude(site0, cell2Site, thetaFlat)
    print 'theta amplitude for basin {0}: {1} K'.format(site0, amp)
    
    area = calc_area(site0, cell2Site, areaLatCell, nLat, nLon)
    #print 'Area for basin {0}: {1} km^2\n'.format(site, area/1.e6)
    rEquiv = np.sqrt(area/np.pi) #pi r^2 = A
    print 'rEquiv for basin {0}: {1} km'.format(site0, rEquiv/1.e3)
    print 'average vorticity for basin {0}: {1:e} 1/s'.format(site0, circ/area)
    
    #aspectRatio = calc_aspectRatio(site0, cell2Site, lat, lon)
    #print 'max/min bbox length for basin {0}: {1}'.format(site0, aspectRatio)
    
    asymmLen = calc_majorAsymmetry(site0, cell2Site, lat, lon, areaLatCell)
    print 'Major assymetry: {0}'.format(asymmLen)
    
    #boundary measures ---------------------------------------------
    nBoundary = find_boundaryCells_basin(site0, cell2Site, nLat, nLon)
    nBoundaryFlat = segment_ll_flat.flatten_2dTo1d(nBoundary, nLat, nLon)
    
    minD, maxD = calc_minMaxDistToBoundary(site0, cell2Site, nBoundaryFlat, lat, lon)
    print 'min/max distance for basin {0} to boundary: {1},{2} km'.format(site0, minD/1.e3, maxD/1.e3)
    
    ampBoundary = calc_amplitudeBoundary(site0, cell2Site, nBoundaryFlat, thetaFlat)
    print 'median PT amplitude for basin {0}: {1} K'.format(site0, ampBoundary)
    
    #metrics for extremum environment -----------------------------
    d = calc_distToSites(site0, sitesAtTime, lat, lon)
    minSiteD = np.min(d)
    print 'Min distance to another site for site {0}: {1} km'.format(site0, minSiteD/1.e3)
    
    thetaVol = calc_fieldVolume(site0, cell2Site, areaLatCell, thetaFlat, thetaFlat[site0], nLat, nLon)
    print 'theta volume for basin {0}: {1:e} K m^2'.format(site0, thetaVol)
    ampArea = thetaVol/(area/1.e6)
    print 'area weighted theta amplitude for basin {0}: {1} K'.format(site0, ampArea)
    
    vals.append(ampArea)
  info = 'Area weighted PT Amplitude, K'
  return (np.array(vals), info)
    
def get_correspondingSites(site0, sites0, sites1, isMatch):
  #isMatch[site0, site1]>0 if basins correspond between frames.
  #return list of sites1 corresponding to site0
  
  ind0 = sites0==site0
  isMatch01 = isMatch[ind0,:].squeeze()>0
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

def basinSimilarity(site0, site1, lat, lon, r):
  #cost function approach where min value over basins is most similar.
  #many options for what to include in cost when comparing to basin at previous time:
  #-distance, distance from advected location, min(theta), ...
  
  nLon = len(lon)
  iLat0, iLon0 = segment_ll_flat.index_1dTo2d(site0, nLon)
  iLat1, iLon1 = segment_ll_flat.index_1dTo2d(site1, nLon)
  
  d = segment_ll_flat.calc_distSphere_multiple(r, lat[iLat0], lon[iLon0], lat[iLat1], lon[iLon1])
  
  return d

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

def demo_trackBasin():
  
  #fDir = '/data02/cases/2006/cfsr_anl/track/'
  #fMesh = '/data02/cases/2006/cfsr_anl/seg/fields_pgbhnl.gdas.2006072012.grb2.trop.nc.npz'
  fDir = '/data02/cases/summer2006/eraI/'
  filesMetr = sorted(glob.glob(fDir+'/seg/fields_*.npz'), key=os.path.getmtime)
  filesSeg = sorted(glob.glob(fDir+'/seg/seg_*.npz'), key=os.path.getmtime)
  filesTrack = sorted(glob.glob(fDir+'/track/track_seg*.npz'), key=os.path.getmtime)
  
  #site0 = 13771 #for cfsr 2006/07/20 12Z
  site0 = 5339 #for ERA-I 2006/07/20 0Z
  
  #gather possible basin correspondences over time ---------
  filesTrack = filesTrack[0:40]; print "Not tracking over full range of files"
  basinCorr = gather_timeCorrespondencesBasin(filesTrack, site0)
  
  #to track an individual tpv, want more than just correspondences.
  #so, we filter choices a la similarity cost function.
  fMesh = filesMetr[0]
  dataMesh = np.load(fMesh)
  lat = dataMesh['lat']; lon = dataMesh['lon']; #nLat = len(lat); nLon = len(lon)
  dataMesh.close()
  #areaLatCell = track.calc_areaLatStrips(lat, r)/nLon
  
  basinTrack = filterTimeCorrespondence(basinCorr, lat, lon, r)
  
  #vals, yInfo = driver_basinMetrics(filesSeg, filesMetr, filesTrack, basinTrack)
  iTime0 = 10
  vals, yInfo = driver_basinMetrics(filesSeg[iTime0:], filesMetr[iTime0:], filesTrack[iTime0:], basinTrack[iTime0:])
  
  tBase_era = dt.datetime(2006, 7, 20, 0); deltaT = dt.timedelta(hours=6)
  dateList_era = [tBase_era+i*deltaT for i in xrange(len(basinTrack))]
  plotTimeSeries(dateList_era, vals, yInfo)
  
  plot_track(np.array(basinTrack), lat, lon)

def trackBasin(fMesh, filesTrack, site0):
  #return array of the most similar basin track over time
  
  #r = 6370.e3
  
  #gather possible basin correspondences over time ---------
  basinCorr = gather_timeCorrespondencesBasin(filesTrack, site0)
  
  #to track an individual tpv, want more than just correspondences.
  #so, we filter choices a la similarity cost function.
  dataMesh = np.load(fMesh)
  lat = dataMesh['lat']; lon = dataMesh['lon']; #nLat = len(lat); nLon = len(lon)
  dataMesh.close()
  
  basinTrack = filterTimeCorrespondence(basinCorr, lat, lon, r)
  return np.array(basinTrack)

def demo_compareTPVs():
  #compare some tpv property between cases
  
  #cases info -----------------------
  deltaT = dt.timedelta(hours=6)
  
  fDir_cfsr= '/data02/cases/2006/cfsr_anl/track/'
  fMesh_cfsr = '/data02/cases/2006/cfsr_anl/seg/fields_pgbhnl.gdas.2006072012.grb2.trop.nc.npz'
  filesTrack_cfsr = sorted(glob.glob(fDir_cfsr+'track_seg*.npz'),key=os.path.getmtime)
  filesSeg_cfsr = sorted(glob.glob('/data02/cases/2006/cfsr_anl/seg/seg_*.npz'), key=os.path.getmtime)
  filesMetr_cfsr = sorted(glob.glob('/data02/cases/2006/cfsr_anl/seg/fields_*.npz'), key=os.path.getmtime)
  site0_cfsr = 13771 #for cfsr 2006/07/20 12Z
  tBase_cfsr = dt.datetime(2006, 7, 20, 12)
  
  fDir_era = '/data02/cases/2006/eraI/track/'
  fMesh_era = '/data02/cases/2006/eraI/seg/fields_2006-07-20_00.npz'
  filesTrack_era = sorted(glob.glob(fDir_era+'track_seg*.npz'), key=os.path.getmtime)
  filesSeg_era = sorted(glob.glob('/data02/cases/2006/eraI/seg/seg_*.npz'), key=os.path.getmtime)
  filesMetr_era = sorted(glob.glob('/data02/cases/2006/eraI/seg/fields*.npz'), key=os.path.getmtime)
  site0_era = 5339 #for ERA-I 2006/07/20 0Z
  tBase_era = dt.datetime(2006, 7, 20, 0)
  
  #tpv tracks ------------------------------
  tpv_cfsr = trackBasin(fMesh_cfsr, filesTrack_cfsr, site0_cfsr);
  tpv_era = trackBasin(fMesh_era, filesTrack_era, site0_era)
  
  dateList_cfsr = [tBase_cfsr+i*deltaT for i in xrange(len(tpv_cfsr))]
  dateList_era = [tBase_era+i*deltaT for i in xrange(len(tpv_era))]
  
  #derive some property about each ---------------------------
  r2d = 180./np.pi;
  
  if (True): #latitude
    dataMesh = np.load(fMesh_cfsr)
    lat = dataMesh['lat']; lon = dataMesh['lon']; #nLat = len(lat); nLon = len(lon)
    dataMesh.close()
    nLon = len(lon)
    iLat, iLon = segment_ll_flat.index_1dTo2d(tpv_cfsr, nLon)
    tpvVals_cfsr = lat[iLat]*r2d

    dataMesh = np.load(fMesh_era)
    lat = dataMesh['lat']; lon = dataMesh['lon']; #nLat = len(lat); nLon = len(lon)
    dataMesh.close()
    nLon = len(lon)
    iLat, iLon = segment_ll_flat.index_1dTo2d(tpv_era, nLon)
    tpvVals_era = lat[iLat]*r2d

    #plot the comparison ---------------------------
    plot_basinCompare(dateList_cfsr, tpvVals_cfsr, 'CFSR', dateList_era, tpvVals_era, 'ERA-I', 'latitude')
    
  if (True): #areal measures
    dataMesh = np.load(fMesh_cfsr)
    lat_cfsr = dataMesh['lat']; lon_cfsr = dataMesh['lon']; 
    nLat_cfsr = len(lat_cfsr); nLon_cfsr = len(lon_cfsr)
    dataMesh.close()
    
    dataMesh = np.load(fMesh_era)
    lat_era = dataMesh['lat']; lon_era = dataMesh['lon']; 
    nLat_era = len(lat_era); nLon_era = len(lon_era)
    dataMesh.close()
    
    areaLatCell_cfsr = track.calc_areaLatStrips(lat_cfsr, r)/nLon_cfsr
    areaLatCell_era = track.calc_areaLatStrips(lat_era, r)/nLon_era
    
    nCfsr = len(tpv_cfsr); nEra = len(tpv_era)
    tpvVals_cfsr = np.empty(nCfsr, dtype=float)
    tpvVals_era = np.empty(nEra, dtype=float)
    
    for iTpv in xrange(nCfsr):
      fSeg = filesSeg_cfsr[iTpv]
      data_seg = np.load(fSeg)
      cell2Site = data_seg['cell2Site'][:]
      data_seg.close()
      
      area = calc_area(tpv_cfsr[iTpv], cell2Site, areaLatCell_cfsr, nLat_cfsr, nLon_cfsr)/1.e6
      rEquiv = np.sqrt(area/np.pi) #pi r^2 = A
      tpvVals_cfsr[iTpv] = rEquiv
      
    for iTpv in xrange(nEra):
      fSeg = filesSeg_era[iTpv]
      data_seg = np.load(fSeg)
      cell2Site = data_seg['cell2Site'][:]
      data_seg.close()
      
      area = calc_area(tpv_era[iTpv], cell2Site, areaLatCell_era, nLat_era, nLon_era)/1.e6
      rEquiv = np.sqrt(area/np.pi) #pi r^2 = A
      tpvVals_era[iTpv] = rEquiv
      
    plot_basinCompare(dateList_cfsr, tpvVals_cfsr, 'CFSR', dateList_era, tpvVals_era, 'ERA-I', 'equivalent radius')
    
  if (True): #theta measures
    dataMesh = np.load(fMesh_cfsr)
    lat_cfsr = dataMesh['lat']; lon_cfsr = dataMesh['lon']; 
    nLat_cfsr = len(lat_cfsr); nLon_cfsr = len(lon_cfsr)
    dataMesh.close()
    
    dataMesh = np.load(fMesh_era)
    lat_era = dataMesh['lat']; lon_era = dataMesh['lon']; 
    nLat_era = len(lat_era); nLon_era = len(lon_era)
    dataMesh.close()
    
    areaLatCell_cfsr = track.calc_areaLatStrips(lat_cfsr, r)/nLon_cfsr
    areaLatCell_era = track.calc_areaLatStrips(lat_era, r)/nLon_era
    
    nCfsr = len(tpv_cfsr); nEra = len(tpv_era)
    tpvVals_cfsr = np.empty(nCfsr, dtype=float)
    tpvVals_era = np.empty(nEra, dtype=float)
    
    for iTpv in xrange(nCfsr):
      fSeg = filesSeg_cfsr[iTpv]
      data_seg = np.load(fSeg)
      cell2Site = data_seg['cell2Site'][:]
      data_seg.close()
      fMetr = filesMetr_cfsr[iTpv]
      data = np.load(fMetr)
      thetaFlat = data['theta']
      data.close()
      
      amp = calc_amplitude(tpv_cfsr[iTpv], cell2Site, thetaFlat)
      tpvVals_cfsr[iTpv] = amp
      
    for iTpv in xrange(nEra):
      fSeg = filesSeg_era[iTpv]
      data_seg = np.load(fSeg)
      cell2Site = data_seg['cell2Site'][:]
      data_seg.close()
      fMetr = filesMetr_era[iTpv]
      data = np.load(fMetr)
      thetaFlat = data['theta']
      data.close()
      
      amp = calc_amplitude(tpv_era[iTpv], cell2Site, thetaFlat)
      tpvVals_era[iTpv] = amp
      
    plot_basinCompare(dateList_cfsr, tpvVals_cfsr, 'CFSR', dateList_era, tpvVals_era, 'ERA-I', 'Potential temperature amplitude, K')
    
def plot_track(flatInds, lat, lon):
  
  nLon = len(lon)
  iLat, iLon = segment_ll_flat.index_1dTo2d(flatInds, nLon)
  r2d = 180./np.pi;
  lats = lat[iLat]*r2d; lons = lon[iLon]*r2d
  
  plt.figure()
  m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  x, y = m(lons, lats)

  m.drawcoastlines()
  m.drawmapboundary()
  
  m.plot(x, y)
  #m.scatter(x[::8], y[::8], marker='o')
  m.scatter(x[0], y[0], s=50, c='g', marker='^')
  m.scatter(x[-1], y[-1], s=50, c='r', marker='v')
  
  plt.show()

def plotTimeSeries(dates0, vals0, yInfo):
  
  plt.figure()
  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
  plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=4))
  plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
  
  plt.plot(dates0, vals0)
  
  plt.xlabel('Date')
  plt.ylabel(yInfo)
  
  plt.gcf().autofmt_xdate()

def plot_basinCompare(dates0, vals0, label0, dates1, vals1, label1, yInfo):
  
  plt.figure() #---------------------------
  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
  plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=4))
  plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
  
  plt.plot(dates0, vals0, 'b', label=label0)
  plt.plot(dates1, vals1, 'r', label=label1)
  
  plt.xlabel('Date')
  plt.ylabel(yInfo)
  plt.legend(loc='best')
  
  plt.gcf().autofmt_xdate()
  
  plt.show()
  
if __name__=='__main__':
  
  '''
  fDir = '/data02/cases/summer2006/cfsr_anl/seg/'
  fMetr = fDir+'fields_pgbhnl.gdas.2006072012.grb2.trop.nc.npz'
  fSeg = fDir+'seg_pgbhnl.gdas.2006072012.grb2.trop.nc.npz'
  
  demo(fMetr, fSeg)
  '''
  
  demo_trackBasin()
  #demo_compareTPVs()
