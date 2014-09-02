import numpy as np
import netCDF4
import glob
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
#from scipy import stats #for mode

import segment_ll_flat

print "latitudes in [-pi/2, pi/2] and longitudes in [0, 2pi)"

def findClosestCell2Pt_diff(latPt, lonPt, lat, lon):
  #closest pt is (closest lat, closest lon)
  iLat = np.argmin(np.abs(lat-latPt)); 
  iLon = np.argmin(np.abs(lon-lonPt));
  
  return (iLat, iLon)

def findClosestCell2Pt_ind(latPt, lonPt, nLat, nLon):
  #given lat [-pi/2, pi/2] and lon [0, 2pi), return 
  #lat, lon indices of specified ll coordinate
  
  #lon = iLon*2pi/nLon
  iLon = int(np.round(lonPt*nLon/(2.*np.pi)))
  if (iLon<0):
    iLon = 0
  if (iLon>nLon-1):
    iLon = nLon-1
    
  #lat = pi/2-iLat*pi/(nLat-1) since inclusive endpts
  iLat = int(np.round((np.pi/2-latPt)*(nLat-1)/np.pi))
  if (iLat<0):
    iLat = 0
  if (iLat>nLat-1):
    iLat = nLat-1
    
  return (iLat, iLon)
  
def demo_plotSeg_from_npz(fName_seg, fName_metr, fSave):
  
  data_seg = np.load(fName_seg)
  cell2Site = data_seg['cell2Site'][:]
  cellIsMin = data_seg['cellIsMin'][:]
  cellIsMax = data_seg['cellIsMax'][:]
  data_seg.close()
  
  #pull in metr fields...should we fill missing?
  metrExt = fName_metr.split('.')[-1]
  if (metrExt=='nc'):
    data = netCDF4.Dataset(fName_metr,'r')
    lat, lon, u, v, tmp, press = segment_ll_flat.get_segmentVars_file(data)
    data.close()

    nLat = len(lat); nLon = len(lon)
    theta = segment_ll_flat.calc_potentialTemperature(tmp, press)
    thetaFlat = segment_ll_flat.flatten_2dTo1d(theta, nLat, nLon)
  elif (metrExt=='npz'):
    data = np.load(fName_metr)
    lat = data['lat']; lon = data['lon']; nLat = len(lat); nLon = len(lon)
    thetaFlat = data['theta']; #u0 = data['u']; v0 = data['v'];
  else:
    print "Unrecognized file extension for fMetr: "+fName_metr
    
  #print thetaFlat[cellIsMin>0]; print thetaFlat[cellIsMax>0]
  isMin = segment_ll_flat.unflatten_1dTo2d(cellIsMin, nLat, nLon)
  isMax = segment_ll_flat.unflatten_1dTo2d(cellIsMax, nLat, nLon)
  vals = thetaFlat[cell2Site];
  vals = segment_ll_flat.unflatten_1dTo2d(vals, nLat, nLon)
  segment_ll_flat.plot_segment_save(fSave, lat, lon, vals, isMin, isMax)

def demo_advectFeature_from_npz(fName_seg, fName_metr, latFs, lonFs):
  #advect a basin id'd by cell2Site[latF,lonF]
  #pass in array of latF, lonF
  
  dt = 3.*3600 #s
  r = 6371.e3
  
  data_seg = np.load(fName_seg)
  cell2Site = data_seg['cell2Site'][:]
  #cellIsMin = data_seg['cellIsMin'][:]
  #cellIsMax = data_seg['cellIsMax'][:]
  data_seg.close()
  
  #pull in metr fields...should we fill missing?
  data = netCDF4.Dataset(fName_metr,'r')
  lat, lon, u, v, tmp, press = segment_ll_flat.get_segmentVars_file(data)
  data.close()
  nLat = len(lat); nLon = len(lon)
  #u = segment_ll_flat.fill_missingVals_region(u, nLat, nLon, isMissing, inRegionHalo)
  #v = segment_ll_flat.fill_missingVals_region(v, nLat, nLon, isMissing, inRegionHalo)
  
  #find basin and advect
  nBasins = len(latFs)
  plt.figure()
  for iBasin in xrange(nBasins):
    latF = latFs[iBasin]; lonF = lonFs[iBasin]
    iLat = np.argmin(np.abs(lat-latF)); iLon = np.argmin(np.abs(lon-lonF));
    ind = segment_ll_flat.index_2dTo1d(iLat, iLon, nLon)

    basinInd = cell2Site[ind]
    newLat, newLon = advect_feature(basinInd, cell2Site, lat, lon, u, v, dt, r)

    #plot it up
    r2d = 180./np.pi
    m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
    #advected points
    x,y = m(newLon*r2d, newLat*r2d)
    #original basin extremum
    iLat, iLon = segment_ll_flat.index_1dTo2d(basinInd, nLon)
    latb = lat[iLat]; lonb = lon[iLon]
    #x0, y0 = m(lonF*r2d, latF*r2d)
    x0, y0 = m(lonb*r2d, latb*r2d)
    #print x.shape, y.shape

    m.drawcoastlines()
    m.drawmapboundary()

    m.scatter(x, y, c='b', marker='x')
    m.scatter(x0, y0, c='r', s=40, marker='o')
    
  plt.show()
  
def advect_LatLon(u, v, lat, lon, dt, r ):
  #return new lat/lon coordinates based on:
  #u,v in m/s, lat/lon in radians, dt in s, rSphere in m
  
  #u = r cos(lat) dLon/dt, v = r dLat/dt
  #constrain latitudes to be in [-pi/2, pi/2] and longitudes in [0, 2pi)
  
  latOut = np.copy(lat); lonOut = np.copy(lon)
  
  dLon_dt = u/(r*np.cos(lat)) #if lat=+- pi/2, this will be big
  dLat_dt = v/r
  
  latOut += dLat_dt*dt; lonOut += dLon_dt*dt
  #want to bound lat/lon. Note that for silly big velocities, the following
  #will adjust to outside of bounds, but then the whole tracking idea is screwed anyway.
  
  #bound lat. imagine crossing north pole so 89->93 is really 89->87N
  crossedN = latOut>np.pi/2
  latOut[crossedN] = np.pi-latOut[crossedN]; #pi/2 - (lat-pi/2)
  lonOut[crossedN] += np.pi #other half of hemisphere
  crossedS = latOut<-np.pi/2
  latOut[crossedS] = -np.pi-latOut[crossedS]; #-pi/2 - (lat-(-pi/2))
  lonOut[crossedS] += np.pi #other half of hemisphere
  #bound lon
  lonOut = lonOut%(2.*np.pi) #[0,2pi)
  
  return (latOut, lonOut)

def advect_feature(siteInd, cell2Site, lat, lon, u, v, dt, r):
  #return advected lat, lon points of feature cells
  
  nCells = len(cell2Site)
  inFeature = cell2Site==siteInd
  cellInds = np.arange(nCells)[inFeature]
  
  latInds, lonInds = segment_ll_flat.index_1dTo2d(cellInds, len(lon))
  latFeature = lat[latInds]; lonFeature = lon[lonInds]
  uFeature = u[latInds, lonInds]; vFeature = v[latInds, lonInds]
  
  newLat, newLon = advect_LatLon(uFeature, vFeature, latFeature, lonFeature, dt, r )
  
  return (newLat, newLon)

def area_latLonCell(latN, latS, dLon, r):
  #input angles in radians of northern bound, southern bound, and width of rectangle in radians
  #r is radius of circle in meters
  
  solidAngle = dLon*(np.sin(latN)-np.sin(latS))
  area = solidAngle*r*r
  return area
  
def calc_areaLatStrips(lat, r):
  #return the areas of the latitude strips centered around the specified points.
  #lat[0] is north and lat[nLat-1] is south.
  #for multiple evenly spaced cells at a given latitude, areaCell = areaLats[iLat]/nLon
  
  #caps of sphere for poles and latitude strips for rest.
  dLon = 2.*np.pi
  nLat = len(lat)
  areaLats = np.empty(nLat, dtype=float)
  nonPole = np.arange(1,nLat-1) #[1,nlat-2]
  latsN = .5*(lat[nonPole]+lat[nonPole-1]) #face is halfway between centers
  latsS = .5*(lat[nonPole]+lat[nonPole+1])
  areaLats[nonPole] = area_latLonCell(latsN, latsS, dLon, r)
  
  #assuming symmetric latitudes for poles, area of cap is just residual.
  #this probably has more roundoff error than just computing area of cap zone, but
  #this doesn't matter too much for weighting since we normalize by sum(weights) anyway.
  nonPoleArea = np.sum(areaLats[nonPole]); totArea = 4.*np.pi*r*r
  #print nonPoleArea, totArea
  areaCap = (totArea-nonPoleArea)/2
  areaLats[0] = areaCap; areaLats[nLat-1] = areaCap
  
  return areaLats

def calc_candidateCells_basin(siteInd, cell2Site, lat, lon, u, v, dt, r):
  #return lat/lon inds of cells that basin advects to next time (+ or - dt)
  #this is effectively a map site_t0 to cell_t1
  
  #coordinates of advected cell centers. 
  #so, we're assuming that mesh is dense enough wrt feature size that cell centers are sufficient?
  #or, overlap is sufficient that this will id it
  latPts, lonPts = advect_feature(siteInd, cell2Site, lat, lon, u, v, dt, r)
  
  #cells corresponding to those coordinates
  nPts = len(latPts); 
  iLat = np.empty(nPts, dtype=int); iLon = np.empty(nPts, dtype=int)
  nLat = len(lat); nLon = len(lon)
  for iPt in xrange(nPts):
    #l0, l1 = findClosestCell2Pt_ind(latPts[iPt], lonPts[iPt], nLat, nLon);
    l0, l1 = findClosestCell2Pt_diff(latPts[iPt], lonPts[iPt], lat, lon)
    iLat[iPt] = l0; iLon[iPt] = l1
  
  #to get unique cells: flatten, unique, then unflatten
  return (iLat, iLon)

def getUnique_llCells(iLat, iLon, nLon):
  #given possibly repeated lat,lon pairs, return lat,lon w/o copies.  
  #this can come up, eg when advecting points around, can fall in multiple copies of same cell
  
  inds1d = segment_ll_flat.index_2dTo1d(iLat, iLon, nLon)
  inds1d = np.unique(inds1d)
  indsLat, indsLon = segment_ll_flat.index_1dTo2d(inds1d, nLon)
  
  return (indsLat, indsLon)

def getCommon_1dInd(inds0, inds1):
  #return common values between lists
  
  return np.intersect1d(inds0, inds1, assume_unique=False)
  
def make_candidateCorrespondence_overlap(siteInds0, cell2Site0, u0, v0,
                            siteInds1, cell2Site1, u1, v1, lat, lon, dt, r):
  #Given fields+basins at t0 and t0+dt,
  #-create candidate matches by overlapping advection
  #in principle, overlap could mean: 
  #-# of common cells, min threshold for common area, min threshold for convex hulls overlapping area,...
  #return matrix [sites0, sites1] with 0 if don't match and 1 if match
  
  #possible siteInds0 = cell2Site0[cellIsMin>0] for tracking mins
  nLon = len(lon); nSites0 = len(siteInds0); nSites1 = len(siteInds1)
  doOverlap = np.zeros((nSites0, nSites1), dtype=int)
  
  #store advection of t1 sites -dt/2
  sites2Cells_t1 = [None]*nSites1
  for iSite1 in xrange(nSites1):
    siteInd = siteInds1[iSite1]
    #advect basin -dt/2
    iLat, iLon = calc_candidateCells_basin(siteInd, cell2Site1, lat, lon, u1, v1, -.5*dt, r)
    sites2Cells_t1[iSite1] = segment_ll_flat.index_2dTo1d(iLat, iLon, nLon)
    sites2Cells_t1[iSite1] = np.unique(sites2Cells_t1[iSite1])
    
  #see which t0 sites advected dt/2 overlap with future sites advected back
  for iSite0 in xrange(nSites0):
    siteInd = siteInds0[iSite0]
    #advect basin +dt/2
    iLat, iLon = calc_candidateCells_basin(siteInd, cell2Site0, lat, lon, u0, v0, .5*dt, r)
    site2Cells_t0 = segment_ll_flat.index_2dTo1d(iLat, iLon, nLon)
    site2Cells_t0 = np.unique(site2Cells_t0)
    
    for iSite1 in xrange(nSites1):
      commonCells = getCommon_1dInd(site2Cells_t0, sites2Cells_t1[iSite1])
      if (len(commonCells)>0):
        doOverlap[iSite0, iSite1] = 1
  
  #print out some quick diagnostics
  print "For overlapping advection", doOverlap
  
  return doOverlap

def make_candidateCorrespondence_fracOverlap(siteInds0, cell2Site0, u0, v0,
                            siteInds1, cell2Site1, u1, v1, lat, lon, dt, r, areaLatCell):
  #Given fields+basins at t0 and t0+dt,
  #-create candidate matches by overlapping advection
  
  #return matrix with fraction of overlap by area(nCellsMatch)/area(nCellsPossible)
  
  #possible siteInds0 = cell2Site0[cellIsMin>0] for tracking mins
  nLon = len(lon); nSites0 = len(siteInds0); nSites1 = len(siteInds1)
  doOverlap = np.zeros((nSites0, nSites1), dtype=float)
  
  #store advection of t1 sites -dt/2
  sites2Cells_t1 = [None]*nSites1
  areaBasin1 = np.empty(nSites1,dtype=float)
  for iSite1 in xrange(nSites1):
    siteInd = siteInds1[iSite1]
    #advect basin -dt/2
    iLat, iLon = calc_candidateCells_basin(siteInd, cell2Site1, lat, lon, u1, v1, -.5*dt, r)
    sites2Cells_t1[iSite1] = segment_ll_flat.index_2dTo1d(iLat, iLon, nLon)
    sites2Cells_t1[iSite1] = np.unique(sites2Cells_t1[iSite1])
    
    iLat, iLon  = segment_ll_flat.index_1dTo2d(sites2Cells_t1[iSite1], nLon)
    areaBasin1[iSite1] = np.sum(areaLatCell[iLat])
    
  #see which t0 sites advected dt/2 overlap with future sites advected back
  for iSite0 in xrange(nSites0):
    siteInd = siteInds0[iSite0]
    #advect basin +dt/2
    iLat, iLon = calc_candidateCells_basin(siteInd, cell2Site0, lat, lon, u0, v0, .5*dt, r)
    site2Cells_t0 = segment_ll_flat.index_2dTo1d(iLat, iLon, nLon)
    site2Cells_t0 = np.unique(site2Cells_t0)
    
    iLat, iLon  = segment_ll_flat.index_1dTo2d(site2Cells_t0, nLon) 
    areaBasin0 = np.sum(areaLatCell[iLat])
    
    for iSite1 in xrange(nSites1):
      #for frac overlap, tough choice for what to use as #possible cells.
      #consider candidate big t0 and small t1.
      #-for min(cellsInBasin), fraction will be high w/ few cells from big needed
      #-for max(cellsInBasin), fraction will be low even if lots of small covered
      commonCells = getCommon_1dInd(site2Cells_t0, sites2Cells_t1[iSite1])
      iLat, iLon = segment_ll_flat.index_1dTo2d(commonCells, nLon)
      areaCommon = np.sum(areaLatCell[iLat])
      
      minBasinArea = min(areaBasin0, areaBasin1[iSite1])
      frac = areaCommon/minBasinArea
      doOverlap[iSite0, iSite1] = frac
  
  #print out some quick diagnostics
  print "For overlapping advection", doOverlap
  
  return doOverlap

#given candidate matches (eg, by advection), can further filter 
#these possible matches by 
#(1) feature properties a la cost function:
#reasonable changes in 
#-thetaMin, circulation, feature velocity, location,...
#over "short" time.
#for feature velocity, can see that dLoc< max(fac*vel*dt, minDist~100km)
#(2) event consistency:
#propagate, split, merge, genesis, lysis should have characteristics
#-propagate: similar intensity and area
#-split: similar sumArea
#-merge:

def filterCorrespondence_thetaMin(siteInd0, cell2Site0, theta0, 
                        siteInd1, cell2Site1, theta1, dt):
  #feature can "intensify" only so quickly a la pv budgets.
  #this isn't true if (1) a super deep, narrow feature (un)tilts or (2) a piece
  #of a basin branches off
  #Inputs: all flat arrays
  #return 1 if reasonable, 0 if not
  
  dtHours = dt/3600
  dtheta_dt_max = 2/dtHours #K/hr
  
  inBasin0 = cell2Site0==siteInd0
  theta_t0 = theta0[inBasin0]; thetaMin0 = np.min(theta_t0)
  
  inBasin1 = cell2Site1==siteInd1
  theta_t1 = theta1[inBasin1]; thetaMin1 = np.min(theta_t1)
  
  dtheta_dt = np.abs((thetaMin1-thetaMin0)/dtHours)
  #print "Calculated and max dtheta/dt: {0}, {1}".format(dtheta_dt, dtheta_dt_max)
  if (dtheta_dt>dtheta_dt_max):
    return 0
  else:
    return 1

def filterCorrespondence_velDistance_object(siteInd0, cell2Site0, u0, v0, 
                                     siteInd1, cell2Site1, u1, v1, 
                                     lat, lon, areaLatCell, r, dt):
  #
  #feature moves at a certain velocity, not each cell advecting independently.
  #Inputs: cell2Site flat. u,v [lat,lon]
  #return 1 if reasonable, 0 if not
  
  minPossibleDist = 120.e3; velFac = 1.
  
  nLat = len(lat); nLon = len(lon)
  nCells = nLat*nLon
  
  #lat/lon cells in each basin
  inBasin0 = cell2Site0==siteInd0
  indsBasin0 = np.arange(nCells, dtype=int)[inBasin0]
  iLat0, iLon0 = segment_ll_flat.index_1dTo2d(indsBasin0, nLon)
  
  inBasin1 = cell2Site1==siteInd1
  indsBasin1 = np.arange(nCells, dtype=int)[inBasin1]
  iLat1, iLon1 = segment_ll_flat.index_1dTo2d(indsBasin1, nLon)
  
  #area weighted mean velocity as substitute for mass weighted
  uMean0 = np.sum(np.abs(u0[iLat0, iLon0])*areaLatCell[iLat0])/np.sum(areaLatCell[iLat0])
  vMean0 = np.sum(np.abs(v0[iLat0, iLon0])*areaLatCell[iLat0])/np.sum(areaLatCell[iLat0])
  #print "u,v mean (m/s): {0}, {1}".format(uMean0, vMean0)
  speed0 = np.sqrt(uMean0*uMean0+vMean0*vMean0)
  
  uMean1 = np.sum(np.abs(u1[iLat1, iLon1])*areaLatCell[iLat1])/np.sum(areaLatCell[iLat1])
  vMean1 = np.sum(np.abs(v1[iLat1, iLon1])*areaLatCell[iLat1])/np.sum(areaLatCell[iLat1])
  #print "u,v mean: {0}, {1}".format(uMean1, vMean1)
  speed1 = np.sqrt(uMean1*uMean1+vMean1*vMean1)
  
  #maximum distance travelled
  dist0 = velFac*speed0*dt
  dist1 = velFac*speed1*dt
  #if net velocity is 0, still want to have an acceptable region
  maxDist = max(dist0, dist1); maxDist = max(maxDist, minPossibleDist)
  
  #many options for feature location.
  #averaging for centroid of feature doesn't work for longitudes that cross 0 and 2pi.
  #modes of lats and lons is 1 alternative but not robust if dense, distant arm rotating.
  '''
  #mode
  iLatBasin0, freq = stats.mode(iLat0); iLatBasin0 = int(iLatBasin0)
  iLonBasin0, freq = stats.mode(iLon0); iLonBasin0 = int(iLonBasin0)
  iLatBasin1, freq = stats.mode(iLat1); iLatBasin1 = int(iLatBasin1)
  iLonBasin1, freq = stats.mode(iLon1); iLonBasin1 = int(iLonBasin1)
  '''
  #extremum loc
  iLatBasin0, iLonBasin0 = segment_ll_flat.index_1dTo2d(siteInd0, nLon)
  iLatBasin1, iLonBasin1 = segment_ll_flat.index_1dTo2d(siteInd1, nLon)
  
  #print iLonBasin1, freq
  dist = segment_ll_flat.calc_distSphere_multiple(r, lat[iLatBasin0], lon[iLonBasin0],
                                                  lat[iLatBasin1], lon[iLonBasin1])
  #
  #print "Actual and max distance between basin modes (km): {0}, {1}".format(dist/1.e3, maxDist/1.e3)
  
  if (dist>maxDist):
    return 0
  else:
    return 1

def filterCorrespondence_velDistance_extremum(siteInd0, u0, v0, 
                                     siteInd1, u1, v1, 
                                     lat, lon, r, dt):
  #
  #feature moves at a certain velocity, not each cell advecting independently.
  #Inputs: cell2Site flat. u,v [lat,lon]
  #return 1 if reasonable, 0 if not
  
  minPossibleDist = 120.e3; velFac = 1.
  nLon = len(lon)
  
  #extremum loc
  iLat0, iLon0 = segment_ll_flat.index_1dTo2d(siteInd0, nLon)
  iLat1, iLon1 = segment_ll_flat.index_1dTo2d(siteInd1, nLon)
  
  #velocities of locations
  uMean0 = u0[iLat0, iLon0]; vMean0 = v0[iLat0, iLon0];
  speed0 = np.sqrt(uMean0*uMean0+vMean0*vMean0)
  uMean1 = u1[iLat0, iLon0]; vMean1 = v1[iLat0, iLon0]
  speed1 = np.sqrt(uMean1*uMean1+vMean1*vMean1)
  
  #farthest distance travelled
  maxSpeed = max(speed0, speed1)
  maxDist = velFac*maxSpeed*dt
  maxDist = max(maxDist, minPossibleDist)  
  
  #print iLonBasin1, freq
  dist = segment_ll_flat.calc_distSphere_multiple(r, lat[iLat0], lon[iLon0],
                                                  lat[iLat1], lon[iLon1])
  #
  #print "Actual and max distance between basin modes (km): {0}, {1}".format(dist/1.e3, maxDist/1.e3)
  
  if (dist>maxDist):
    return 0
  else:
    return 1

def demo_testTrack():
  
  fnames_metr = sorted(glob.glob('/data02/cases/2014/gfs_4_20140101_00*.nc'))
  fnames_seg = sorted(glob.glob('/data02/cases/2014/segment/seg_dFilter_300km/seg_gfs_4_20140101*.npz'))
  
  fDir = '/data02/cases/2014/segment/track/' #to save
  
  dt = 3*3600.; r = 6370.e3
  
  nFiles = len(fnames_seg);
  for iFile in xrange(nFiles-1): #xrange(4,5):
  #for iFile in xrange(7,14):
    fSeg0 = fnames_seg[iFile]; fMetr0 = fnames_metr[iFile]
    fSeg1 = fnames_seg[iFile+1]; fMetr1 = fnames_metr[iFile+1]
    
    #need u,v fields for advection. if for real, fill missing values w/in region too.
    #theta if filtering on min(theta)
    data = netCDF4.Dataset(fMetr0,'r')
    isMissing0 = segment_ll_flat.get_missingCells_file(data)
    lat, lon, u0, v0, tmp0, press0 = segment_ll_flat.get_segmentVars_file(data)
    data.close()
    data = netCDF4.Dataset(fMetr1,'r')
    isMissing1 = segment_ll_flat.get_missingCells_file(data)
    lat, lon, u1, v1, tmp1, press1 = segment_ll_flat.get_segmentVars_file(data)
    data.close()
    
    nLat = len(lat); nLon = len(lon)
    latThreshHalo = 48.*np.pi/180.
    latThresh = 50.*np.pi/180.
    inRegionHalo = np.zeros((nLat,nLon), dtype=int); inRegionHalo[lat>latThreshHalo,:] = 1
    inRegion = np.zeros((nLat,nLon), dtype=int); inRegion[lat>latThresh,:] = 1
    
    u0 = segment_ll_flat.fill_missingVals_region(u0, nLat, nLon, isMissing0, inRegionHalo)
    v0 = segment_ll_flat.fill_missingVals_region(v0, nLat, nLon, isMissing0, inRegionHalo)
    tmp0 = segment_ll_flat.fill_missingVals_region(tmp0, nLat, nLon, isMissing0, inRegionHalo)
    press0 = segment_ll_flat.fill_missingVals_region(press0, nLat, nLon, isMissing0, inRegionHalo)
    
    u1 = segment_ll_flat.fill_missingVals_region(u1, nLat, nLon, isMissing1, inRegionHalo)
    v1 = segment_ll_flat.fill_missingVals_region(v1, nLat, nLon, isMissing1, inRegionHalo)
    tmp1 = segment_ll_flat.fill_missingVals_region(tmp1, nLat, nLon, isMissing1, inRegionHalo)
    press1 = segment_ll_flat.fill_missingVals_region(press1, nLat, nLon, isMissing1, inRegionHalo)
    
    theta0 = segment_ll_flat.calc_potentialTemperature(tmp0, press0)
    theta1 = segment_ll_flat.calc_potentialTemperature(tmp1, press1)
    
    #segmentation into basins at each time
    data_seg = np.load(fSeg0)
    cell2Site0 = data_seg['cell2Site'][:]
    cellIsMin0 = data_seg['cellIsMin'][:]
    cellIsMax0 = data_seg['cellIsMax'][:]
    data_seg.close()
    data_seg = np.load(fSeg1)
    cell2Site1 = data_seg['cell2Site'][:]
    cellIsMin1 = data_seg['cellIsMin'][:]
    cellIsMax1 = data_seg['cellIsMax'][:]
    data_seg.close()
    
    #choose which basins we want to track
    trackMin = True; trackMax=False #pick one or other
    if (trackMin):
      siteInds0 = cell2Site0[cellIsMin0>0]
      siteInds1 = cell2Site1[cellIsMin1>0]
    elif (trackMax):
      siteInds0 = cell2Site0[cellIsMax0>0]
      siteInds1 = cell2Site1[cellIsMax1>0]
    
    areaLatCell = calc_areaLatStrips(lat, r)/nLon
    if(True):
      doOverlapFrac = make_candidateCorrespondence_fracOverlap(siteInds0, cell2Site0, u0, v0,
                              siteInds1, cell2Site1, u1, v1, lat, lon, dt, r, areaLatCell)
      matchFracThresh = .5; print "Fractional area match if > {0}".format(matchFracThresh)
      overlapFracThresh = 0.0; print "Fractional area overlap if > {0}".format(overlapFracThresh)
      doMatch = (doOverlapFrac>matchFracThresh).astype(int)
      doOverlap = (doOverlapFrac>overlapFracThresh).astype(int)
    else:
      doOverlap = make_candidateCorrespondence_overlap(siteInds0, cell2Site0, u0, v0,
                            siteInds1, cell2Site1, u1, v1, lat, lon, dt, r)
      doMatch = 0*doOverlap
    
    print "Number of correspondences from advection: {0}".format(np.sum(doOverlap))
    print "Number of matches from overlapping advection: {0}".format(np.sum(doMatch))
    
    #filter correspondences
    #nLon = len(lon); nLat = len(lat)
    nSites0 = len(siteInds0); nSites1 = len(siteInds1)
    
    theta0 = segment_ll_flat.flatten_2dTo1d(theta0, nLat, nLon)
    theta1 = segment_ll_flat.flatten_2dTo1d(theta1, nLat, nLon)
    
    if (trackMax):
      theta0 = -theta0; theta1 = -theta1
    
    #area overlaps are matches and we can further add on more ----------------
    isMatch = doMatch.copy()
    
    #if (True):
    if (False):
      print "Filtering by thetaMin"
      for iSite0, site0 in enumerate(siteInds0):
        for iSite1, site1 in enumerate(siteInds1):
          if (doOverlap[iSite0,iSite1]>0 and not doMatch[iSite0,iSite1]):
            isMatch[iSite0,iSite1] += filterCorrespondence_thetaMin(site0, cell2Site0, theta0, 
                          site1, cell2Site1, theta1, dt)
      print "Number of correspondences after thetaMin threshold: {0}".format(np.sum(isMatch))
    
    if (True):
    #if (False):
      print "Filtering by velDistance_object"
      for iSite0, site0 in enumerate(siteInds0):
        for iSite1, site1 in enumerate(siteInds1):
          if (doOverlap[iSite0,iSite1]>0 and not doMatch[iSite0,iSite1]):
            isMatch[iSite0,iSite1] += filterCorrespondence_velDistance_object(site0, cell2Site0, u0, v0, 
                                           site1, cell2Site1, u1, v1, 
                                           lat, lon, areaLatCell, r, dt)
      print "Number of correspondences after velDistance_obj threshold: {0}".format(np.sum(isMatch))
    
    #if (True):
    if (False):
      print "Filtering by velDistance_extremum"
      
      areaLatCell = calc_areaLatStrips(lat, r)/nLon
      
      for iSite0, site0 in enumerate(siteInds0):
        for iSite1, site1 in enumerate(siteInds1):
          if (doOverlap[iSite0,iSite1]>0 and not doMatch[iSite0,iSite1]):
            isMatch[iSite0,iSite1] += filterCorrespondence_velDistance_extremum(site0, u0, v0, site1, u1, v1, 
                                                      lat, lon, r, dt)
      print "Number of correspondences after velDistance_extr threshold: {0}".format(np.sum(isMatch))
    
    #plot correspondences
    iLat0, iLon0 = segment_ll_flat.index_1dTo2d(siteInds0, nLon)
    r2d = 180./np.pi
    lat0 = lat[iLat0]*r2d; lon0 = lon[iLon0]*r2d
    iLat1, iLon1 = segment_ll_flat.index_1dTo2d(siteInds1, nLon)
    lat1 = lat[iLat1]*r2d; lon1 = lon[iLon1]*r2d
    
    plt.figure()
    m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
    x0,y0 = m(lon0, lat0); x1,y1 = m(lon1, lat1);

    m.drawcoastlines()
    m.drawmapboundary()
    
    for iSite0 in xrange(nSites0):
      for iSite1 in xrange(nSites1):
        '''
        #blue-correspond, red-not
        col = 'r'
        if (doOverlap[iSite0,iSite1]>0):
          col = 'b'
        m.drawgreatcircle(lon0[iSite0], lat0[iSite0], lon1[iSite1], lat1[iSite1], del_s=100.0, color=col)
        '''
        if (isMatch[iSite0,iSite1]>0):
          col = 'b'
          m.drawgreatcircle(lon0[iSite0], lat0[iSite0], lon1[iSite1], lat1[iSite1], del_s=100.0, color=col)
          m.scatter(x0[iSite0], y0[iSite0], marker='+', color='g', s=40)
          m.scatter(x1[iSite1], y1[iSite1], marker='o', color='r')
    fInfo = fSeg0.split('/')[-1]
    fNameSave = fDir+fInfo+'_4'+'.png'
    plt.savefig(fNameSave, bbox_inches='tight'); plt.close()
    #plt.show()
    
    #save track information
    #matches 0-1

def demo_trackCFSR():
  
  #dirData = '/data02/cases/2006/cfsr_anl/seg/'
  dirData = '/data02/cases/2006/eraI/seg/'
  fnames_metr = sorted(glob.glob(dirData+'fields_*.npz'))
  fnames_seg = sorted(glob.glob(dirData+'seg_*.npz'))
  
  #fDir = '/data02/cases/2006/cfsr_anl/track/'
  fDir = '/data02/cases/2006/eraI/track/'
  
  dt = 6*3600.; r = 6370.e3
  
  nFiles = len(fnames_seg);
  for iFile in xrange(nFiles-1): #xrange(4,5):
  #for iFile in xrange(0,4):
    fSeg0 = fnames_seg[iFile]; fMetr0 = fnames_metr[iFile]
    fSeg1 = fnames_seg[iFile+1]; fMetr1 = fnames_metr[iFile+1]
    
    #need u,v fields for advection. if for real, fill missing values w/in region too.
    #theta if filtering on min(theta)
    data = np.load(fMetr0)
    lat = data['lat']; lon = data['lon']; nLat = len(lat); nLon = len(lon)
    u0 = data['u']; v0 = data['v']; thetaFlat = data['theta']
    theta0 = segment_ll_flat.unflatten_1dTo2d(thetaFlat, nLat, nLon)
    data.close()
    data = np.load(fMetr1)
    u1 = data['u']; v1 = data['v']; thetaFlat = data['theta']
    theta1 = segment_ll_flat.unflatten_1dTo2d(thetaFlat, nLat, nLon)
    data.close()
        
    #segmentation into basins at each time
    data_seg = np.load(fSeg0)
    cell2Site0 = data_seg['cell2Site'][:]
    cellIsMin0 = data_seg['cellIsMin'][:]
    cellIsMax0 = data_seg['cellIsMax'][:]
    data_seg.close()
    data_seg = np.load(fSeg1)
    cell2Site1 = data_seg['cell2Site'][:]
    cellIsMin1 = data_seg['cellIsMin'][:]
    cellIsMax1 = data_seg['cellIsMax'][:]
    data_seg.close()
    
    #choose which basins we want to track.
    #could just as well track all with:
    #siteInds0 = cell2Site0==range(nCells), ie cell goes to self
    #or np.unique[cell2Site[inRegion]]
    trackMin = True; trackMax=False #pick one or other. 
    #something to do with filtering by thetaMin
    #that I don't use (since tpv's merging/breaking off can have any theta).
    if (trackMin):
      siteInds0 = cell2Site0[cellIsMin0>0]
      siteInds1 = cell2Site1[cellIsMin1>0]
    elif (trackMax):
      siteInds0 = cell2Site0[cellIsMax0>0]
      siteInds1 = cell2Site1[cellIsMax1>0]
    
    areaLatCell = calc_areaLatStrips(lat, r)/nLon
    if(True):
      doOverlapFrac = make_candidateCorrespondence_fracOverlap(siteInds0, cell2Site0, u0, v0,
                              siteInds1, cell2Site1, u1, v1, lat, lon, dt, r, areaLatCell)
      matchFracThresh = .35; print "Fractional area match if > {0}".format(matchFracThresh)
      overlapFracThresh = 0.0; print "Fractional area overlap if > {0}".format(overlapFracThresh)
      doMatch = (doOverlapFrac>matchFracThresh).astype(int)
      doOverlap = (doOverlapFrac>overlapFracThresh).astype(int)
    else:
      doOverlap = make_candidateCorrespondence_overlap(siteInds0, cell2Site0, u0, v0,
                            siteInds1, cell2Site1, u1, v1, lat, lon, dt, r)
      doMatch = 0*doOverlap
    
    print "Number of correspondences from advection: {0}".format(np.sum(doOverlap))
    print "Number of matches from overlapping advection: {0}".format(np.sum(doMatch))
    
    #filter correspondences
    #nLon = len(lon); nLat = len(lat)
    nSites0 = len(siteInds0); nSites1 = len(siteInds1)
    
    theta0 = segment_ll_flat.flatten_2dTo1d(theta0, nLat, nLon)
    theta1 = segment_ll_flat.flatten_2dTo1d(theta1, nLat, nLon)
    
    if (trackMax):
      theta0 = -theta0; theta1 = -theta1
    
    #area overlaps are matches and we can further add on more ----------------
    isMatch = doMatch.copy()
    
    #if (True):
    if (False):
      print "Filtering by thetaMin"
      for iSite0, site0 in enumerate(siteInds0):
        for iSite1, site1 in enumerate(siteInds1):
          if (doOverlap[iSite0,iSite1]>0 and not doMatch[iSite0,iSite1]):
            isMatch[iSite0,iSite1] += filterCorrespondence_thetaMin(site0, cell2Site0, theta0, 
                          site1, cell2Site1, theta1, dt)
      print "Number of correspondences after thetaMin threshold: {0}".format(np.sum(isMatch))
    
    if (True):
    #if (False):
      print "Filtering by velDistance_object"
      for iSite0, site0 in enumerate(siteInds0):
        for iSite1, site1 in enumerate(siteInds1):
          if (doOverlap[iSite0,iSite1]>0 and not doMatch[iSite0,iSite1]):
            isMatch[iSite0,iSite1] += filterCorrespondence_velDistance_object(site0, cell2Site0, u0, v0, 
                                           site1, cell2Site1, u1, v1, 
                                           lat, lon, areaLatCell, r, dt)
      print "Number of correspondences after velDistance_obj threshold: {0}".format(np.sum(isMatch))
    
    #if (True):
    if (False):
      print "Filtering by velDistance_extremum"
      
      areaLatCell = calc_areaLatStrips(lat, r)/nLon
      
      for iSite0, site0 in enumerate(siteInds0):
        for iSite1, site1 in enumerate(siteInds1):
          if (doOverlap[iSite0,iSite1]>0 and not doMatch[iSite0,iSite1]):
            isMatch[iSite0,iSite1] += filterCorrespondence_velDistance_extremum(site0, u0, v0, site1, u1, v1, 
                                                      lat, lon, r, dt)
      print "Number of correspondences after velDistance_extr threshold: {0}".format(np.sum(isMatch))
    
    #plot correspondences
    iLat0, iLon0 = segment_ll_flat.index_1dTo2d(siteInds0, nLon)
    r2d = 180./np.pi
    lat0 = lat[iLat0]*r2d; lon0 = lon[iLon0]*r2d
    iLat1, iLon1 = segment_ll_flat.index_1dTo2d(siteInds1, nLon)
    lat1 = lat[iLat1]*r2d; lon1 = lon[iLon1]*r2d
    
    plt.figure()
    m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
    x0,y0 = m(lon0, lat0); x1,y1 = m(lon1, lat1);

    m.drawcoastlines()
    m.drawmapboundary()
    
    for iSite0 in xrange(nSites0):
      for iSite1 in xrange(nSites1):
        '''
        #blue-correspond, red-not
        col = 'r'
        if (doOverlap[iSite0,iSite1]>0):
          col = 'b'
        m.drawgreatcircle(lon0[iSite0], lat0[iSite0], lon1[iSite1], lat1[iSite1], del_s=100.0, color=col)
        '''
        if (isMatch[iSite0,iSite1]>0):
          col = 'b'
          m.drawgreatcircle(lon0[iSite0], lat0[iSite0], lon1[iSite1], lat1[iSite1], del_s=100.0, color=col)
          m.scatter(x0[iSite0], y0[iSite0], marker='+', color='g', s=40)
          m.scatter(x1[iSite1], y1[iSite1], marker='o', color='r')
    fInfo = fSeg0.split('/')[-1]
    fNameSave = fDir+fInfo+'.png'
    plt.savefig(fNameSave, bbox_inches='tight'); plt.close()
    #plt.show()
    
    #save track information
    #matches 0-1
    fInfo = fSeg0.split('/')[-1]
    fNameSave = 'track_'+fInfo+'.npz'
    f = fDir+fNameSave
    print "Saving track info to file "+f
    np.savez(f, sites0=siteInds0, sites1=siteInds1, isMatch=isMatch)

def demo_trackSteven():
  
  dirData = '/home/scavallo/for_nick/'
  fnames_metr = sorted(glob.glob(dirData+'fields_*.npz'))
  fnames_seg = sorted(glob.glob(dirData+'seg_*.npz'))
  
  fDir = '/data02/cases/test_segment/stevenCase/track/'
  
  dt = 6*3600.; r = 6370.e3
  
  nFiles = len(fnames_seg);
  for iFile in xrange(nFiles-1): #xrange(4,5):
  #for iFile in xrange(0,4):
    fSeg0 = fnames_seg[iFile]; fMetr0 = fnames_metr[iFile]
    fSeg1 = fnames_seg[iFile+1]; fMetr1 = fnames_metr[iFile+1]
    
    #need u,v fields for advection. if for real, fill missing values w/in region too.
    #theta if filtering on min(theta)
    data = np.load(fMetr0)
    lat = data['lat']; lon = data['lon']%(2.*np.pi) #[0,2pi);
    nLat = len(lat); nLon = len(lon)
    u0 = data['u']; v0 = data['v']; thetaFlat = data['theta']
    theta0 = segment_ll_flat.unflatten_1dTo2d(thetaFlat, nLat, nLon)
    data.close()
    data = np.load(fMetr1)
    u1 = data['u']; v1 = data['v']; thetaFlat = data['theta']
    theta1 = segment_ll_flat.unflatten_1dTo2d(thetaFlat, nLat, nLon)
    data.close()
        
    #segmentation into basins at each time
    data_seg = np.load(fSeg0)
    cell2Site0 = data_seg['cell2Site'][:]
    cellIsMin0 = data_seg['cellIsMin'][:]
    cellIsMax0 = data_seg['cellIsMax'][:]
    data_seg.close()
    data_seg = np.load(fSeg1)
    cell2Site1 = data_seg['cell2Site'][:]
    cellIsMin1 = data_seg['cellIsMin'][:]
    cellIsMax1 = data_seg['cellIsMax'][:]
    data_seg.close()
    
    #choose which basins we want to track.
    #could just as well track all with:
    #siteInds0 = cell2Site0==range(nCells), ie cell goes to self
    #or np.unique[cell2Site[inRegion]]
    trackMin = True; trackMax=False #pick one or other. 
    #something to do with filtering by thetaMin
    #that I don't use (since tpv's merging/breaking off can have any theta).
    if (trackMin):
      siteInds0 = cell2Site0[cellIsMin0>0]
      siteInds1 = cell2Site1[cellIsMin1>0]
    elif (trackMax):
      siteInds0 = cell2Site0[cellIsMax0>0]
      siteInds1 = cell2Site1[cellIsMax1>0]
    
    areaLatCell = calc_areaLatStrips(lat, r)/nLon
    if(True):
      doOverlapFrac = make_candidateCorrespondence_fracOverlap(siteInds0, cell2Site0, u0, v0,
                              siteInds1, cell2Site1, u1, v1, lat, lon, dt, r, areaLatCell)
      matchFracThresh = .2; print "Fractional area match if > {0}".format(matchFracThresh)
      overlapFracThresh = 0.0; print "Fractional area overlap if > {0}".format(overlapFracThresh)
      doMatch = (doOverlapFrac>matchFracThresh).astype(int)
      doOverlap = (doOverlapFrac>overlapFracThresh).astype(int)
    else:
      doOverlap = make_candidateCorrespondence_overlap(siteInds0, cell2Site0, u0, v0,
                            siteInds1, cell2Site1, u1, v1, lat, lon, dt, r)
      doMatch = 0*doOverlap
    
    print "Number of correspondences from advection: {0}".format(np.sum(doOverlap))
    print "Number of matches from overlapping advection: {0}".format(np.sum(doMatch))
    
    #filter correspondences
    #nLon = len(lon); nLat = len(lat)
    nSites0 = len(siteInds0); nSites1 = len(siteInds1)
    
    theta0 = segment_ll_flat.flatten_2dTo1d(theta0, nLat, nLon)
    theta1 = segment_ll_flat.flatten_2dTo1d(theta1, nLat, nLon)
    
    if (trackMax):
      theta0 = -theta0; theta1 = -theta1
    
    #area overlaps are matches and we can further add on more ----------------
    isMatch = doMatch.copy()
    
    #if (True):
    if (False):
      print "Filtering by thetaMin"
      for iSite0, site0 in enumerate(siteInds0):
        for iSite1, site1 in enumerate(siteInds1):
          if (doOverlap[iSite0,iSite1]>0 and not doMatch[iSite0,iSite1]):
            isMatch[iSite0,iSite1] += filterCorrespondence_thetaMin(site0, cell2Site0, theta0, 
                          site1, cell2Site1, theta1, dt)
      print "Number of correspondences after thetaMin threshold: {0}".format(np.sum(isMatch))
    
    if (True):
    #if (False):
      print "Filtering by velDistance_object"
      for iSite0, site0 in enumerate(siteInds0):
        for iSite1, site1 in enumerate(siteInds1):
          if (doOverlap[iSite0,iSite1]>0 and not doMatch[iSite0,iSite1]):
            isMatch[iSite0,iSite1] -= filterCorrespondence_velDistance_object(site0, cell2Site0, u0, v0, 
                                           site1, cell2Site1, u1, v1, 
                                           lat, lon, areaLatCell, r, dt)
      print "Number of correspondences after velDistance_obj threshold: {0}".format(np.sum(isMatch))
    
    #if (True):
    if (False):
      print "Filtering by velDistance_extremum"
      
      areaLatCell = calc_areaLatStrips(lat, r)/nLon
      
      for iSite0, site0 in enumerate(siteInds0):
        for iSite1, site1 in enumerate(siteInds1):
          if (doOverlap[iSite0,iSite1]>0 and not doMatch[iSite0,iSite1]):
            isMatch[iSite0,iSite1] += filterCorrespondence_velDistance_extremum(site0, u0, v0, site1, u1, v1, 
                                                      lat, lon, r, dt)
      print "Number of correspondences after velDistance_extr threshold: {0}".format(np.sum(isMatch))
    
    #plot correspondences
    iLat0, iLon0 = segment_ll_flat.index_1dTo2d(siteInds0, nLon)
    r2d = 180./np.pi
    lat0 = lat[iLat0]*r2d; lon0 = lon[iLon0]*r2d
    iLat1, iLon1 = segment_ll_flat.index_1dTo2d(siteInds1, nLon)
    lat1 = lat[iLat1]*r2d; lon1 = lon[iLon1]*r2d
    
    plt.figure()
    m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
    x0,y0 = m(lon0, lat0); x1,y1 = m(lon1, lat1);

    m.drawcoastlines()
    m.drawmapboundary()
    
    for iSite0 in xrange(nSites0):
      for iSite1 in xrange(nSites1):
        '''
        #blue-correspond, red-not
        col = 'r'
        if (doOverlap[iSite0,iSite1]>0):
          col = 'b'
        m.drawgreatcircle(lon0[iSite0], lat0[iSite0], lon1[iSite1], lat1[iSite1], del_s=100.0, color=col)
        '''
        if (isMatch[iSite0,iSite1]>0):
          col = 'b'
          m.drawgreatcircle(lon0[iSite0], lat0[iSite0], lon1[iSite1], lat1[iSite1], del_s=100.0, color=col)
          m.scatter(x0[iSite0], y0[iSite0], marker='+', color='g', s=40)
          m.scatter(x1[iSite1], y1[iSite1], marker='o', color='r')
    fInfo = fSeg0.split('/')[-1]
    fNameSave = fDir+fInfo+'.png'
    plt.savefig(fNameSave, bbox_inches='tight'); plt.close()
    #plt.show()
    
    #save track information
    #matches 0-1
    fInfo = fSeg0.split('/')[-1]
    fNameSave = 'track_'+fInfo+'.npz'
    f = fDir+fNameSave
    print "Saving track info to file "+f
    np.savez(f, sites0=siteInds0, sites1=siteInds1, isMatch=isMatch)

def demo_plotSegHistory():
  
  #fDirSave = '/data02/cases/2014/segment/seg_dFilter_300km/'
  #fnames_metr = sorted(glob.glob('/data02/cases/2014/gfs_4_20140101_00*.nc'))
  #fnames_seg = sorted(glob.glob('/data02/cases/2014/segment/seg_dFilter_300km/seg_gfs_4_20140101*.npz'))
  
  #dirData = '/data02/cases/2006/cfsr_anl/seg/'
  #fDirSave = '/data02/cases/2006/cfsr_anl/basins/' #to save
  dirData = '/data02/cases/2006/eraI/seg/'
  fDirSave = dirData
  fnames_metr = sorted(glob.glob(dirData+'fields_*.npz'))
  fnames_seg = sorted(glob.glob(dirData+'seg_*.npz'))
  
  nFiles = len(fnames_seg)
  print "Saving seg history plots in "+fDirSave
  for iFile in xrange(nFiles):
    fSeg = fnames_seg[iFile]; fMetr = fnames_metr[iFile]
    fSave = 'test_'+fSeg.split('/')[-1]+'.png'
    demo_plotSeg_from_npz(fSeg, fMetr, fDirSave+fSave)  

def demo_plotSegAdvect():
  
  fnames_metr = sorted(glob.glob('/data02/cases/2014/gfs_4_20140101_00*.nc'))
  fnames_seg = sorted(glob.glob('/data02/cases/2014/segment/seg_dFilter_300km/seg_gfs_4_20140101*.npz'))
  
  nFiles = len(fnames_seg); d2r = np.pi/180.
  for iFile in xrange(4,5): #xrange(nFiles):
    fSeg = fnames_seg[iFile]; fMetr = fnames_metr[iFile]
    
    latF = np.array([71., 70.1, 82.5])*d2r; lonF = np.array([130., 360-42.9,360-85])*d2r;
    
    demo_advectFeature_from_npz(fSeg, fMetr, latF, lonF)

def demo_testCandidateCells():
  
  fnames_metr = sorted(glob.glob('/data02/cases/2014/gfs_4_20140101_00*.nc'))
  fnames_seg = sorted(glob.glob('/data02/cases/2014/segment/seg_dFilter_300km/seg_gfs_4_20140101*.npz'))
  
  dt = 3*3600.; r = 6370.e3
  
  nFiles = len(fnames_seg); d2r = np.pi/180.
  for iFile in xrange(1): #xrange(nFiles):
    fSeg = fnames_seg[iFile]; fMetr = fnames_metr[iFile]
    
    data = netCDF4.Dataset(fMetr,'r')
    lat, lon, u, v, tmp, press = segment_ll_flat.get_segmentVars_file(data)
    data.close()
    nLat = len(lat); nLon = len(lon)
    
    data_seg = np.load(fSeg)
    cell2Site = data_seg['cell2Site'][:]
    #cellIsMin = data_seg['cellIsMin'][:]
    #cellIsMax = data_seg['cellIsMax'][:]
    data_seg.close()
    
    #latF = np.array([71., 70.1, 82.5])*d2r; lonF = np.array([130., 360-42.9,360-76.5])*d2r;
    latF = 70.1*d2r; lonF = (360-42.9)*d2r;
    latInd, lonInd = findClosestCell2Pt_ind(latF, lonF, nLat, nLon)
    cellInd = segment_ll_flat.index_2dTo1d(latInd, lonInd, nLon)
    siteInd = cell2Site[cellInd]
  
    latc, lonc = calc_candidateCells_basin(siteInd, cell2Site, lat, lon, u, v, dt, r)
    var = np.zeros((nLat, nLon),dtype=int)
    var[latc, lonc] = 1
    
    segment_ll_flat.plot_var_ll(lat, lon, var, [0, 1])

if __name__ == '__main__':
  #demo_plotSegHistory()
  #demo_plotSegAdvect()
  #demo_testTrack()
  #demo_trackCFSR()
  demo_trackSteven()
  
  
