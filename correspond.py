import numpy as np
import netCDF4
from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.pyplot as plt
import cPickle as pickle; pickleProtocol = 2

import basinMetrics

def advect_LatLon(u, v, latIn, lonIn, dt, r ):
  """
  Return new lat/lon coordinates based on simple advection
  
  Arguments:
  u - zonal velocity
  v - meridional velocity
  latIn - latitude [-pi/2, pi/2]
  lonIn - longitude [0, 2pi)
  r - radius of sphere (m)
  """
  
  #u = r cos(lat) dLon/dt, v = r dLat/dt
  #constrain latitudes to be in [-pi/2, pi/2] and longitudes in [0, 2pi)
  
  dLon_dt = u/(r*np.cos(latIn)) #if lat=+- pi/2, this will be big
  dLat_dt = v/r
  
  lat = latIn+ dLat_dt*dt; lon = lonIn+ dLon_dt*dt
  #want to bound lat/lon. Note that for silly big meridional velocities, the following
  #can adjust to outside of poles, but then the whole tracking idea is screwed anyway.
  
  #bound lat. 
  #imagine crossing north pole so 89->93 is really 89->87N and 180degree switch in longitude
  crossedN = lat>np.pi/2
  lat[crossedN] = np.pi-lat[crossedN]; #pi/2 - (lat-pi/2)
  lon[crossedN] += np.pi #other half of hemisphere
  crossedS = lat<-np.pi/2
  lat[crossedS] = -np.pi-lat[crossedS]; #-pi/2 - (lat-(-pi/2))
  lon[crossedS] += np.pi #other half of hemisphere
  #bound lon
  lon = lon%(2.*np.pi) #[0,2pi)
  
  return (lat, lon)

def advect_feature(siteInd, cell2Site, mesh, u, v, dt):
  """return advected lat, lon points of feature cells"""
  
  nCells = len(cell2Site)
  inFeature = cell2Site==siteInd
  cellInds = np.arange(nCells)[inFeature]
  
  latFeature, lonFeature = mesh.get_latLon_inds(cellInds)
  uFeature = u[cellInds]; vFeature = v[cellInds]
  
  newLat, newLon = advect_LatLon(uFeature, vFeature, latFeature, lonFeature, dt, mesh.r )
  
  if (False):
    #plot the original points and the advected
    print latFeature; print lonFeature;
    print newLat; print newLon
    
    m = Basemap(projection='ortho',lon_0=0,lat_0=89.5, resolution='l')
    r2d = 180./np.pi
    plt.figure()
    m.drawcoastlines()
    x0,y0 = m(lonFeature*r2d, latFeature*r2d)
    x1,y1 = m(newLon*r2d, newLat*r2d)
    m.scatter(x0,y0,color='b')
    m.scatter(x1,y1,color='r')
    plt.show()
  
  return (newLat, newLon)

def advect_basin(siteInd, cell2Site, mesh, u, v, dt):
  """return cells inds of cells that basin advects to next time (+ or - dt)"""
  
  #coordinates of advected cell centers. 
  #so, we're assuming that mesh is dense enough wrt feature size that cell centers are sufficient?
  #or, overlap is sufficient that this will id it
  latPts, lonPts = advect_feature(siteInd, cell2Site, mesh, u, v, dt)
  
  #cells corresponding to those coordinates
  nPts = len(latPts);
  advCells = np.empty(nPts, dtype=int)
  guessCell = siteInd #used for mpas and wrf meshes
  inDomain = [] #used for LAM meshes
  if ('wrf' in mesh.info):
    inDomain = np.ones(nPts,dtype=int)
    
  for iPt in xrange(nPts):
    if (mesh.info == 'mpas'):
      advCells[iPt] = mesh.get_closestCell2Pt(latPts[iPt], lonPts[iPt], guessCell=guessCell)
      guessCell = advCells[iPt]
    elif ('wrf' in mesh.info):
      advCells[iPt] = mesh.get_closestCell2Pt(latPts[iPt], lonPts[iPt], guessCell=guessCell)
      guessCell = advCells[iPt]
      inDomain[iPt] = mesh.isPointInDomain(latPts[iPt], lonPts[iPt], advCells[iPt])
    else:
      advCells[iPt] = mesh.get_closestCell2Pt(latPts[iPt], lonPts[iPt])
  
  if ('wrf' in mesh.info):
    #we'll just ignore points that advect outside the domain
    #print "Number of points advected outside of domain: ", np.sum(inDomain==0)
    advCells = advCells[inDomain>0]
  
  if (False):
    latCell, lonCell = mesh.get_latLon_inds(advCells)
    print latPts; print lonPts
    print latCell; print lonCell
    print advCells
    m = Basemap(projection='ortho',lon_0=0,lat_0=89.5, resolution='l')
    r2d = 180./np.pi
    plt.figure()
    m.drawcoastlines()
    x0,y0 = m(lonPts*r2d, latPts*r2d)
    x1,y1 = m(lonCell*r2d, latCell*r2d)
    m.scatter(x0,y0,color='b')
    m.scatter(x1,y1,color='r')
    plt.show()
  
  #multiple cells can advect into the same cell
  basinCells = np.unique(advCells) #apparently there's a sorting bug fixed in numpy 1.6.2...
  #basinCells = np.array( list(set(advCells)), dtype=int )
  
  return basinCells

def getCommon_1dInd(inds0, inds1):
  """return common values between lists with unique values"""
  return np.intersect1d(inds0, inds1, assume_unique=True)
  
def calc_fracOverlap_advection(sites0, cell2Site0, u0, v0, dt,
                               sites1, cell2Site1, u1, v1, mesh, doMaxArea=False):
  """
  #Given fields+basins at t0 and t0+dt,
  -create candidate matches by overlapping advection
  in principle, overlap could mean: 
  -# of common cells, min threshold for common area, min threshold for convex hulls overlapping area,...
  return matrix with fraction of overlap by area(nCellsMatch)/area(nCellsPossible)
  """
  
  nSites0 = len(sites0); nSites1 = len(sites1)
  fracOverlap = np.zeros((nSites0, nSites1), dtype=float)
  if (doMaxArea):
    fracOverlapMax = np.zeros((nSites0, nSites1), dtype=float)
  
  #store advection of t1 sites -dt/2
  sites2Cells_t1 = [None]*nSites1
  areaBasin1 = np.empty(nSites1,dtype=float)
  for iSite1 in xrange(nSites1):
    siteInd = sites1[iSite1]
    #advect basin -dt/2
    sites2Cells_t1[iSite1] = advect_basin(siteInd, cell2Site1, mesh, u1, v1, -.5*dt)
    areaBasin1[iSite1] = np.sum( mesh.get_area_inds(sites2Cells_t1[iSite1]) )
  #print areaBasin1
    
  #see which t0 sites advected dt/2 overlap with future sites advected back
  for iSite0 in xrange(nSites0):
    siteInd = sites0[iSite0]
    #advect basin +dt/2
    site2Cells_t0 = advect_basin(siteInd, cell2Site0, mesh, u0, v0, .5*dt)
    areaBasin0 = np.sum( mesh.get_area_inds(site2Cells_t0) )
    #print areaBasin0
    
    for iSite1 in xrange(nSites1):
      #for frac overlap, there's a choice for what cells to use.
      #consider candidate big t0 and small t1.
      #-for min(cellsInBasin), fraction will be high w/ few cells from big needed
      #-for max(cellsInBasin), fraction will be low even if lots of small covered
      commonCells = getCommon_1dInd(site2Cells_t0, sites2Cells_t1[iSite1])
      areaCommon = np.sum( mesh.get_area_inds(commonCells) )
      
      #Note that areas can be zero for LAM where point advects out of domain.
      potentialArea = min(areaBasin0, areaBasin1[iSite1])
      #Avoid divide by 0
      frac = 0
      if (potentialArea>0):
        frac = areaCommon/potentialArea
      fracOverlap[iSite0, iSite1] = frac
      if (doMaxArea):
        potentialArea = max(areaBasin0, areaBasin1[iSite1])
        #Avoid divide by 0
        frac = 0
        if (potentialArea>0):
          frac = areaCommon/potentialArea
        fracOverlapMax[iSite0, iSite1] = frac
  
  if (True):
    #print out some quick diagnostics
    print "Fractional area of overlapping advection\n", fracOverlap
  
  if (doMaxArea):
    return fracOverlap, fracOverlapMax
  else:
    return fracOverlap

def calc_fracOverlap_PT(sites0, sites1, cell2Site0, cell2Site1, theta0, theta1, siteIsMin):
  """Calculate fractional "vertical" overlap in terms of potential temperature """
  
  #If the "air mass" persists, the PT range should overlap between corresponding TPVs.
  
  #calculate bounds for each tpv that potentially matches another
  nSites0 = len(sites0); nSites1 = len(sites1)
  fracOverlap = np.zeros((nSites0, nSites1), dtype=float)
  
  min0 = np.empty(nSites0, dtype=float); max0 = np.empty(nSites0, dtype=float); len0 = np.empty(nSites0, dtype=float);
  min1 = np.empty(nSites1, dtype=float); max1 = np.empty(nSites1, dtype=float); len1 = np.empty(nSites1, dtype=float);
  for ind in xrange(nSites0):
    site = sites0[ind]
    minVal, maxVal = basinMetrics.get_minMax_cell2Site(site, cell2Site0, theta0)
    min0[ind] = minVal; max0[ind] = maxVal; len0[ind] = maxVal-minVal
  for ind in xrange(nSites1):
    site = sites1[ind]
    minVal, maxVal = basinMetrics.get_minMax_cell2Site(site, cell2Site1, theta1)
    min1[ind] = minVal; max1[ind] = maxVal; len1[ind] = maxVal-minVal
  
  for iSite0 in xrange(nSites0):    
    for iSite1 in xrange(nSites1):
      #intersection of 2 ordered line segments
      rangeTop = min(max0[iSite0],max1[iSite1]) #upper bound of overlap is smaller of the maxima
      rangeBottom = max(min0[iSite0],min1[iSite1])
      lenRange = rangeTop-rangeBottom
      lenPossible = max(len0[iSite0], len1[iSite1])
      
      #treat boundary case issues of no overlap
      frac = 0.0
      if (lenRange>0 and lenPossible>0):
        frac = lenRange/lenPossible
      
      fracOverlap[iSite0, iSite1] = frac
  
  return fracOverlap

def get_correspondMetrics(dataMetrics, sitesOut, iTime):
  """Read metrics for correspondence (used in old version)"""
  #It's slow to load 1 value at a time from file.
  #So, we can get (load,calculate?) values for all sites at a given time.
  
  #"extremum" metrics may not be too robust:
  #-thetaExtr: value jumps if, say, TPV goes to surface
  #-latExtr: seems okay, right?
  
  #area of non-overlap of advected tpvs would be useful if the "outer" filaments/chunks didn't break off and join willy-nilly
  #(see http://docs.scipy.org/doc/numpy/reference/routines.set.html)
  diffKeys = ['thetaExtr', 'latExtr', 'rEquiv', 'vortMean']; nKeys = len(diffKeys)
  refDiffs = [1.0, 1.0, 100.e3, 2.e-5]
  
  nSites = dataMetrics.variables['nSites'][iTime]
  allSites = dataMetrics.variables['sites'][iTime,:]; allSites = allSites[0:nSites]
  #isSiteReq = np.array([i in sitesOut for i in allSites], dtype=bool)
  isSiteReq = np.in1d(allSites, sitesOut, assume_unique=True)
  
  nSitesOut = len(sitesOut)
  valsOut = np.empty((nKeys,nSitesOut),dtype=float)
  for iKey in xrange(nKeys):
    key = diffKeys[iKey]
    vals = dataMetrics.variables[key][iTime,:]
    vals = vals[0:nSites]
    
    valsOut[iKey,:] = vals[isSiteReq]
    
  return (valsOut, refDiffs)
  
def calc_basinSimilarity(vals0, vals1, refDiffs):
  """Calculate 2 basins' similarity based on metrics' persistence (used in old version)"""
  #Input vals[variable,sites]
  
  nKeys, nSites1 = vals1.shape
  d = np.zeros(nSites1,dtype=float)
  for iKey in xrange(nKeys):
    diffs = np.absolute( vals1[iKey,:]-vals0[iKey] )
    d += diffs/refDiffs[iKey]
  
  return d
  
def correspond(sites0, cell2Site0, u0, v0, dt, 
               sites1, cell2Site1, u1, v1, mesh,
               trackMinMaxBoth, fracOverlapThresh,
               iTime0, dataMetrics, theta0, theta1):
  """
  Identify major and minor correspondences between basins at time0 and time1=time0+deltaT using persistence of basins' metrics 
  (used in old version)
  
  Arguments:
  sites{0,1} - indices of basin extrema at t{0,1}
  cell2Site{0,1} - segmentation map of cells to basins at t{0,1}
  u{0,1} - zonal velocity at t{0,1}
  v{0,1} - meridional velocity at t{0,1}
  dt - timestep t1=t0+dt
  mesh - Mesh instance
  trackMinMaxBoth - 0 tracks minima, 1 tracks maxima, 2 tracks both together (which doesn't make physical sense since lows could correspond with highs)
  fracOverlapThresh - Threshold of minimum fractional horizontal overlap that must exist to consider quantitative similarity of two basins [0,1]
  iTime0 - time index of t0
  dataMetrics - tpvTrack metrics netCDF4 object
  theta{0,1} - tropopause potential temperature at t{0,1}
  """
  
  #additional filters ---------------------
  #(1) feature properties a la cost function:
  #reasonable changes over "short" time in 
  #-thetaMin, circulation, feature velocity, location,...
  #(2) event consistency:
  #propagate, split, merge, genesis, lysis should have characteristics
  #-propagate: similar intensity and area
  #-split: similar sumArea
  #-merge:
  
  #these really end up needing to include matches not found in fracOverlap
  #when consider cases like small TPVs breaking off of a large one
  
  #decide whether sites correspond --------------------------
  #area overlap
  fracOverlap = calc_fracOverlap_advection(sites0, cell2Site0, u0, v0, dt, sites1, cell2Site1, u1, v1, mesh)
  isMatch = fracOverlap>fracOverlapThresh
  print "Number of matches after horizontal overlap: {0}".format(np.sum(isMatch))
  
  #potential temperature overlap
  if (False):
    isMatch = check_overlap_PT(isMatch, sites0, sites1, cell2Site0, cell2Site1, theta0, theta1)
    print "Number of matches after vertical overlap: {0}".format(np.sum(isMatch))
  
  #decide type of site correspondence (major vs. minor) ------------------
  #0-noMatch, 1-minor, 2-major
  nSites0 = len(sites0); nSites1 = len(sites1);
  
  metrics0, refDiffs = get_correspondMetrics(dataMetrics, sites0, iTime0); #print metrics0
  metrics1, refDiffs = get_correspondMetrics(dataMetrics, sites1, iTime0+1); #print metrics1
  
  #major time0->time1
  typeMatch01 = isMatch.copy().astype(int)
  for iSite0 in xrange(nSites0):
    site0 = sites0[iSite0]
    corrSites = sites1[isMatch[iSite0,:]>0]
    
    if (len(corrSites)<1):
      continue
    if (len(corrSites)==1):
      site1 = corrSites[0]
      typeMatch01[iSite0,sites1==site1] = 2
    else:
      vals0 = metrics0[:,iSite0]
      vals1 = metrics1[:,isMatch[iSite0,:]>0]
      d = calc_basinSimilarity(vals0, vals1, refDiffs); #print d
      
      minInd = np.argmin(d)
      similarSite = corrSites[minInd]
      typeMatch01[iSite0,sites1==similarSite] = 2
  
  #major time1<-time0 
  typeMatch10 = isMatch.copy().astype(int)
  for iSite1 in xrange(nSites1):
    site1 = sites1[iSite1]
    corrSites = sites0[isMatch[:,iSite1]>0]
    
    if (len(corrSites)<1):
      continue
    if (len(corrSites)==1):
      site0 = corrSites[0]
      typeMatch10[sites0==site0,iSite1] = 2
    else:
      vals0 = metrics1[:,iSite1]
      vals1 = metrics0[:,isMatch[:,iSite1]>0]
      d = calc_basinSimilarity(vals0, vals1, refDiffs);
      
      minInd = np.argmin(d)
      similarSite = corrSites[minInd]
      typeMatch10[sites0==similarSite,iSite1] = 2
  
  typeMatch = np.minimum(typeMatch01, typeMatch10) #e.g., site0.a-site1 not major if site0.a splits from site0 into site1 but site0.b more similar to site1
  print "Number of {0}s in 0->1 and 1<-0: {1}, {2}".format(2, np.sum(typeMatch01==2), np.sum(typeMatch10==2))
  print "Number of -major- correspondences: {0}".format(np.sum(typeMatch==2))
  return typeMatch

def correspond_overlap(sites0, cell2Site0, u0, v0, dt, 
               sites1, cell2Site1, u1, v1, mesh,
               trackMinMaxBoth, fracOverlapThresh,
               theta0, theta1):
  """
  Identify major and minor correspondences between basins at time0 and time1=time0+deltaT using horizontal and vertical overlap.
  
  Arguments:
  sites{0,1} - indices of basin extrema at t{0,1}
  cell2Site{0,1} - segmentation map of cells to basins at t{0,1}
  u{0,1} - zonal velocity at t{0,1} (m/s)
  v{0,1} - meridional velocity at t{0,1} (m/s)
  dt - timestep t1=t0+dt (s)
  mesh - Mesh instance
  trackMinMaxBoth - 0 tracks minima, 1 tracks maxima, 2 tracks both together (which doesn't make physical sense since lows could correspond with highs)
  fracOverlapThresh - Threshold of minimum fractional horizontal overlap that must exist to consider quantitative similarity of two basins [0,1]
  theta{0,1} - tropopause potential temperature at t{0,1}
  """
  #horizontal overlap under advection, and
  #overlap fraction defines similarity
  
  #maybe the PT range is more "accurate" than our estimated advection+horizOverlap, so we can weight it more...?
  #doing d = 1.0*fracOverlapHoriz+1.5*fracOverlapPT doesn't continue 2006/9/9 track. 
  #Maybe use similarity that penalizes not having both...overlapHoriz*PT?
  #Note that calc_fracOverlap_advection() uses potentialArea = min(areaBasin0, areaBasin1[iSite1]). For similarity, we want big to correspond to big
  #so edited advectionOverlap code to return this option as well.
  
  #decide whether sites correspond --------------------------
  #area overlap
  fracOverlap, fracOverlapMax = calc_fracOverlap_advection(sites0, cell2Site0, u0, v0, dt, sites1, cell2Site1, u1, v1, mesh, doMaxArea=True)
  isMatch = fracOverlap>fracOverlapThresh
  print "Number of matches after horizontal overlap: {0}".format(np.sum(isMatch))
  fracOverlapPT = calc_fracOverlap_PT(sites0, sites1, cell2Site0, cell2Site1, theta0, theta1, trackMinMaxBoth==0)
  #tanh(2)~.96 so can construct profile s.t. areaOverlapThresh->1
  #tanhGoal = .66; tanhFac = 2./tanhGoal
  
  #decide type of site correspondence (major vs. minor) ------------------
  #0-noMatch, 1-minor, 2-major
  nSites0 = len(sites0); nSites1 = len(sites1);
  
  #major time0->time1
  typeMatch01 = isMatch.copy().astype(int)
  for iSite0 in xrange(nSites0):
    site0 = sites0[iSite0]
    corrSites = sites1[isMatch[iSite0,:]>0]
    
    if (len(corrSites)<1):
      continue
    if (len(corrSites)==1):
      site1 = corrSites[0]
      typeMatch01[iSite0,sites1==site1] = 2
    else:
      #d = wtHoriz*fracOverlap[iSite0, isMatch[iSite0,:]>0]+wtVert*fracOverlapPT[iSite0, isMatch[iSite0,:]>0]
      #d = np.tanh(tanhFac*fracOverlapMax[iSite0, isMatch[iSite0,:]>0])*fracOverlapPT[iSite0, isMatch[iSite0,:]>0]
      d = fracOverlapMax[iSite0, isMatch[iSite0,:]>0] + fracOverlapPT[iSite0, isMatch[iSite0,:]>0]
      #d = fracOverlapMax[iSite0, isMatch[iSite0,:]>0]*fracOverlapPT[iSite0, isMatch[iSite0,:]>0]
      minInd = np.argmax(d); #print d,'\n',d[minInd]
      
      similarSite = corrSites[minInd]
      typeMatch01[iSite0,sites1==similarSite] = 2
  
  #major time1<-time0 
  typeMatch10 = isMatch.copy().astype(int)
  for iSite1 in xrange(nSites1):
    site1 = sites1[iSite1]
    corrSites = sites0[isMatch[:,iSite1]>0]
    
    if (len(corrSites)<1):
      continue
    if (len(corrSites)==1):
      site0 = corrSites[0]
      typeMatch10[sites0==site0,iSite1] = 2
    else:
      #d = wtHoriz*fracOverlap[isMatch[:,iSite1]>0, iSite1]+wtVert*fracOverlapPT[isMatch[:,iSite1]>0, iSite1]
      #d = np.tanh(tanhFac*fracOverlapMax[isMatch[:,iSite1]>0, iSite1])*fracOverlapPT[isMatch[:,iSite1]>0, iSite1]
      d = fracOverlapMax[isMatch[:,iSite1]>0, iSite1] + fracOverlapPT[isMatch[:,iSite1]>0, iSite1]
      #d = fracOverlapMax[isMatch[:,iSite1]>0, iSite1]*fracOverlapPT[isMatch[:,iSite1]>0, iSite1]
      minInd = np.argmax(d); #print d,'\n',d[minInd]
      
      similarSite = corrSites[minInd]
      typeMatch10[sites0==similarSite,iSite1] = 2
  
  typeMatch = np.minimum(typeMatch01, typeMatch10) #e.g., site0.a-site1 not major if site0.a splits from site0 into site1 but site0.b more similar to site1
  print "Number of {0}s in 0->1 and 1<-0: {1}, {2}".format(2, np.sum(typeMatch01==2), np.sum(typeMatch10==2))
  print "Number of -major- correspondences: {0}".format(np.sum(typeMatch==2))
  return typeMatch

def run_correspond(fNameOut, dataMetr, dataSeg, mesh, dt, 
                   trackMinMaxBoth, fracOverlapThresh, iTimeStart, iTimeEnd, dataMetrics):
  """
  Driver for correspondence module
  
  Arguments:
  fNameOut - Output file path
  dataMetr - tpvTrack preprocessing netCDF4 object
  dataSeg - tpvTrack segmentation netCDF4 object
  mesh - Mesh instance
  dt - timestep (s)
  trackMinMaxBoth - 0 tracks minima, 1 tracks maxima, 2 tracks both together (which doesn't make physical sense since lows could correspond with highs)
  fracOverlapThresh - Threshold of minimum fractional horizontal overlap that must exist to consider quantitative similarity of two basins [0,1]
  iTimeStart - Global time index to start correspondences
  iTimeEnd - Global time index to end correspondences (last correspondence is between iTimeEnd-1 <---> iTimeEnd)
  dataMetrics - tpvTrack metrics netCDF4 object (used in old version with similarity based on basins' metrics persistence)
  """
  #file for correspondences
  maxNSites = max(np.max(dataSeg.variables['nSitesMin'][:]), np.max(dataSeg.variables['nSitesMax'][:])); print "Maximum of {0} sites at any time".format(maxNSites)
  nTimes = iTimeEnd+1
  dataCorr = write_corr_netcdf_header(fNameOut, 'test', maxNSites, nTimes)
  
  for iTime in xrange(iTimeStart,iTimeEnd): #iTimeEnd will be the end of the correspondences
    #segmentation data
    cell2Site0 = dataSeg.variables['cell2Site'][iTime,:]
    sitesMin0 = dataSeg.variables['sitesMin'][iTime,:];
    nMin0 = dataSeg.variables['nSitesMin'][iTime]; sitesMin0 = sitesMin0[0:nMin0]
    sitesMax0 = dataSeg.variables['sitesMax'][iTime,:]; 
    nMax0 = dataSeg.variables['nSitesMax'][iTime]; sitesMax0 = sitesMax0[0:nMax0]
    
    cell2Site1 = dataSeg.variables['cell2Site'][iTime+1,:]
    sitesMin1 = dataSeg.variables['sitesMin'][iTime+1,:];
    nMin1 = dataSeg.variables['nSitesMin'][iTime+1]; sitesMin1 = sitesMin1[0:nMin1]
    sitesMax1 = dataSeg.variables['sitesMax'][iTime+1,:]; 
    nMax1 = dataSeg.variables['nSitesMax'][iTime+1]; sitesMax1 = sitesMax1[0:nMax1]
    
    #metr data
    u0 = dataMetr.variables['u'][iTime,:]; v0 = dataMetr.variables['v'][iTime,:]
    theta0 = dataMetr.variables['theta'][iTime,:]
    u1 = dataMetr.variables['u'][iTime+1,:]; v1 = dataMetr.variables['v'][iTime+1,:]
    theta1 = dataMetr.variables['theta'][iTime+1,:]
    
    #which basins we want to track ----------------------------
    sites0 = []; sites1 = []
    if (trackMinMaxBoth == 0): #just minima
      sites0 = sitesMin0
      sites1 = sitesMin1
    elif (trackMinMaxBoth == 1): #just maxima
      sites0 = sitesMax0
      sites1 = sitesMax1
    else: #track min+max
      print "Do you really want minima to be able to correspond to maxima?"
      sites0 = np.concatenate((sitesMin0,sitesMax0))
      sites1 = np.concatenate((sitesMin1,sitesMax1))
      
    #time correspondence ------------------
    if (False):
      typeMatch = correspond(sites0, cell2Site0, u0, v0, dt,
                           sites1, cell2Site1, u1, v1, mesh,
                           trackMinMaxBoth, fracOverlapThresh,
                           iTime, dataMetrics, theta0, theta1)
    else:
      typeMatch = correspond_overlap(sites0, cell2Site0, u0, v0, dt,
                           sites1, cell2Site1, u1, v1, mesh,
                           trackMinMaxBoth, fracOverlapThresh, theta0, theta1)
                           
    write_corr_iTime_netcdf(dataCorr, iTime, sites0, sites1, typeMatch)
  
  dataCorr.close()

def write_corr_netcdf_header(fName, info, maxNSites, nTimes):
  '''
  Make and write header for correspondence output file
  
  Arguments:
  fName - output file path
  info - additional description added to output metadata
  maxNSites - Maximum number of sites/basins at any time (since I don't know how to make ragged arrays)
  nTimes - number of timesteps
  
  Previous version used pickle.
  Unpickling objects gets expensive for long track times since we have to go sequentially load the relevant times.
  
  netcdf allows more direct/quicker access
  '''
  data = netCDF4.Dataset(fName, 'w', format='NETCDF4')
  data.description = info
  
  # dimensions
  data.createDimension('maxNSites', maxNSites)
  data.createDimension('nTimes', nTimes)
  
  # variables
  nSites0_data = data.createVariable('nSites0', 'i4', ('nTimes',))
  sites0_data = data.createVariable('sites0', 'i4', ('nTimes','maxNSites',))
  nCorrSites_data = data.createVariable('nCorrSites', 'i4', ('nTimes','maxNSites',))
  corrSites_data = data.createVariable('corrSites', 'i4', ('nTimes','maxNSites','maxNSites',))
  corrTypes_data = data.createVariable('corrTypes', 'i4', ('nTimes','maxNSites','maxNSites',))
  
  #Descriptions
  nSites0_data.description = 'Number of sites'
  sites0_data.description = 'Indices of sites'
  nCorrSites_data.description='Number of corresponding sites'
  corrSites_data.description = 'Indices of corresponding sites'
  corrTypes_data.description = 'Correspondence category [0-none,1-minor,2-major]'
  
  return data

def write_corr_iTime_netcdf(data, iTime, sites0, sites1, typeMatch):
  '''
  Write one time into correspondence output file
  
  Arguments:
  data - output netCDF4 object
  iTime - time index
  sites{0,1] - basin extrema at t{0,1}
  typeMatch - the type of correspondence (0-none, 1-minor, 2-major)
  
  variables to store are:
  {iTime, sites0, correspondingSites1[iSite0][iCorrespondingSites1], correspondenceType[iSite0][iCorrespondingSites1]}
  '''
  
  nSites0 = len(sites0); nSites1 = len(sites1)
  
  data.variables['nSites0'][iTime] = nSites0
  
  for iSite0 in xrange(nSites0):    
    corr1 = typeMatch[iSite0,:]
    corrSites = sites1[corr1>0]; nCorr = len(corrSites)
    typeCorr = corr1[corr1>0]
    
    data.variables['sites0'][iTime,iSite0] = sites0[iSite0]
    data.variables['nCorrSites'][iTime,iSite0] = nCorr
    data.variables['corrSites'][iTime,iSite0,0:nCorr] = corrSites
    data.variables['corrTypes'][iTime,iSite0,0:nCorr] = typeCorr
    
def plot_correspondences(fDirSave, fCorr, nTimes, mesh, iTimeStart=0):
  """
  Example plot of correspondences on map.
  Basins at t0 are plotted as green +. Basins at t1 are plotted as red o.
  Minor correspondence are denoted by red, thinner lines. Major are denoted by blue, thicker lines.
  """
  
  dataCorr = netCDF4.Dataset(fCorr,'r')
  
  m = Basemap(projection='ortho',lon_0=0,lat_0=89.5, resolution='l')
  r2d = 180./np.pi
  
  for iTime in xrange(iTimeStart,nTimes):
    plt.figure()
    m.drawcoastlines()
    
    allSites0, corrSites, typesCorr = read_corr_iTime(dataCorr, iTime)
    nSites0 = len(allSites0)
    for iSite in xrange(nSites0):
      site0 = allSites0[iSite]
      sites1 = corrSites[iSite]; nSites1 = len(sites1)
      minorMajor = typesCorr[iSite]
      if (nSites1<1):
        continue
      
      lat0, lon0 = mesh.get_latLon_inds(site0)
      lat1, lon1 = mesh.get_latLon_inds(np.array(sites1,dtype=int))
      lat0 = lat0*r2d; lon0 = lon0*r2d;
      lat1 = lat1*r2d; lon1 = lon1*r2d
      
      x0,y0 = m(lon0,lat0)
      m.scatter(x0,y0, marker='+', color='g', s=55)
      for iSite1 in xrange(nSites1):
        c = 'r'; lw=.5
        if (minorMajor[iSite1]>1):
          c='b'; lw=2
        m.drawgreatcircle(lon0, lat0, lon1[iSite1], lat1[iSite1], del_s=50.0, color=c, lw=lw)
        x1,y1 = m(lon1[iSite1], lat1[iSite1])
        #maybe change fillstyle of one's non-major markers
        markerStyle = 'x'
        if (minorMajor[iSite1]>1):
          markerStyle = 'o'
        m.scatter(x1,y1, marker=markerStyle, color='r', s=20)
    
    if (False):
      plt.show()
    else:
      fName = 'corr_debug_{0}.png'.format(iTime)
      fSave = fDirSave+fName
      print "Saving file to: "+fSave
      plt.savefig(fSave); plt.close()
      
def read_corr_iTime(data, iTime):
  """
  Read/return correspondences for the specified time index
  """
  #later accessed as correspondingSites1[iSite0][iCorrespondingSites1], correspondenceType[iSite0][iCorrespondingSites1]
  
  nSites0 = data.variables['nSites0'][iTime]
  sites0 = data.variables['sites0'][iTime,0:nSites0]
  nCorr = data.variables['nCorrSites'][iTime,0:nSites0]
  paddedCorrSites = data.variables['corrSites'][iTime,:,:]
  paddedTypeCorr = data.variables['corrTypes'][iTime,:,:]
  
  corrSites = []
  corrTypes = []
  for iSite0 in xrange(nSites0):
    corrSites.append(paddedCorrSites[iSite0, 0:nCorr[iSite0]])
    corrTypes.append(paddedTypeCorr[iSite0, 0:nCorr[iSite0]])
  
  return (sites0, corrSites, corrTypes)
  
def get_correspondingSites(dataCorr, iTime, site):
  """
  Return the sites that correspond to specified site at time index iTime
  """
  allSites0, corrSites, typeCorr = read_corr_iTime(dataCorr, iTime)
  if (site not in allSites0):
    print "Uhoh, site doesn't correspond to another..."
    print site, allSites0
  iSite = np.where(allSites0==site)[0][0]
  return (corrSites[iSite], typeCorr[iSite])
