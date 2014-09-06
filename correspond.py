import numpy as np
import netCDF4
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
  
def advect_LatLon(u, v, latIn, lonIn, dt, r ):
  #return new lat/lon coordinates based on:
  #u,v in m/s, lat/lon in radians, dt in s, rSphere in m
  
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
  #return advected lat, lon points of feature cells
  
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
  #return cells inds of cells that basin advects to next time (+ or - dt)
  
  #coordinates of advected cell centers. 
  #so, we're assuming that mesh is dense enough wrt feature size that cell centers are sufficient?
  #or, overlap is sufficient that this will id it
  latPts, lonPts = advect_feature(siteInd, cell2Site, mesh, u, v, dt)
  
  #cells corresponding to those coordinates
  nPts = len(latPts);
  advCells = np.empty(nPts, dtype=int)
  
  for iPt in xrange(nPts):
    advCells[iPt] = mesh.get_closestCell2Pt(latPts[iPt], lonPts[iPt])
  
  if (False):
    latCell, lonCell = mesh.get_latLon_inds(advCells)
    print latPts; print lonPts
    print latCell; print lonCell
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
  #return common values between lists with unique values
  return np.intersect1d(inds0, inds1, assume_unique=True)
  
def calc_fracOverlap_advection(sites0, cell2Site0, u0, v0, dt,
                               sites1, cell2Site1, u1, v1, mesh):
  #Given fields+basins at t0 and t0+dt,
  #-create candidate matches by overlapping advection
  #in principle, overlap could mean: 
  #-# of common cells, min threshold for common area, min threshold for convex hulls overlapping area,...
  #return matrix with fraction of overlap by area(nCellsMatch)/area(nCellsPossible)
  
  nSites0 = len(sites0); nSites1 = len(sites1)
  fracOverlap = np.zeros((nSites0, nSites1), dtype=float)
  
  #store advection of t1 sites -dt/2
  sites2Cells_t1 = [None]*nSites1
  areaBasin1 = np.empty(nSites1,dtype=float)
  for iSite1 in xrange(nSites1):
    siteInd = sites1[iSite1]
    #advect basin -dt/2
    sites2Cells_t1[iSite1] = advect_basin(siteInd, cell2Site1, mesh, u1, v1, -.5*dt)
    areaBasin1[iSite1] = np.sum( mesh.get_area_inds(sites2Cells_t1[iSite1]) )
    
  #see which t0 sites advected dt/2 overlap with future sites advected back
  for iSite0 in xrange(nSites0):
    siteInd = sites0[iSite0]
    #advect basin +dt/2
    site2Cells_t0 = advect_basin(siteInd, cell2Site0, mesh, u0, v0, .5*dt)
    areaBasin0 = np.sum( mesh.get_area_inds(site2Cells_t0) )
    
    for iSite1 in xrange(nSites1):
      #for frac overlap, there's a choice for what cells to use.
      #consider candidate big t0 and small t1.
      #-for min(cellsInBasin), fraction will be high w/ few cells from big needed
      #-for max(cellsInBasin), fraction will be low even if lots of small covered
      commonCells = getCommon_1dInd(site2Cells_t0, sites2Cells_t1[iSite1])
      areaCommon = np.sum( mesh.get_area_inds(commonCells) )
      
      potentialArea = min(areaBasin0, areaBasin1[iSite1])
      frac = areaCommon/potentialArea
      fracOverlap[iSite0, iSite1] = frac
  
  if (True):
    #print out some quick diagnostics
    print "For overlapping advection", fracOverlap
  
  return fracOverlap

def correspond(sites0, cell2Site0, u0, v0, dt, 
               sites1, cell2Site1, u1, v1, mesh,
               trackMinMaxBoth, fracOverlapThresh):
  
  #area overlap -------------------
  fracOverlap = calc_fracOverlap_advection(sites0, cell2Site0, u0, v0, dt, sites1, cell2Site1, u1, v1, mesh)
  
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
  
  #decide correspondence --------------------------
  isMatch = fracOverlap>fracOverlapThresh
  print "Number of matches from correspondence: {0}".format(np.sum(isMatch))
  
  return isMatch

def run_correspond(fNameOut, dataMetr, dataSeg, mesh, dt, 
                   trackMinMaxBoth, fracOverlapThresh, iTimeStart, iTimeEnd):
  
  #write tracks to a text file?
  fCorr = open(fNameOut,'w')
  
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
    u1 = dataMetr.variables['u'][iTime+1,:]; v1 = dataMetr.variables['v'][iTime+1,:]
    
    #which basins we want to track ----------------------------
    sites0 = []; sites1 = []
    if (trackMinMaxBoth == 0): #just minima
      sites0 = sitesMin0
      sites1 = sitesMin1
    elif (trackMinMaxBoth == 1): #just maxima
      sites0 = sitesMax0
      sites1 = sitesMax1
    else: #track min+max
      sites0 = np.concatenate((sitesMin0,sitesMax0))
      sites1 = np.concatenate((sitesMin1,sitesMax1))
      
    #time correspondence ------------------
    isMatch = correspond(sites0, cell2Site0, u0, v0, dt,
                         sites1, cell2Site1, u1, v1, mesh,
                         trackMinMaxBoth, fracOverlapThresh)
                         
    write_seg_iTime(fCorr, iTime, sites0, sites1, isMatch)
  
  fCorr.close()
  
def write_seg_iTime(f, iTime, sites0, sites1, isMatch):
  '''
  format is:
  Time iTime0 nSites0
  site0[0] nCorrespondingSites1 : correspondingSites1
  site0[1] nCorrespondingSites1 : correspondingSites1
  ...
  Time iTime1 nSites0
  ...
  '''
  
  nSites0 = len(sites0); nSites1 = len(sites1)
  
  s = 'Time {0} {1}\n'.format(iTime, nSites0)
  f.write(s)
  for iSite0 in xrange(nSites0):
    corrSites = sites1[isMatch[iSite0,:]>0]
    nCorr = len(corrSites)
    s = '{0} {1} : '.format(sites0[iSite0], nCorr)
    sCorr = '\n'
    if (nCorr>0):
      sCorr = np.array_str(corrSites, max_line_width=100000000000000)[1:-1]+'\n' #to skip the brackets
    f.write(s+sCorr)
    
def plot_correspondences(fDirSave, fCorr, nTimes, mesh):

  f = open(fCorr,'r')
  
  m = Basemap(projection='ortho',lon_0=0,lat_0=89.5, resolution='l')
  r2d = 180./np.pi
  
  while(nTimes>0):
    plt.figure()
    m.drawcoastlines()
    
    s = f.readline(); line = s.strip().split()
    iTime = int(line[1]); nSites0 = int(line[2])
    for iSite in xrange(nSites0):
      s = f.readline(); line = s.strip().split()
      site0 = int(line[0]); nSites1 = int(line[1])
      if (nSites1<1):
        continue
      sites1 = [int(i) for i in line[3:]]
      
      lat0, lon0 = mesh.get_latLon_inds(site0)
      lat1, lon1 = mesh.get_latLon_inds(np.array(sites1,dtype=int))
      lat0 = lat0*r2d; lon0 = lon0*r2d;
      lat1 = lat1*r2d; lon1 = lon1*r2d
      
      x0,y0 = m(lon0,lat0)
      m.scatter(x0,y0, marker='+', color='g', s=40)
      for iSite1 in xrange(nSites1):
        m.drawgreatcircle(lon0, lat0, lon1[iSite1], lat1[iSite1], del_s=100.0, color='b')
        x1,y1 = m(lon1[iSite1], lat1[iSite1])
        m.scatter(x1,y1, marker='o', color='r', s=15)
    
    if (False):
      plt.show()
    else:
      fName = 'corr_{0}.png'.format(iTime)
      fSave = fDirSave+fName
      print "Saving file to: "+fSave
      plt.savefig(fSave); plt.close()
      
    nTimes = nTimes-1
  
  f.close()



