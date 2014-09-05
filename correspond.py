import numpy as np
import netCDF4
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
  
def advect_LatLon(u, v, lat, lon, dt, r ):
  #return new lat/lon coordinates based on:
  #u,v in m/s, lat/lon in radians, dt in s, rSphere in m
  
  #u = r cos(lat) dLon/dt, v = r dLat/dt
  #constrain latitudes to be in [-pi/2, pi/2] and longitudes in [0, 2pi)
  
  dLon_dt = u/(r*np.cos(lat)) #if lat=+- pi/2, this will be big
  dLat_dt = v/r
  
  lat = lat+ dLat_dt*dt; lon = lon+ dLon_dt*dt
  #want to bound lat/lon. Note that for silly big meridional velocities, the following
  #can adjust to outside of poles, but then the whole tracking idea is screwed anyway.
  
  #bound lat. 
  #imagine crossing north pole so 89->93 is really 89->87N and 180degree switch in longitude
  crossedN = latOut>np.pi/2
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
  
  latInds, lonInds = segment_ll_flat.index_1dTo2d(cellInds, len(lon))
  latFeature = lat[latInds]; lonFeature = lon[lonInds]
  uFeature = u[latInds, lonInds]; vFeature = v[latInds, lonInds]
  
  newLat, newLon = advect_LatLon(uFeature, vFeature, latFeature, lonFeature, dt, mesh.r )
  
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
  mesh.get_closestCell2Pt(latPt, lonPt)
  
  for iPt in xrange(nPts):
    advCells[iPt] = mesh.get_closestCell2Pt(latPts[iPt], lonPts[iPt])
  
  #multiple cells can advect into the same cell
  advCells = np.unique(advCells)
  
  return advCells

def getCommon_1dInd(inds0, inds1):
  #return common values between lists with unique values
  return np.intersect1d(inds0, inds1, assume_unique=True)

def make_candidateCorrespondence_fracOverlap(siteInds0, cell2Site0, u0, v0,
                            siteInds1, cell2Site1, u1, v1, lat, lon, dt, r, areaLatCell):
  #Given fields+basins at t0 and t0+dt,
  #-create candidate matches by overlapping advection
  #in principle, overlap could mean: 
  #-# of common cells, min threshold for common area, min threshold for convex hulls overlapping area,...
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


