import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic
import glob
import time
import datetime as dt
import os

#watershed has a few options for implementation:
#-for every cell, walk down steepest gradient to the basin
#-every local min is a voronoi map with distance as the sum of ?
#along path to min/voronoi site. benefit of this is ability to remove voronoi sites for iterative merging.
#could remove based on min size of watershed,...

#active contours to identify regions:
#-minimize: internal energy + external energy + curveCost(eg smoothness)
#can create logically square say nearest nbr resample of area N of 75N and ship over to image processing+

def calc_distSphere_multiple(r, lat1, lon1, lat2, lon2):
  '''
  #return the distance between 1 ll1 point and >=1 ll2 points.
  on a sphere.
  input lat/lon in radians!!
  '''
  
  dlat = lat2-lat1
  dlon = lon2-lon1
  latTerm = np.sin(.5*dlat); latTerm = latTerm*latTerm;
  lonTerm = np.sin(.5*dlon); lonTerm = lonTerm*lonTerm*np.cos(lat1)*np.cos(lat2);
  dAngle = np.sqrt(latTerm+lonTerm)
  dist = 2.*r*np.arcsin(dAngle)
  
  return(dist)

def gatherCells_region(iLat0, iLon0, nLat, nLon, latCell, lonCell, 
                       r, distRegion):
  #return list of cell indices within specified spatial region of given point  
  lat0 = latCell[iLat0]; lon0 = lonCell[iLon0]
  
  inRegion_lat = []; inRegion_lon = []
  nbrInds_lat, nbrInds_lon = nbrInds_ll(iLat0, iLon0, nLat, nLon)
  candidateRgn_lat = nbrInds_lat; candidateRgn_lon = nbrInds_lon;
  
  beenAdded = np.zeros((nLat,nLon),dtype=int)
  
  while (len(candidateRgn_lat)>0):
    iLat = candidateRgn_lat.pop(); iLon = candidateRgn_lon.pop()
    beenAdded[iLat,iLon] = 1
    d = calc_distSphere_multiple(r, lat0, lon0, latCell[iLat], lonCell[iLon])
    if (d<distRegion):
      inRegion_lat.append(iLat); inRegion_lon.append(iLon)
      
      nbrInds_lat, nbrInds_lon = nbrInds_ll(iLat, iLon, nLat, nLon)
      nNbrs = len(nbrInds_lat)
      for iNbr in xrange(nNbrs):
        iLat = nbrInds_lat[iNbr]; iLon = nbrInds_lon[iNbr]
        #add nbr if not previously added
        if (beenAdded[iLat,iLon]<1):
          candidateRgn_lat.append(iLat); candidateRgn_lon.append(iLon)
          beenAdded[iLat,iLon]=1
  
  return (inRegion_lat, inRegion_lon)

def calc_lonIndicesWithinLength(lats, nLon, r, distRegion):
  #return the number of longitude indices within given length at each specified latitude.
  #imagine this is like a pyramid with more at the poles and less at the equator.
  #since radius, take nLons[iLat] points both east and west for disk
  dRadLon = 2.*np.pi/nLon #[0,2pi)
  dLons = r*np.cos(lats)*dRadLon; #arc lengths of pieces of latitude circles
  dLons = np.maximum(dLons, 1.e-12) #elementwise to avoid /0
  nLons = np.floor(distRegion/dLons).astype(int);
  #since made to index left and right, as max:
  #want len(range(left, right))=nLons (remember, not w/ right endpt)
  nLons = np.minimum(nLons, nLon/2)
  
  #indices are: np.arange(-nLons/2, nLons/2)%nLons...maybe np.sort(...)
  return nLons

def calc_latIndicesWithinLength(nLats, r, distRegion):
  #return the number of latitude indices within given length
  dRadLat = np.pi/(nLats-1) #[-pi/2, pi/2]
  distSN = r*dRadLat #arc length South-North
  nNorth = int(np.floor(distRegion/distSN))
  
  #actual indices are:
  #indN = np.min((iLat+nNorth, nLat-1)); indS = np.max((iLat-nNorth, 0))
  return nNorth

def gatherInds_region_latBox_1AtPole(iLat0, iLon0, nLat, nLon, latCell, lonCell,
                             nLatIndsLength, nLonIndsLength, r, distRegion):
  #return list of lat,lon indices within specified spatial region of given point.
  #for efficiency, call 1x per latitude and just shift lon indices.
  #nLonIndsLength = calc_lonIndicesWithinLength(lats, nLon, r, distRegion)...1x per mesh
  #Note that self=[iLat0,iLon0] pair will be in the returned region.
  #for added efficiency, only return index of 1 value at pole since all really same point.
  
  #we'll do this by searching within the input bounding box.
  lat0 = latCell[iLat0]; lon0 = lonCell[iLon0]
  indNPole = np.max((iLat0-nLatIndsLength, 0)); indSPole = np.min((iLat0+nLatIndsLength, nLat-1))
  indN = np.max((indNPole, 1)); indS = np.min((indSPole, nLat-2))
  candLats = np.arange(indN,indS+1) #include endpoint
  
  #go east at each latitude until too far
  inRegion_lat = []; inRegion_lon = []
  for iLat in candLats:
    #since symmetric about longitude, just have to plus distance and include
    # 0 and -plus indices in the returned result.
    candLons = np.arange(iLon0+1,iLon0+nLonIndsLength[iLat]+1)%nLon
    d = calc_distSphere_multiple(r, lat0, lon0, latCell[iLat], lonCell[candLons])
    closeLons = candLons[d<distRegion]
    lonInds = np.concatenate((closeLons,iLon0+iLon0-closeLons)).tolist();
    lonInds.append(iLon0)
    
    inRegion_lon.extend(lonInds)
    inRegion_lat.extend([iLat]*len(lonInds))
  
  #add in 1 value at pole if w/in region
  if (indNPole==0):
    inRegion_lon.append(iLon0); inRegion_lat.append(indNPole)
  if (indSPole==nLat-1):
    inRegion_lon.append(iLon0); inRegion_lat.append(indSPole)
    
  inRegion_lon = np.array(inRegion_lon)%nLon
  inRegion_lat = np.array(inRegion_lat)
  return (inRegion_lat, inRegion_lon)                  
                         
def find_minCells_flat(vals, nLat, nLon):
  #return array[nCells] with 1 if cell is min

  isMin = np.zeros((nLat*nLon),dtype=int)
  for iLat  in xrange(nLat):
    for iLon in xrange(nLon):
      nbrs = nbrInds_ll_flat(iLat, iLon, nLat, nLon)
      ind0 = index_2dTo1d(iLat, iLon, nLon)
      val0 = vals[ind0]
      valNbrs = vals[nbrs]
      if (np.all(val0<=valNbrs)): #is site if can't descend from it
        isMin[iCell] = 1

  nMin = np.sum(isMin)
  print "Number of local min: ", nMin
  return isMin

def find_minCells_region_flat(vals, nLat, nLon, inRegion):
  #return array[nCells] with 1 if cell is min and cell in region
  
  isMin = np.zeros((nLat*nLon),dtype=int)
  for iLat  in xrange(nLat):
    for iLon in xrange(nLon):
      ind0 = index_2dTo1d(iLat, iLon, nLon)
      if (inRegion[ind0]<1):
        continue
      
      nbrs = nbrInds_ll_flat(iLat, iLon, nLat, nLon)
      val0 = vals[ind0]
      valNbrs = vals[nbrs]
      if (np.all(val0<=valNbrs)): #is site if can't descend from it
        isMin[ind0] = 1

  nMin = np.sum(isMin)
  print "Number of local min in region: ", nMin
  return isMin

def watershed_region(vals, cellIsMin, nLat, nLon, r, dRegion, latCell, lonCell, inRegion):
  #to make adding/deleting basins simple, follow gradient until reach a site.
  #map every cell to follow local steepest gradient. basins go to self.
  #filter basins so have to be a min within specified region (disk of radius dRegion)
  #return map of cell to basin.
  
  '''
  #to adapt global watershed to region, make values outside of region huge so 
  don't steepest descend that way. since we pass in minCells, do this before call:
  bigVal = 1.e10
  vals = np.copy(valsIn)
  vals[inRegion<1] = bigVal
  '''
  
  nCells = nLat*nLon
  cell2Site = -np.ones(nCells,dtype=int) #so no cell2Site[iCell]=iCell
  for iCell in xrange(nCells):
    if (inRegion[iCell]<1):
      continue
    if (cellIsMin[iCell]>0):
      cell2Site[iCell]= iCell

  #get local steepest path
  for iCell in xrange(nCells):
    if (inRegion[iCell]<1):
      continue
    if (cellIsMin[iCell]>0):
      continue
    
    iLat0, iLon0 = index_1dTo2d(iCell, nLon)
    nbrInds_lat, nbrInds_lon = nbrInds_ll(iLat0, iLon0, nLat, nLon)
    nbrInds = index_2dTo1d(np.array(nbrInds_lat), np.array(nbrInds_lon), nLon)
    #print iLat0, iLon0, nbrInds_lat; print nbrInds_lon; print nbrInds
    #nNbrs = len(nbrInds_lat)
    
    val0 = vals[iCell]
    valNbrs = vals[nbrInds]

    #correspondence is towards minimum gradient.
    lat0 = latCell[iLat0]; lon0 = lonCell[iLon0]
    latNbrs = latCell[nbrInds_lat]; lonNbrs = lonCell[nbrInds_lon]
    dNbrs = calc_distSphere_multiple(r, lat0, lon0, latNbrs, lonNbrs)
    dMin = r/1.e16; dNbrs[dNbrs<dMin]=dMin #avoid divide by 0
    #print valNbrs, dNbrs, val0
    valNbrs = (valNbrs-val0)/dNbrs
    iNbr = np.argmin(valNbrs)
    cell2Site[iCell] = nbrInds[iNbr]
  
  nRedirect = 0
  
  #Filter local extrema by area to limit high (spatial) frequency "noise".
  #For multiple close mins, the smallest counts as min for that region.
  #dRegion = 300.e3 #radius in meters of disk of filtering region
  nLatIndsLength = calc_latIndicesWithinLength(nLat, r, dRegion)
  nLonIndsLength = calc_lonIndicesWithinLength(latCell, nLon, r, dRegion)
    
  for iLat in xrange(nLat):
    iLonRef = 0
    inDiskLat, inDiskLon_ref = gatherInds_region_latBox_1AtPole(iLat, iLonRef, nLat, nLon, latCell, lonCell,
                                                nLatIndsLength, nLonIndsLength, r, dRegion)
    #
    for iLon in xrange(nLon):
      iCell = index_2dTo1d(iLat, iLon, nLon)
      if (inRegion[iCell]<1):
        continue
      if (cellIsMin[iCell]>0):
        #see if cell is min in region, not just neighbors.
        #if not regional min, update cell2Site so local min goes to another basin
        diffLonInd = iLon-iLonRef
        inDiskLon = (inDiskLon_ref+diffLonInd)%nLon

        cellsRegion = index_2dTo1d(inDiskLat, inDiskLon, nLon)
        valsRegion = vals[cellsRegion]
        minInd = np.argmin(valsRegion)
        minVal = valsRegion[minInd]; minCell = cellsRegion[minInd];
        val0 = vals[iCell]; #print val0, minVal
        if (minVal < val0):
          #print "Redirecting cell {0} to {1}".format(iCell, minCell)
          cellIsMin[iCell] = 0
          cell2Site[iCell] = minCell
          nRedirect = nRedirect+1
  print "Number of redirects for regional min: ", nRedirect
  
  #follow local steepest path (and any redirections from, say, regional thresholds) to site
  for iCell in xrange(nCells):
    if (inRegion[iCell]<1):
      #what should these correspond to?
      continue
    nextCell = cell2Site[iCell]
    nCount = 0
    while (not cellIsMin[nextCell]>0):
      nextCell = cell2Site[nextCell]
      #print "Cell {0} going to {1}".format(iCell,nextCell); print vals[iCell], vals[nextCell]
      nCount=nCount+1
      if (nCount>nCells):
        iLat, iLon = index_1dTo2d(iCell, nLon)
        print "Uhoh, stuck in while loop for iLat,iLon ({0},{1}) with value {2}".format(iLat, iLon, vals[iCell])
        break

    cell2Site[iCell] = nextCell

  return (cell2Site, cellIsMin)

def segment_high_low_watershed_region(lat, lon, theta, vort, r, inRegion, dRegion):
  #get high and low basin seeds, associate cells to both high and low basins if not extrema.
  #to decide whether "really" part of high or low basin, we have options:
  #-(anti-)cyclonic for (high) low...is local vorticity noisy?
  #-closer theta value to maxima a la color scale grouping...huge min or max value now matters
  #-whether steeper gradient is to high or low
  #-physical distance
  #-concavity of surface a la last closed contour
  
  #theta, vort, inRegion come in as 1d arrays
  #use inRegion for (1) region of interest (2)ignore missing values
  
  nLat = len(lat); nLon = len(lon); nCells = nLat*nLon
  
  #segment ------------------------------------------
  #dRegion = 300.e3 #radius in same units as r of disk for min basin
  
  #mins
  print "Finding minima"
  #to adapt global watershed to region, make values outside of region huge so don't steepest descend that way
  bigVal = 1.e10
  vals = np.copy(theta) #so don't affect variable passed in
  vals[inRegion<1] = bigVal
  
  cellIsMin = find_minCells_region_flat(vals, nLat, nLon, inRegion)
  cell2SiteMin, cellIsMin = watershed_region(vals, cellIsMin, nLat, nLon, r, dRegion, lat, lon, inRegion)
  
  #maxs: perform min on an inverted surface
  print "Finding maxima"
  #adapt global watershed to region
  vals = -np.copy(theta)
  vals[inRegion<1] = bigVal
  
  cellIsMax = find_minCells_region_flat(vals, nLat, nLon, inRegion)
  cell2SiteMax, cellIsMax = watershed_region(vals, cellIsMax, nLat, nLon, r, dRegion, lat, lon, inRegion)
  
  #"voting" procedure for low/high classification ------
  print "Associating to max or min"
  cell2Site = -np.ones(nCells, dtype=int)
  
  classifyMethod = 0 #see if-checks below for corresponding number for method
  if (classifyMethod==0):
    #local vorticity
    print("Classifying by local vorticity")
    for iCell in xrange(nCells):
      if (inRegion[iCell]<1):
        continue
      #cyclonic depends on hemisphere
      if (cellIsMin[iCell]>0 or cellIsMax[iCell]>0): #allows for cyclonic max. is that right?
        cell2Site[iCell] = iCell
      else:
        signHem = 1 #sign function is problematic since sign(0)=0
        iLat0, iLon0 = index_1dTo2d(iCell, nLon)
        if (lat[iLat0]<0): #lat=0 gets put in NH
          signHem = -signHem
        if (signHem*vort[iCell]<0): #anticyclonic
          cell2Site[iCell] = cell2SiteMax[iCell]
        else: #0 or cyclonic
          cell2Site[iCell] = cell2SiteMin[iCell]
  else:
    print "Error. Classification method not recognized/specified."
          
  return (cell2Site, cellIsMin, cellIsMax)

def get_missingCells_file(data):
  #if search vertical column down from top, there can be no 2pvu value if:
  #-entire column is above 2pvu: DT below sfc
  #-entire column is below 2pvu: wrong hemisphere or low pv column (tropics, anticyclonic,...)
  
  pvInd = 1
  isMissing = data.variables['TMP_P0_L109_GLL0'][pvInd,:,:].mask
  return isMissing
  
def fill_missingVals_region(valsIn, nLat, nLon, isMissing, inRegion):
  #fill value is average of non-missing neighbors
  
  vals = np.copy(valsIn)
  needFill = isMissing*inRegion;
  nNeedFill = np.sum(needFill); print "Filling {0} values".format(nNeedFill)
  
  while (np.sum(needFill)>0):
    for iLat in xrange(nLat):
      for iLon in xrange(nLon):
        if (needFill[iLat,iLon]>0): #True>0 is True, False==0 is True
          nbrInds_lat, nbrInds_lon = nbrInds_ll(iLat, iLon, nLat, nLon)
          nbrsNeedFill = needFill[nbrInds_lat, nbrInds_lon]
          if (False in nbrsNeedFill): #have neighbor with value
            #fill value is average of valid nbrs
            validNbrs = nbrsNeedFill==False;
            valsNbrs = vals[nbrInds_lat, nbrInds_lon]
            vals[iLat, iLon] = np.mean(valsNbrs[validNbrs])+1.e-10 #so don't have same value
            needFill[iLat,iLon]=False
  
  return vals
  
def get_segmentVars_file(data):
  #return data of variables needed from file, no time index.
  #return SI units
  
  #fname = '/data02/cases/2014/gfs_4_20140101_0000_123.nc'
  #data = netCDF4.Dataset(fname,'r')
  lat = data.variables['lat_0'][:] * np.pi/180. #deg 2 radians
  lon = data.variables['lon_0'][:] * np.pi/180.
  
  #for pvu levels, lv_PVL4 = [-2, 2]*E-6
  pvInd = 1
  tmp = data.variables['TMP_P0_L109_GLL0'][pvInd,:,:].data #K
  press = data.variables['PRES_P0_L109_GLL0'][pvInd,:,:].data #Pa
  
  u = data.variables['UGRD_P0_L109_GLL0'][pvInd,:,:].data #m/s
  v = data.variables['VGRD_P0_L109_GLL0'][pvInd,:,:].data
  
  return (lat, lon, u, v, tmp, press)

def calc_vertVorticity_ll(u, v, nLat, nLon, lat, r):
  '''
  Pulled from: http://www.ncl.ucar.edu/Document/Functions/Built-in/uv2vr_cfd.shtml :
  According to H.B. Bluestein [Synoptic-Dynamic Meteorology in Midlatitudes, 1992, 
  Oxford Univ. Press p113-114], 
  let D represent the partial derivative, a the radius of the earth, 
  phi the latitude and dx2/dy2 the appropriate longitudinal and latitudinal spacing, 
  respectively. Then, letting j be the latitude y-subscript, and i be the longitude x-subscript:

    rv = Dv/Dx - Du/Dy + (u/a)*tan(phi)


    rv(j,i) = (v(j,i+1)-v(j,i-1))/dx2(j)
              - (u(j+1,i)-u(j-1,i))/dy2(j)
              + (u(j,i)/a)*tan(phi(j)) #since meridians aren't parallel

  The last terms accounts for the convergence of the meridians on a sphere. 
  '''
  
  #we'll do the above centered finite differencing for the non-pole latitudes.
  #for the poles, we'll do a finite volume \int gradxu dA = \int u.n dS since 
  #trying to finite difference it confuses me. remember that the poles are really
  #nLon copies of the same point
  
  vort = np.empty((nLat, nLon), dtype=float)
  
  dRadLat = np.pi/(nLat-1) #[-pi/2, pi/2], ie with values at both poles
  dRadLon = 2.*np.pi/nLon #[0,2pi)
  dy = r*dRadLat; dy2 = 2.*dy #arc length on a sphere
  
  #calc values for non poles
  for iLat in xrange(1,nLat-1):
    tanphi = np.tan(lat[iLat])/r
    dx = r*np.cos(lat[iLat])*dRadLon; dx2 = 2.*dx
    for iLon in xrange(nLon):
      iWest = (iLon-1)%nLon # -1%4=3 so don't worry about negatives
      iEast = (iLon+1)%nLon
      iSouth = iLat+1; iNorth = iLat-1
      dv_dx = (v[iLat, iEast]-v[iLat, iWest])/dx2
      du_dy = (u[iNorth, iLon]-u[iSouth, iLon])/dy2
      meridTerm = u[iLat, iLon]*tanphi
      vort[iLat, iLon] = dv_dx-du_dy+meridTerm
      
  #calc values for north and south poles with finite volume approach
  #around next latitude equatorward of pole
  iLat = nLat-1;
  dx = r*np.cos(lat[iLat-1])*dRadLon #for evenly spaced lats, same dx for south and north
  #for area of the cap, http://mathworld.wolfram.com/SphericalCap.html
  a = dx/dRadLon
  h = r-np.sqrt(r*r-a*a)
  areaCap = 2.*np.pi*r*h
  #around south pole, remember integrate with domain on left
  undS = np.sum(u[iLat-1,:])*-dx #since +dx has domain on right
  vort[iLat,:] = undS/areaCap
  iLat = 0
  undS = np.sum(u[iLat+1,:])*dx
  vort[iLat,:] = undS/areaCap
  
  return vort

def calc_potentialTemperature(tmp, press):
    
  Cp = 1004.5; Rd = 287.04;
  Rd_cp = Rd/Cp; p0 = 1.e5
  theta = tmp*((p0/press)**Rd_cp)
  
  return theta

def plot_var_ll(lat, lonIn, varIn, useBounds):
  #Input lat and lon as 1d arrays in radians, var as latxlon
  #if len(useBounds)==2, plot with the given [min, max]
  var, lons = addcyclic(varIn, lonIn)
  lats,lons = np.meshgrid(lat, lons)
  var = np.transpose(varIn) #to be x,y for plotting instead of lat x lon
  plt.figure()

  #m = Basemap(projection='ortho',lon_0=100,lat_0=60, resolution='l')
  r2d = 180./np.pi
  m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  x,y = m(lons*r2d, lats*r2d)

  m.drawcoastlines()
  m.drawmapboundary()
  
  if (len(useBounds)==2):
    m.pcolor(x,y,var,shading='flat',edgecolors='none',cmap=plt.cm.RdBu_r, vmin=useBounds[0], vmax=useBounds[1])
  else:
    m.pcolor(x,y,var,shading='flat',edgecolors='none',cmap=plt.cm.RdBu_r)
  #cbar = m.colorbar(plt1,location='bottom',pad="5%")
  #cbar.set_label('\Delta m')
  plt.colorbar()
  plt.show()

def plot_segment(lat, lonIn, varIn, isMin, isMax):
  #Input lat and lon as 1d arrays in radians, var as latxlon
  #if len(useBounds)==2, plot with the given [min, max]
  var, lons = addcyclic(varIn, lonIn)
  lats,lons = np.meshgrid(lat, lons)
  #var = np.transpose(varIn) #to be x,y for plotting instead of lat x lon
  var = np.transpose(var)
  plt.figure()

  #m = Basemap(projection='ortho',lon_0=100,lat_0=60, resolution='l')
  r2d = 180./np.pi
  m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  x,y = m(lons*r2d, lats*r2d)
  #print x.shape, y.shape

  m.drawcoastlines()
  m.drawmapboundary()
  
  pPlot = m.pcolor(x,y,var,shading='flat',edgecolors='none',cmap=plt.cm.jet, vmin=280, vmax=360)
  
  isMin, lons = addcyclic(isMin, lonIn); isMin = np.transpose(isMin)
  isMax, lons = addcyclic(isMax, lonIn); isMax = np.transpose(isMax)
  xMin = x[isMin>0]; yMin = y[isMin>0]; m.scatter(xMin, yMin, c='m', marker="o")
  xMax = x[isMax>0]; yMax = y[isMax>0]; m.scatter(xMax, yMax, c='y', marker="o")
  
  plt.colorbar(pPlot)
  plt.show()

def plot_segment_save(fNameSave, lat, lonIn, varIn, isMin, isMax):
  #Input lat and lon as 1d arrays in radians, var as latxlon
  #if len(useBounds)==2, plot with the given [min, max]
  var, lons = addcyclic(varIn, lonIn)
  lats,lons = np.meshgrid(lat, lons)
  #var = np.transpose(varIn) #to be x,y for plotting instead of lat x lon
  var = np.transpose(var) #to be x,y for plotting instead of lat x lon
  plt.figure()

  #m = Basemap(projection='ortho',lon_0=100,lat_0=60, resolution='l')
  r2d = 180./np.pi
  m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  x,y = m(lons*r2d, lats*r2d)
  #print x.shape, y.shape

  m.drawcoastlines()
  m.drawmapboundary()
  
  pPlot = m.pcolor(x,y,var,shading='flat',edgecolors='none',cmap=plt.cm.jet, vmin=280, vmax=360)
  
  isMin, lons = addcyclic(isMin, lonIn); isMin = np.transpose(isMin)
  isMax, lons = addcyclic(isMax, lonIn); isMax = np.transpose(isMax)
  xMin = x[isMin>0]; yMin = y[isMin>0]; m.scatter(xMin, yMin, c='m', marker="o")
  xMax = x[isMax>0]; yMax = y[isMax>0]; m.scatter(xMax, yMax, c='y', marker="o")
  
  plt.colorbar(pPlot)
  plt.savefig(fNameSave, bbox_inches='tight'); plt.close()

def plot_segment_inds(lat, lonIn, cell2Site, isMin, isMax):
  
  nLat = len(lat); nLon = len(lonIn)
  
  var = -isMin[cell2Site]+isMax[cell2Site]
  var = unflatten_1dTo2d(var, nLat, nLon)
  
  var, lons = addcyclic(var, lonIn)
  lats,lons = np.meshgrid(lat, lons)
  var = np.transpose(var) #to be x,y for plotting instead of lat x lon
  plt.figure()

  #m = Basemap(projection='ortho',lon_0=100,lat_0=60, resolution='l')
  r2d = 180./np.pi
  m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  x,y = m(lons*r2d, lats*r2d)
  #print x.shape, y.shape

  m.drawcoastlines()
  m.drawmapboundary()
  
  pPlot = m.pcolor(x,y,var,shading='flat',edgecolors='none') #,cmap=plt.cm.jet, vmin=280, vmax=360)
  
  #for scatter plot, can also change marker styles
  #markerString = '${0}$'.format(ind)
  #m.scatter(x0, y0, s=100, marker=markerString)
  
  if (True):
    #sites = np.array([0], dtype=int);
    sites = np.unique(cell2Site[isMin])
    nSites = len(sites)
    iLat, iLon = index_1dTo2d(sites, nLon)
    for iSite in xrange(nSites):
      latInd = iLat[iSite]; lonInd = iLon[iSite]
      x0 = x[lonInd, latInd]; y0 = y[lonInd, latInd]; ind = sites[iSite]
      plt.text(x0, y0, '{0}'.format(ind))
  
  #plt.colorbar(pPlot)
  plt.show()

def demo_plot_segInds():
  
  fDir = '/data02/cases/2006/cfsr_anl/seg/'
  fMetr = fDir+'fields_pgbhnl.gdas.2006072012.grb2.trop.nc.npz'
  fSeg = fDir+'seg_pgbhnl.gdas.2006072012.grb2.trop.nc.npz'
  
  data = np.load(fMetr)
  lat = data['lat']; lon = data['lon']; nLat = len(lat); nLon = len(lon)
  data.close()
  
  data_seg = np.load(fSeg)
  cell2Site = data_seg['cell2Site'][:]
  cellIsMin = data_seg['cellIsMin'][:]
  cellIsMax = data_seg['cellIsMax'][:]
  data_seg.close()
  
  plot_segment_inds(lat, lon, cell2Site, cellIsMin, cellIsMax)

def plot_segmentSites(lat, lon, cell2Site, isMin, isMax):
  
  nLat = len(lat); nLon = len(lon); nCells = len(cell2Site)
  
  r2d = 180./np.pi
  m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  
  plt.figure()
  m.drawcoastlines()
  m.drawmapboundary()
  
  if (True):
    #sites = np.unique(cell2Site[isMin])
    sites = np.array([i for i in xrange(nCells) if cell2Site[i]==i])
    nSites = len(sites); print nSites
    iLats, iLons = index_1dTo2d(sites, nLon)
    for iSite in xrange(nSites):
      iLat = iLats[iSite]; iLon = iLons[iSite]
      x0, y0 = m(lon[iLon]*r2d, lat[iLat]*r2d)
      ind = sites[iSite]
      markerString = '${0}$'.format(ind)
      #markerString = 'o'
      m.scatter(x0, y0, s=500, marker=markerString)
  
  if (False):
    sites = np.unique(cell2Site[isMax])
    nSites = len(sites)
    iLats, iLons = index_1dTo2d(sites, nLon)
    for iSite in xrange(nSites):
      iLat = iLats[iSite]; iLon = iLons[iSite]
      x0 = x[iLon, iLat]; y0 = y[iLon, iLat]; ind = sites[iSite]
      markerString = '${0}$'.format(ind)
      m.scatter(x0, y0, s=100, marker=markerString)
  
  plt.show()
  
def demo_plot_segSites():
  
  #fDir = '/data02/cases/2006/cfsr_anl/seg/'
  #fMetr = fDir+'fields_pgbhnl.gdas.2006072012.grb2.trop.nc.npz'
  #fSeg = fDir+'seg_pgbhnl.gdas.2006072012.grb2.trop.nc.npz'
  fDir = '/data02/cases/2006/eraI/seg/'
  fMetr = fDir+'fields_2006-07-20_00.npz'
  fSeg = fDir+'seg_2006-07-20_00.npz'
  
  data = np.load(fMetr)
  lat = data['lat']; lon = data['lon']; nLat = len(lat); nLon = len(lon)
  data.close()
  
  data_seg = np.load(fSeg)
  cell2Site = data_seg['cell2Site'][:]
  cellIsMin = data_seg['cellIsMin'][:]
  cellIsMax = data_seg['cellIsMax'][:]
  data_seg.close()
  
  plot_segmentSites(lat, lon, cell2Site, cellIsMin, cellIsMax)

def demo_segment():
  
  rEarth = 6370.e3; dFilter = 300.e3 #radius in same units as r of disk for min basin
  fDirSave = '/data02/cases/2014/segment/2006/'
  
  #fnames = sorted(glob.glob('/data02/cases/2014/gfs_4_20140101_0000_000.n*'))
  #fnames = sorted(glob.glob('/data02/cases/2014/gfs_4_20140102_1200_000.n*'))
  #fnames = sorted(glob.glob('/data02/cases/2013/gfs_*.nc'))
  fnames = sorted(glob.glob('/arctic1/nick/cases/cfsr/2006-anl/pgbhnl.gdas*.nc'), key=os.path.getmtime)
  #fnames = sorted(glob.glob('/data02/cases/2014/gfs_4_20140101_00*.nc'))
  print fnames
  for iFile in xrange(len(fnames)):
    #pull in info from file and derive necessary -----------------
    print "Reading data from file"; t0 = time.clock()
    fname = fnames[iFile]; 
    data = netCDF4.Dataset(fname,'r')
    
    isMissing = get_missingCells_file(data)
    lat, lon, u, v, tmp, press = get_segmentVars_file(data)
    
    data.close()
    t1 = time.clock(); print "Done reading data from file in {0} s".format(t1-t0)
    nLat = len(lat); nLon = len(lon); nCells = nLat*nLon
    
    #specify region for segmentation, ----------------------------
    #including halo so make sure nbr values are decent
    latThreshHalo = 43.*np.pi/180.
    latThresh = 45.*np.pi/180.
    inRegionHalo = np.zeros((nLat,nLon), dtype=int); inRegionHalo[lat>latThreshHalo,:] = 1
    inRegion = np.zeros((nLat,nLon), dtype=int); inRegion[lat>latThresh,:] = 1
    
    #fill missing values in region of segmentation --------------
    t3 = time.clock();
    u = fill_missingVals_region(u, nLat, nLon, isMissing, inRegionHalo)
    v = fill_missingVals_region(v, nLat, nLon, isMissing, inRegionHalo)
    tmp = fill_missingVals_region(tmp, nLat, nLon, isMissing, inRegionHalo)
    press = fill_missingVals_region(press, nLat, nLon, isMissing, inRegionHalo)
    
    t4 = time.clock(); print "Time for filling missing values: {0} s".format(t4-t3)
    
    theta = calc_potentialTemperature(tmp, press)
    vort = calc_vertVorticity_ll(u, v, nLat, nLon, lat, rEarth)
    t5 = time.clock(); print "Time for calculating theta, vVort: {0} s".format(t5-t4)
    #vort[inRegion==False] = 0.
    #plot_var_ll(lat, lon, vort, [-5.e-4, 5.e-4])
    #theta[inRegionHalo==False] = 320.; plot_var_ll(lat, lon, theta,[])# [280, 360])
    #print theta[5,:]; #lat[0] = pi/2
    
    #segment ------------------------------------------
    theta = flatten_2dTo1d(theta, nLat, nLon)
    vort = flatten_2dTo1d(vort, nLat, nLon)
    inRegion = flatten_2dTo1d(inRegion, nLat, nLon)
    cell2Site, cellIsMin, cellIsMax = segment_high_low_watershed_region(lat, lon, theta, vort, rEarth, inRegion, dFilter)
    t6 = time.clock(); print "Time for segmenting watersheds: {0} s".format(t6-t5)
    
    #save segmented result
    fInfo = fname.split('/')[-1]
    fNameSave = 'seg_'+fInfo+'.npz'
    f = fDirSave+fNameSave
    print "Saving cell2Site to file "+f
    np.savez(f, cell2Site=cell2Site, cellIsMin=cellIsMin, cellIsMax=cellIsMax)
    fNameSave = 'fields_'+fInfo+'.npz'
    f = fDirSave+fNameSave
    print "Saving cell2Site data to file "+f
    np.savez(f, lat=lat, lon=lon, u=u, v=v, theta=theta, inRegion=inRegion)
    
    #plot segmented result ----------
    #if (True):
    if (False):
      print "Plotting segmentation"
      thetaFlat = theta; #flatten_2dTo1d(theta, nLat, nLon)
      #print thetaFlat[cellIsMin>0]; print thetaFlat[cellIsMax>0]
      isMin = unflatten_1dTo2d(cellIsMin, nLat, nLon)
      isMax = unflatten_1dTo2d(cellIsMax, nLat, nLon)
      vals = thetaFlat[cell2Site];
      vals = unflatten_1dTo2d(vals, nLat, nLon)
      plot_segment(lat, lon, vals, isMin, isMax)

def eraTimeToCalendarTime(hrs):
  #in the ERA file, time is stored as "hours since 1900-01-01 00:00:0.0"
  #we'll convert that to a datetime object and return a nice looking string
  
  tBase = dt.datetime(1900, 1, 1, 0)
  #note that TypeError: unsupported type for timedelta hours component: numpy.int32
  tNew = tBase + dt.timedelta(hours=hrs)
  tTuple = dt.datetime.timetuple(tNew);
  s = time.strftime('%Y-%m-%d_%H', tTuple)
  return s

def demo_segment_era():
  #ERA-I data ordered as netcdf from ecmwf comes as all times in 1 file with
  #missing data filled on the variable read...ie there's no missing data array.
  #so, we'll just make a different function for dealing with this rather than making
  #a more abstract (ie difficult to read) version for general use.
  
  rEarth = 6370.e3; dFilter = 300.e3 #radius in same units as r of disk for min basin
  fDirSave = '/data02/cases/2006/eraI/seg/'
  
  fDirData = '/data02/cases/2006/eraI/pv/'
  fnames = sorted(glob.glob(fDirData+'eraI_theta-u-v_2pvu_2006-07-20*.nc'), key=os.path.getmtime)
  print fnames
  for iFile in xrange(len(fnames)):
    #gather persistent info like mesh, times,... ----------------------
    fname = fnames[iFile]; 
    data = netCDF4.Dataset(fname,'r')
    
    times = data.variables['time'][:]; nTimes = len(times)
    d2r = np.pi/180.; 
    lat = data.variables['latitude'][:]*d2r; lon = data.variables['longitude'][:]*d2r
    nLat = len(lat); nLon = len(lon)
    
    #specify region for segmentation, ----------------------------
    #including halo so make sure nbr values are decent
    latThreshHalo = 43.*np.pi/180.
    latThresh = 45.*np.pi/180.
    inRegionHalo = np.zeros((nLat,nLon), dtype=int); inRegionHalo[lat>latThreshHalo,:] = 1
    inRegion = np.zeros((nLat,nLon), dtype=int); inRegion[lat>latThresh,:] = 1
    
    #loop over individual times ------------------------------
    for iTime in xrange(nTimes):
      theta = data.variables['pt'][iTime,:,:]
      u = data.variables['u'][iTime,:,:]; v = data.variables['v'][iTime,:,:]
      vort = calc_vertVorticity_ll(u, v, nLat, nLon, lat, rEarth)
    
      #segment ------------------------------------------
      theta = flatten_2dTo1d(theta, nLat, nLon)
      vort = flatten_2dTo1d(vort, nLat, nLon)
      inRegion = flatten_2dTo1d(inRegion, nLat, nLon)
      cell2Site, cellIsMin, cellIsMax = segment_high_low_watershed_region(lat, lon, theta, vort, rEarth, inRegion, dFilter)
    
      #save segmented result
      t0 = int(times[iTime])
      fInfo = eraTimeToCalendarTime(t0)
      fNameSave = 'seg_'+fInfo+'.npz'
      f = fDirSave+fNameSave
      print "Saving cell2Site to file "+f
      np.savez(f, cell2Site=cell2Site, cellIsMin=cellIsMin, cellIsMax=cellIsMax)
      fNameSave = 'fields_'+fInfo+'.npz'
      f = fDirSave+fNameSave
      print "Saving cell2Site data to file "+f
      np.savez(f, lat=lat, lon=lon, u=u, v=v, theta=theta, inRegion=inRegion)
    
      #plot segmented result ----------
      #if (True):
      if (False):
        print "Plotting segmentation"
        thetaFlat = theta; #flatten_2dTo1d(theta, nLat, nLon)
        #print thetaFlat[cellIsMin>0]; print thetaFlat[cellIsMax>0]
        isMin = unflatten_1dTo2d(cellIsMin, nLat, nLon)
        isMax = unflatten_1dTo2d(cellIsMax, nLat, nLon)
        vals = thetaFlat[cell2Site];
        vals = unflatten_1dTo2d(vals, nLat, nLon)
        plot_segment(lat, lon, vals, isMin, isMax)
    #looped over all times in this file
    data.close()

def demo_segment_steven():
  #Seeing what's happening with Steven's NH file
  #Latitude is reversed
  
  rEarth = 6370.e3; dFilter = 300.e3 #radius in same units as r of disk for min basin
  fDirSave = '/data02/cases/test_segment/stevenCase/'
  
  fDirData = '/home/scavallo/'
  fnames = sorted(glob.glob(fDirData+'erainterim_pv_201309*.nc'), key=os.path.getmtime)
  print fnames
  for iFile in xrange(len(fnames)):
    #gather persistent info like mesh, times,... ----------------------
    fname = fnames[iFile]; 
    data = netCDF4.Dataset(fname,'r')
    
    nTimes = len(data.dimensions['time']); times = range(nTimes)
    d2r = np.pi/180.; 
    lat = data.variables['xlat'][:,0]*d2r; lon = data.variables['xlong'][0,:]*d2r
    nLat = len(lat); nLon = len(lon)
    lat = lat[::-1] # reverse so north pole is lat[0]
    
    iLevel = 3; #for 2 pvu
    print "PVU level: ", data.variables['levels'][iLevel]
    
    #specify region for segmentation, ----------------------------
    #including halo so make sure nbr values are decent
    latThreshHalo = 43.*np.pi/180.
    latThresh = 45.*np.pi/180.
    inRegionHalo = np.zeros((nLat,nLon), dtype=int); inRegionHalo[lat>latThreshHalo,:] = 1
    inRegion = np.zeros((nLat,nLon), dtype=int); inRegion[lat>latThresh,:] = 1
    
    #loop over individual times ------------------------------
    for iTime in xrange(nTimes):
      theta = data.variables['th'][iTime,iLevel,:,:]
      u = data.variables['u'][iTime,iLevel,:,:]; v = data.variables['v'][iTime,iLevel,:,:]
      
      #reorient latitudes
      theta = theta[::-1,:]; u = u[::-1,:]; v = v[::-1,:]
      
      #fill missing values
      isMissing = np.isnan(theta)
      u = fill_missingVals_region(u, nLat, nLon, isMissing, inRegionHalo)
      v = fill_missingVals_region(v, nLat, nLon, isMissing, inRegionHalo)
      theta = fill_missingVals_region(theta, nLat, nLon, isMissing, inRegionHalo)
      
      #calc vorticity
      vort = calc_vertVorticity_ll(u, v, nLat, nLon, lat, rEarth)
    
      #segment ------------------------------------------
      theta = flatten_2dTo1d(theta, nLat, nLon)
      vort = flatten_2dTo1d(vort, nLat, nLon)
      inRegion = flatten_2dTo1d(inRegion, nLat, nLon)
      cell2Site, cellIsMin, cellIsMax = segment_high_low_watershed_region(lat, lon, theta, vort, rEarth, inRegion, dFilter)
    
      #save segmented result
      t0 = int(times[iTime])
      fInfo = eraTimeToCalendarTime(t0)
      fNameSave = 'seg_'+fInfo+'.npz'
      f = fDirSave+fNameSave
      print "Saving cell2Site to file "+f
      np.savez(f, cell2Site=cell2Site, cellIsMin=cellIsMin, cellIsMax=cellIsMax)
      fNameSave = 'fields_'+fInfo+'.npz'
      f = fDirSave+fNameSave
      print "Saving cell2Site data to file "+f
      np.savez(f, lat=lat, lon=lon, u=u, v=v, theta=theta, inRegion=inRegion)
    
      #plot segmented result ----------
      if (True):
      #if (False):
        print "Plotting segmentation"
        thetaFlat = theta; #flatten_2dTo1d(theta, nLat, nLon)
        #print thetaFlat[cellIsMin>0]; print thetaFlat[cellIsMax>0]
        isMin = unflatten_1dTo2d(cellIsMin, nLat, nLon)
        isMax = unflatten_1dTo2d(cellIsMax, nLat, nLon)
        vals = thetaFlat[cell2Site];
        vals = unflatten_1dTo2d(vals, nLat, nLon)
        #plot_segment(lat, lon, vals, isMin, isMax)
        fNameSave = fDirSave+fNameSave+'.png'
        print "Saving segmentation plot to: "+fNameSave
        plot_segment_save(fNameSave, lat, lon, vals, isMin, isMax)
    #looped over all times in this file
    data.close()

if __name__ == '__main__':
  #example_segment()
  #demo_plot_segInds()
  #demo_plot_segSites()
  #demo_segment()
  #demo_segment_era()
  demo_segment_steven()
