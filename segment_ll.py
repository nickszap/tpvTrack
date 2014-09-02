import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#watershed has a few options for implementation:
#-for every cell, walk down steepest gradient to the basin
#-every local min is a voronoi map with distance as the sum of ?
#along path to min/voronoi site. benefit of this is ability to remove voronoi sites for iterative merging.
#could remove based on min size of watershed,...

#active contours to identify regions:
#-minimize: internal energy + external energy + curveCost(eg smoothness)
#can create logically square say nearest nbr resample of area N of 75N and ship over to image processing+

def index_2dTo1d(iLat, iLon, nLon):
  return iLat*nLon+iLon
  
def index_1dTo2d(ind, nLon):
  iLon = ind%nLon
  iLat = (ind-ind%nLon)/nLon
  return (iLat, iLon)

def flatten_2dTo1d(vals, nLat, nLon): #untested!
  #vals[lat][lon] goes to vals[iLat*nLon+iLon]
  
  valsOut = np.reshape(vals, nLat*nLon)
  '''
  gives same values as:
  for iLat  in xrange(nLat):
    for iLon in xrange(nLon):
      ind1d = index_2dTo1d(iLat, iLon, nLon)
      valsOut[ind1d] = vals[iLat, iLon]
  plus, retains mask
  '''
  return valsOut

def nbrInds_ll(iLat, iLon, nLat, nLon):
  #return list of 8-conn nbr indices: [(latNbr1, lonNbr1), (latNbr2, lonNbr2),...]
  #always have east, west neighbors. not north/south at respective pole
  
  iWest = (iLon-1)%nLon
  iEast = (iLon+1)%nLon
  iSouth = iLat-1
  iNorth = iLat+1
  
  haveSouth = iSouth>0
  haveNorth = iNorth<nLat
  
  nbrs = [(iLat, iWest), (iLat, iEast)]
  if (haveSouth):
    nbrs.extend([(iSouth, iWest), (iSouth, iLon), (iSouth, iEast)])
    
  if (haveNorth):
    nbrs.extend([(iNorth, iWest), (iNorth, iLon), (iNorth, iEast)])

  return nbrs

def calcDistance_dLat(r, dAngle):
  #return distance given angle in radians
  s = np.abs(r*dAngle)
  return s

def calcDistance_dLon(r, dAngle, lat):
  rLat = r*np.cos(lat)
  s = np.abs(rLat*dAngle)
  return s

def calcDistance_ll(r, dLat, dLon, lat):
  #"cartesian spherical" distance
  distLat = calcDistance_dLat(r, dLat)
  distLon = calcDistance_dLon(r, dLon, lat)
  return np.sqrt(distLat*distLat+distLon*distLon)

def gatherCells_region(iLat0, iLon0, nLat, nLon, latCell, lonCell, 
                       refDist, distRegion):
  #return list of cell indices within specified spatial region of given point  
  
  nbrsInRegion = []
  candidateNbrs = nbrInds_ll(iLat0, iLon0, nLat, nLon)
  while (len(candidateNbrs)>0):
    nbrInd = candidateNbrs.pop(); iLat = nbrInd[0]; iLon = nbrInd[1]
    dLat = latCell[iLat]-latCell[iLat0]
    dLon = lonCell[iLon]-lonCell[iLon0]
    
    d = calcDistance_ll(refDist, dLat, dLon, iLat0)
    if (d<distRegion):
      nbrsInRegion.append(nbrInd)
      nbrs = nbrInds_ll(iLat0, iLon0, nLat, nLon)
      for nbrInd in nbrs:
        if ((nbrInd not in candidateNbrs) and (nbrInd not in nbrsInRegion)):
          candidateNbrs.append(nbrInd)
  
  return nbrsInRegion
  
def nbrDist_ll(iLat0, iLon0, nbrInds, distLat, distLon):
  #"cartesian spherical" distance
  #difference in indices can be huge, eg 359-0 for 1 degree
  
  nNbrs = len(nbrInds)
  d = np.empty(nNbrs)
  for iNbr in xrange(nNbrs):
    #sign(0) = 0
    dLat = np.sign(nbrInds[iNbr][0]-iLat0); dLat *= distLat
    dLon = np.sign(nbrInds[iNbr][1]-iLon0); dLon *= distLon
    d[iNbr] = np.sqrt(dLat*dLat+dLon*dLon)
    
  return d

def find_minCells(vals, nLat, nLon):
  #return array[nCells] with 1 if cell is min

  isMin = np.zeros((nLat, nLon),dtype=int)
  for iLat  in xrange(nLat):
    for iLon in xrange(nLon):
      nbrs = nbrInds_ll(iLat, iLon, nLat, nLon)
      
      val0 = vals[iLat, iLon]
      valNbrs = vals[nbrs]
      if (False not in (val0<=valNbrs)): #is site if can't descend from it
        isMin[iCell] = 1

  nMin = np.sum(isMin)
  print "Number of local min: ", nMin
  return isMin

def find_minCells_region(vals, nLat, nLon, inRegion):
  #return array[nCells] with 1 if cell is min

  isMin = np.zeros((nLat, nLon),dtype=int)
  for iLat  in xrange(nLat):
    for iLon in xrange(nLon):
      if (inRegion[iLat, iLon]<1):
        continue
      nbrs = nbrInds_ll(iLat, iLon, nLat, nLon)
      
      val0 = vals[iLat, iLon]
      valNbrs = vals[nbrs]
      if (False not in (val0<=valNbrs)): #is site if can't descend from it
        isMin[iCell] = 1

  nMin = np.sum(isMin)
  print "Number of local min: ", nMin
  return isMin

#def watershed(vals, cellIsMin, cellsOnCell, nEdgesOnCell, edgesOnCell, dcEdge, nCells, dRegion, latCell, lonCell):
def watershed(vals, cellIsMin, nLat, nLon, refDist, distRegion, latCell, lonCell):
  #to make adding/deleting basins simple, follow gradient until reach a site.
  #map every cell to follow local steepest gradient. basins go to self.
  #filter basins so have to be a min within specified region (disk of radius dRegion)
  #return map of cell to basin
  
  indMap_2dTo1d = np.array(nLat*nLon).reshape(nLat,nLon) #iLat*nLon+iLon

  #min goes to self
  cell2Site = np.zeros((nLat, nLon),dtype=int)
  for iLat  in xrange(nLat):
    for iLon in xrange(nLon):
      if (cellIsMin[iLat, iLon]>0):
        ind = indMap_2dTo1d[iLat, iLon]
        cell2Site[iLat, iLon] = ind

  #get local steepest path
  dAngleLat = np.pi/nLat; dAngleLon = 2*np.pi/nLon
  #refDist = 6371.e3
  for iLat  in xrange(nLat):
    distLat = calcDistance_dLat(refDist, dAngleLat)
    distLon = calcDistance_dLon(refDist, dAngleLon, latCell[iLat])
    
    for iLon in xrange(nLon):
      if (cellIsMin[iLat, iLon]>0):
        continue
        
      val0 = vals[iLat, iLon]
      nbrs = nbrInds_ll(iLat, iLon, nLat, nLon)
      valNbrs = vals[nbrs]
      dNbrs = nbrDist_ll(iLat, iLon, nbrs, distLat, distLon)
      
      valNbrs = (valNbrs-val0)/dNbrs
      iNbr = np.argmin(valNbrs)
      
      latNbr = nbrs[iNbr][0]; lonNbr = nbrs[iNbr][1]
      ind = indMap_2dTo1d[latNbr, lonNbr]
      cell2Site[iLat, iLon] = ind
    
  #Filter local extrema by area to limit high (spatial) frequency "noise".
  #For multiple close mins, the smallest counts as min for that region.
  #dRegion = 300.e3 #radius in meters of disk of filtering region
  nRedirect = 0
  
  for iLat  in xrange(nLat):
    for iLon in xrange(nLon):
      if (cellIsMin[iLat, iLon]>0):
        indsRegion = gatherCells_region(iLat, iLon, nLat, nLon, latCell, lonCell, 
                                        refDist, distRegion)
        
        val0 = vals[iLat, iLon]
        valsRegion = vals[indsRegion]
        minInd = np.argmin(valsRegion); minVal = valsRegion[minInd]
        if (minVal<valsRegion):
          cellIsMin[iLat, iLon] = 0
          minInd = indsRegion[minInd]
          cell2Site[iLat, iLon] = indMap_2dTo1d[minInd]
          nRedirect = nRedirect+1
  print "Number of redirects for regional min: ", nRedirect
  
  #follow local steepest path (and any redirections from, say, regional thresholds) to site
  for iLat  in xrange(nLat):
    for iLon in xrange(nLon):
      nextCell = cell2Site[iLat, iLon]
      nextLat, nextLon = index_1dTo2d(nextCell, nLon)
      while (not cellIsMin[nextLat, nextLon]>0):
        nextCell = cell2Site[nextLat, nextLon]
        nextLat, nextLon = index_1dTo2d(nextCell, nLon)
      
      cell2Site[iLat, iLon] = nextCell
      
  return cell2Site

def watershed_region(vals, cellIsMin, cellsOnCell, nEdgesOnCell, edgesOnCell, dcEdge, nCells, dRegion, latCell, lonCell, inRegion):
  #to make adding/deleting basins simple, follow gradient until reach a site.
  #map every cell to follow local steepest gradient. basins go to self.
  #filter basins so have to be a min within specified region (disk of radius dRegion)
  #return map of cell to basin.
  
  '''
  #to adapt global watershed to region, make values outside of region huge so don't steepest descend that way
  bigVal = 1.e10
  vals = np.copy(valsIn)
  vals[inRegion<1] = bigVal
  '''
  
  cell2Site = -np.ones(nCells,dtype=int)
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
    nNbrs = nEdgesOnCell[iCell]
    nbrs = cellsOnCell[iCell,0:nNbrs]
    valNbrs = vals[nbrs]
    val0 = vals[iCell]

    #correspondence is towards minimum gradient.
    for iNbr in xrange(nNbrs):
      iEdge = edgesOnCell[iCell,iNbr]
      dx = dcEdge[iEdge]
      valNbrs[iNbr] = (valNbrs[iNbr]-val0)/dx

    iNbr = np.argmin(valNbrs)
    cell2Site[iCell] = nbrs[iNbr]
    
  #Filter local extrema by area to limit high (spatial) frequency "noise".
  #For multiple close mins, the smallest counts as min for that region.
  #dRegion = 300.e3 #radius in meters of disk of filtering region
  pt_ll = np.empty(2, dtype=float)
  nRedirect = 0
  for iCell in xrange(nCells):
    if (inRegion[iCell]<1):
      continue
    if (cellIsMin[iCell]>0):
      #see if cell is min in region, not just neighbors.
      #if not regional min, update cell2Site so local min goes to another basin
      pt_ll[0] = latCell[iCell]; pt_ll[1] = lonCell[iCell]
      cellsRegion = conn.gatherCells_radius(pt_ll, dRegion, iCell, cellsOnCell, nEdgesOnCell, latCell, lonCell)
      valsRegion = vals[cellsRegion]
      minInd = np.argmin(valsRegion);
      minCell = cellsRegion[minInd]; minRegion = valsRegion[minInd]
      if (minRegion<vals[iCell]): #have smaller value somewhere else in disk
        cellIsMin[iCell] = 0
        cell2Site[iCell] = minCell
        nRedirect = nRedirect+1
  print "Number of redirects for regional min: ", nRedirect
  
  #follow local steepest path (and any redirections from, say, regional thresholds) to site
  for iCell in xrange(nCells):
    if (inRegion[iCell]<1):
      continue
    nextCell = cell2Site[iCell]
    while (not cellIsMin[nextCell]>0):
      nextCell = cell2Site[nextCell]
      #print "Cell {0} going to {1}".format(iCell,nextCell)

    cell2Site[iCell] = nextCell

  return cell2Site

def plotSegment(lat, lon, var, isMarker):  
  #map = Basemap(projection='ortho',lon_0=100,lat_0=60, resolution='l')
  map = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  x,y = map(lon, lat)
  
  fig1 = plt.figure()
  map.drawcoastlines()
  map.drawmapboundary()
  pPlot = map.pcolor(x,y,var,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet) #cmap=plt.cm.hot_r) #,vmin=100,vmax=1000)
  
  xMarker = x[isMarker>0]
  yMarker = y[isMarker>0]
  map.scatter(xMarker,yMarker,marker="o")
  
  plt.colorbar(pPlot)
  #plt.show()
  return plt

def plotSegment_time(fnameMPAS, fnameSegment):
  
  data = netCDF4.Dataset(fnameMPAS,'r')
  r2d = 180./np.pi
  latCell = data.variables['latCell'][:]*r2d;
  lonCell = data.variables['lonCell'][:]*r2d;
  
  dataSeg = np.load(fnameSegment)
  cell2SiteHistory = dataSeg['cell2SiteHistory'][:]
  nTimes, nCells = cell2SiteHistory.shape
  
  for iTime in xrange(nTimes):
    tInd = 0+iTime
    cell2Site = cell2SiteHistory[iTime,:]
    
    isSite = cell2Site==range(nCells) #site goes to self
    #plotSegment(latCell, lonCell, cell2Site, isSite) #pcolor value is of cell index
    vals = data.variables['theta_pv'][tInd, :]
    for iCell in xrange(nCells):
      vals[iCell] = vals[cell2Site[iCell]]
    
    plotSegment(latCell,lonCell,vals,isSite) #pcolor value is of metr value
    
  plt.show()

def demo_3depv_basin_vtk():
  
  fnameMPAS = '/arctic1/nick/cases/epv-budget/testRun/x4.t.output.2006-08-01_00.00.00.nc'
  data = netCDF4.Dataset(fnameMPAS,'r')
  r2d = 180./np.pi
  #latCell = data.variables['latCell'][:]*r2d;
  #lonCell = data.variables['lonCell'][:]*r2d;
  
  fnameSegment = '/data01/tracks/algo/test/segment_cell2Site_conn.npz'
  dataSeg = np.load(fnameSegment)
  cell2SiteHistory = dataSeg['cell2SiteHistory'][:]
  nTimes, nCells = cell2SiteHistory.shape
  
  vortCells = [56355, 14136, 62441] #vortex min tracked over time
  
  for iTime in xrange(nTimes):
    #tInd = 0+iTime
    cell2Site = cell2SiteHistory[iTime,:]
    
    inBasin = cell2Site==vortCells[iTime]
    basinCells = np.array(np.arange(nCells))[inBasin]
    
    output_data.example_smallRegion_3dplane(data, basinCells)
    
  data.close()  

#filter field
def smooth_mean(cellsOnCell, nEdgesOnCell, valsIn, nCells):
  #val[c0] is average of neighbor values
  
  vals = np.empty_like (valsIn)
  #np.copyto(vals, valsIn)
  vals[:] = valsIn[:]
  
  for iCell in xrange(nCells):
    nNbrs = nEdgesOnCell[iCell]
    nbrs = cellsOnCell[iCell,0:nNbrs]
    valNbrs = valsIn[nbrs]
    val0 = valsIn[iCell]
    
    vals[iCell] = (nNbrs*np.mean(valNbrs)+val0)/(nNbrs+1)
    
  return vals

#reconstruction

def example_segment():
  #we want to associate features with cyclones and anticyclones.
  #we'll create basins for local minima and maxima separately so each cell
  #can be part of 2 basins. To avoid over-segmentation, smooth field first.
  #Then, we map the cell to the proper type of basin based on additional information (eg vorticity)

  ncfname = '/arctic1/nick/cases/v1.0/x4/longer/x4.kf.output.2006-08-15_00.00.00.nc'
  #ncfname = '/arctic1/nick/cases/2011/integrate/x7.kf.o3.output.2011-12-25_00.00.00.nc'
  data = netCDF4.Dataset(ncfname,'r')

  nCells = len(data.dimensions['nCells'])
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  edgesOnCell = data.variables['edgesOnCell'][:]-1;
  dcEdge = data.variables['dcEdge'][:]
  
  #field
  tInd = 0
  #tInd = 12
  epv_ht, theta_trop = compare.calc_height_theta_2PVU(data, tInd)
  #vort_trop, theta_trop = calc_vorticity_theta_2PVU(data, t0)
  #segment based on min or max
  theta_trop = -theta_trop
  
  
  #smooth
  nSmooth = 10
  for iSmooth in xrange(nSmooth):
    theta_trop = smooth_mean(cellsOnCell, nEdgesOnCell, theta_trop, nCells)
  
  #segment
  cellIsMin = find_minCells(cellsOnCell, nEdgesOnCell, theta_trop, nCells)
  cell2Site = watershed(theta_trop, cellIsMin, cellsOnCell, nEdgesOnCell, edgesOnCell, dcEdge, nCells)

  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];

  latCell *= 180./np.pi
  lonCell *= 180./np.pi

  plotSegment(latCell, lonCell, cell2Site, cellIsMin)
  
  data.close()

def segment_high_low_watershed(data, tInd):
  #get high and low basin seeds, associate cells to both high and low basins if not extrema.
  #to decide whether "really" part of high or low basin, we have options:
  #-(anti-)cyclonic for (high) low...is local vorticity noisy?
  #-closer theta value to maxima a la color scale grouping...huge min or max value now matters
  #-whether steeper gradient is to high or low
  #-physical distance
  #-concavity of surface a la last closed contour

  nCells = len(data.dimensions['nCells'])
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  edgesOnCell = data.variables['edgesOnCell'][:]-1;
  dcEdge = data.variables['dcEdge'][:]
  latCell = data.variables['latCell'][:]
  lonCell = data.variables['lonCell'][:]
  
  #field ---------------------------------
  vort_trop, theta_trop = calc_vorticity_theta_2PVU(data, tInd)
  print "Finished reading conn and fields"

  #smooth ------------------------------
  nSmooth = 0
  print "Spatial smoothing passes: ", nSmooth
  for iSmooth in xrange(nSmooth):
    theta_trop = smooth_mean(cellsOnCell, nEdgesOnCell, theta_trop, nCells)
  '''  
  printCells = [0,79440,114005]
  print theta_trop[printCells]
  '''
  #segment ------------------------------------------
  dRegion = 300.e3 #radius in meters of disk of filtering region
  
  #mins
  print "Finding minima"
  cellIsMin = find_minCells(cellsOnCell, nEdgesOnCell, theta_trop, nCells)
  cell2SiteMin = watershed(theta_trop, cellIsMin, cellsOnCell, nEdgesOnCell, edgesOnCell, dcEdge, nCells, dRegion, latCell, lonCell)  
  
  #maxs
  print "Finding maxima"
  theta_trop = -theta_trop
  cellIsMax = find_minCells(cellsOnCell, nEdgesOnCell, theta_trop, nCells)
  cell2SiteMax = watershed(theta_trop, cellIsMax, cellsOnCell, nEdgesOnCell, edgesOnCell, dcEdge, nCells, dRegion,latCell, lonCell)
  theta_trop = -theta_trop
  
  #"voting" procedure for low/high classification ------
  print "Associating to max or min"
  cell2Site = np.empty(nCells, dtype=int)
  
  classifyMethod = 0 #see if-checks below for corresponding number for method
  if (classifyMethod==0):
    #local vorticity
    print("Classifying by local vorticity")
    for iCell in xrange(nCells):
      #cyclonic depends on hemisphere
      if (cellIsMin[iCell]>0 or cellIsMax[iCell]>0):
        cell2Site[iCell] = iCell
      else:
        signHem = 1 #sign function is problematic since sign(0)=0
        if (latCell[iCell]<0): #lat=0 gets put in NH
          signHem = -signHem
        if (signHem*vort_trop[iCell]<0): #anticyclonic
          cell2Site[iCell] = cell2SiteMax[iCell]
        else: #0 or cyclonic
          cell2Site[iCell] = cell2SiteMin[iCell]
  
  elif (classifyMethod==1):      
    #curvature of contour:
    #this should be sensitive to surface noise/undulations so far from robust!!!
    print("Classifying by local curvature on low->cell->high")
    pt_ll = np.empty(2,dtype=float)
    for iCell in xrange(nCells):
      if (cellIsMin[iCell]>0 or cellIsMax[iCell]>0): #extrema goes to self
        cell2Site[iCell] = iCell
      else: #concave up goes to low, concave down to high
        #can form a path low->cell->high.
        #calculate local concavity at cell "along" path
        iCellLow = cell2SiteMin[iCell];
        iCellHigh = cell2SiteMax[iCell];
        
        nbrs = cellsOnCell[iCell, 0:nEdgesOnCell[iCell]]; nNbrs = len(nbrs)
        
        pt_ll[0] = latCell[iCellLow]; pt_ll[1] = lonCell[iCellLow]
        #cells2Low = conn.gatherCellsOnLine(iCell, pt_ll, cellsOnCell, nEdgesOnCell, latCell, lonCell)
        #iCellL = cells2Low[1]; 
        iCellL, dMin = conn.find_closestCellToPoint(pt_ll, nbrs, nNbrs, latCell, lonCell)
        
        pt_ll[0] = latCell[iCellHigh]; pt_ll[1] = lonCell[iCellHigh]
        #cells2High = conn.gatherCellsOnLine(iCell, pt_ll, cellsOnCell, nEdgesOnCell, latCell, lonCell)
        #iCellH = cells2High[1];
        iCellH, dMin = conn.find_closestCellToPoint(pt_ll, nbrs, nNbrs, latCell, lonCell)
        
        #"central" difference for 2nd derivative
        val0 = theta_trop[iCell]; valL = theta_trop[iCellL]; valH = theta_trop[iCellH];
        deriv2 = valL-2.*val0+valH #/dx^2 . constant spacing not insane for smooth mesh
        if (deriv2<0): #concave down
          cell2Site[iCell] = iCellHigh
        else: #inflection or concave up
          cell2Site[iCell] = iCellLow
  
  else:
    print "Error. Classification method not recognized/specified."
          
  return cell2Site

def segment_high_low_watershed_region(data, tInd, inRegion):
  #get high and low basin seeds, associate cells to both high and low basins if not extrema.
  #to decide whether "really" part of high or low basin, we have options:
  #-(anti-)cyclonic for (high) low...is local vorticity noisy?
  #-closer theta value to maxima a la color scale grouping...huge min or max value now matters
  #-whether steeper gradient is to high or low
  #-physical distance
  #-concavity of surface a la last closed contour
  
  #inRegion<1 for cells not in region

  nCells = len(data.dimensions['nCells'])
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  edgesOnCell = data.variables['edgesOnCell'][:]-1;
  dcEdge = data.variables['dcEdge'][:]
  latCell = data.variables['latCell'][:]
  lonCell = data.variables['lonCell'][:]
  
  #field ---------------------------------
  vort_trop, theta_trop = calc_vorticity_theta_2PVU(data, tInd)
  print "Finished reading conn and fields"

  #smooth ------------------------------
  nSmooth = 0
  print "Spatial smoothing passes: ", nSmooth
  for iSmooth in xrange(nSmooth):
    theta_trop = smooth_mean(cellsOnCell, nEdgesOnCell, theta_trop, nCells)
  '''  
  printCells = [0,79440,114005]
  print theta_trop[printCells]
  '''
  #segment ------------------------------------------
  dRegion = 300.e3 #radius in meters of disk of filtering region
  
  #mins
  #to adapt global watershed to region, make values outside of region huge so don't steepest descend that way
  bigVal = 1.e10
  vals = np.copy(theta_trop)
  vals[inRegion<1] = bigVal
  
  print "Finding minima"
  cellIsMin = find_minCells_region(cellsOnCell, nEdgesOnCell, vals, nCells, inRegion)
  cell2SiteMin = watershed_region(vals, cellIsMin, cellsOnCell, nEdgesOnCell, edgesOnCell, dcEdge, nCells, dRegion, latCell, lonCell, inRegion)
  
  #plot
  if (False):
    isSite = cell2SiteMin==range(nCells) #site goes to self
    r2d = 180./np.pi
    plotSegment(latCell*r2d, lonCell*r2d, cell2SiteMin, isSite)
  
  #maxs
  #adapt global watershed to region
  '''
  vals = np.copy(theta_trop)
  vals = -vals
  vals[inRegion<1] = bigVal
  '''
  vals[:] = -theta_trop[:]
  vals[inRegion<1] = bigVal
  
  print "Finding maxima"
  cellIsMax = find_minCells_region(cellsOnCell, nEdgesOnCell, vals, nCells, inRegion)
  cell2SiteMax = watershed_region(vals, cellIsMax, cellsOnCell, nEdgesOnCell, edgesOnCell, dcEdge, nCells, dRegion, latCell, lonCell, inRegion)
  
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
      if (cellIsMin[iCell]>0 or cellIsMax[iCell]>0):
        cell2Site[iCell] = iCell
      else:
        signHem = 1 #sign function is problematic since sign(0)=0
        if (latCell[iCell]<0): #lat=0 gets put in NH
          signHem = -signHem
        if (signHem*vort_trop[iCell]<0): #anticyclonic
          cell2Site[iCell] = cell2SiteMax[iCell]
        else: #0 or cyclonic
          cell2Site[iCell] = cell2SiteMin[iCell]
  
  elif (classifyMethod==1):      
    #curvature of contour:
    #this should be sensitive to surface noise/undulations so far from robust!!!
    print("Classifying by local curvature on low->cell->high")
    pt_ll = np.empty(2,dtype=float)
    for iCell in xrange(nCells):
      if (cellIsMin[iCell]>0 or cellIsMax[iCell]>0): #extrema goes to self
        cell2Site[iCell] = iCell
      else: #concave up goes to low, concave down to high
        #can form a path low->cell->high.
        #calculate local concavity at cell "along" path
        iCellLow = cell2SiteMin[iCell];
        iCellHigh = cell2SiteMax[iCell];
        
        nbrs = cellsOnCell[iCell, 0:nEdgesOnCell[iCell]]; nNbrs = len(nbrs)
        
        pt_ll[0] = latCell[iCellLow]; pt_ll[1] = lonCell[iCellLow]
        #cells2Low = conn.gatherCellsOnLine(iCell, pt_ll, cellsOnCell, nEdgesOnCell, latCell, lonCell)
        #iCellL = cells2Low[1]; 
        iCellL, dMin = conn.find_closestCellToPoint(pt_ll, nbrs, nNbrs, latCell, lonCell)
        
        pt_ll[0] = latCell[iCellHigh]; pt_ll[1] = lonCell[iCellHigh]
        #cells2High = conn.gatherCellsOnLine(iCell, pt_ll, cellsOnCell, nEdgesOnCell, latCell, lonCell)
        #iCellH = cells2High[1];
        iCellH, dMin = conn.find_closestCellToPoint(pt_ll, nbrs, nNbrs, latCell, lonCell)
        
        #"central" difference for 2nd derivative
        val0 = theta_trop[iCell]; valL = theta_trop[iCellL]; valH = theta_trop[iCellH];
        deriv2 = valL-2.*val0+valH #/dx^2 . constant spacing not insane for smooth mesh
        if (deriv2<0): #concave down
          cell2Site[iCell] = iCellHigh
        else: #inflection or concave up
          cell2Site[iCell] = iCellLow
  
  else:
    print "Error. Classification method not recognized/specified."
          
  return (cell2Site, cellIsMin, cellIsMax)

def demo_segment_region():
  #example for segmentation of the arctic
  ncfname = '/arctic1/nick/cases/epv-budget/testRun/x4.t.output.2006-08-01_00.00.00.nc'
  data = netCDF4.Dataset(ncfname,'r')
  
  nCells = len(data.dimensions['nCells'])
  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];
  latCell *= 180./np.pi
  lonCell *= 180./np.pi
  
  tInd = 0
  
  #region for segmentation
  latThresh = 60.
  inRegion = np.ones(nCells, dtype=int)
  inRegion[latCell<latThresh]=0
  
  #segment
  cell2Site = segment_high_low_watershed_region(data, tInd, inRegion)[0]
  
  if (True):
    print "Saving cell2Site to file"
    np.savez('segment_cell2Site', cell2Site=cell2Site)
  
  #plot
  if (True):
    print "Plotting segmentation"
    isSite = cell2Site==range(nCells) #site goes to self
    #plotSegment(latCell, lonCell, cell2Site, isSite) #pcolor value is of cell index
    vals = data.variables['theta_pv'][tInd, :]
    for iCell in xrange(nCells):
      vals[iCell] = vals[cell2Site[iCell]]
    
    plotSegment(latCell,lonCell,vals,isSite) #pcolor value is of metr value
  
  data.close()

def trackFeature_time():
  #to track a given high/low over time:
  #-at every time, have watershed defined cells as highs and lows
  #-each high/low is expected to be in a specified area at a future time
  #-if a high/low exists in the specified area at future time, the corresponding high/low is the one most similar
  
  #the specified area at future time can be a function of:
  #-a ball of specified distance, extrapolation of past feature trajectory,...
  
  #similarity can be some cost function combination of:
  #-area, min value, northernmost,...
  #should feature be connected???
  
  saveDir = '/data01/tracks/algo/test/'
  
  ncfname = '/arctic1/nick/cases/epv-budget/testRun/x4.t.output.2006-08-01_00.00.00.nc'
  data = netCDF4.Dataset(ncfname,'r')
  
  nTimes = len(data.dimensions['Time'])
  nCells = len(data.dimensions['nCells'])
  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  
  #region for segmentation ------------------
  latThreshd = 60.
  latThresh = latThreshd*np.pi/180.
  inRegion = np.ones(nCells, dtype=int)
  inRegion[latCell<latThresh]=0
  
  #segment --------------------------
  print "Watershed'ing over time with latThresh={0}".format(latThreshd)
  time0 = 0; timeEndp1 = 3; nTimes = timeEndp1-time0
  cell2SiteHistory = np.empty((nTimes, nCells), dtype=int)
  cellIsMinHistory = np.empty((nTimes, nCells), dtype=int)
  cellIsMaxHistory = np.empty((nTimes, nCells), dtype=int)
  
  if (True):
  #if (False):
    for tInd in xrange(time0,timeEndp1):
      cell2Site, cellIsMin, cellIsMax = segment_high_low_watershed_region(data, tInd, inRegion)
          
      if (True):
        #keep only connected regions by setting non-connected to -1
        for iCell in xrange(nCells):
          if (cell2Site[iCell]==iCell):
            connCells = conn.gatherConnectedCells(cell2Site, iCell, cellsOnCell, nEdgesOnCell)
            cell2Site[cell2Site==iCell] = -1
            cell2Site[connCells] = iCell
            
      cell2SiteHistory[tInd-time0,:]  = cell2Site[:]
      cellIsMinHistory[tInd-time0,:]  = cellIsMin[:]
      cellIsMaxHistory[tInd-time0,:]  = cellIsMax[:]
    
  if (True):
  #if (False):
    #fName = saveDir+'segment_cell2Site'
    fName = saveDir+'segment_cell2Site_conn'
    print "Saving cell2Site to file "+fName+'.npz'
    np.savez(fName, cell2SiteHistory=cell2SiteHistory, cellIsMinHistory=cellIsMinHistory, cellIsMaxHistory=cellIsMaxHistory)
  
  #track ----------------------------------
  vortTrack = [] #stores cells containing extremum of feature over time
  
  #load in segmentation history if not running watershed and save above
  if (True):
    #fSeg = saveDir+'segment_cell2Site'+'.npz'
    fSeg = saveDir+'segment_cell2Site_conn'+'.npz'
    f = np.load(fSeg)
    cell2SiteHistory = f['cell2SiteHistory'][:]
    cellIsMinHistory = f['cellIsMinHistory'][:]
    cellIsMaxHistory = f['cellIsMaxHistory'][:]
    
  #llVort = np.array([87., 320.])*np.pi/180. #cell in feature at this loc at initial time
  llVort = np.array([78.5, 3.])*np.pi/180.
  isExtr = cellIsMinHistory #change to ..IsMax... for tracking anticyclones
  
  #at initial time, find the vortex
  radiusSearchRegion = 500.e3 #region around lat/lon point for candidate features
  tInd = 0; iCellVortex = 0 #initial guess
  iCellVortex = conn.findOwner_horizNbrs_latLon(llVort, iCellVortex, latCell, lonCell, nEdgesOnCell, cellsOnCell)
  
  cell2Site0 = cell2SiteHistory[tInd]
  vortSite0 = cell2Site0[iCellVortex]; vortTrack.append(vortSite0);
  
  #initial vortex properties
  #extreme val ----
  val0 = data.variables['theta_pv'][time0+tInd,vortSite0]
  #region area -----
  areaCell = data.variables['areaCell'][:]
  isVortCell = cell2Site0==vortSite0
  area0 = np.sum(areaCell[isVortCell]/1.e6); #in km^2
  print "Area equivalent radius, km:", np.sqrt(area0/np.pi)
  #velocity ----
  #since didn't output 2pvu vel use 200hPa
  uZonal = data.variables['uzonal_200hPa'][time0+tInd,:]
  uMerid = data.variables['umeridional_200hPa'][time0+tInd,:]
  uZonal_vort = np.mean(uZonal[isVortCell]); uMerid_vort = np.mean(uMerid[isVortCell]) #m/s
  print "velocity d(lon,lat)/dt in m/s: ", uZonal_vort, uMerid_vort
  #convert to degrees/hr, change in long depends on radius(latitude)
  latMean = np.mean(latCell[isVortCell])
  uv = velocityToLatLon(uZonal_vort, uMerid_vort, latMean) #radians/hr
  print "velocity d(lon/lat)/dt in deg/hr: ", uv[0]*180/np.pi, uv[1]*180/np.pi
  
  #at future times, track the found vortex
  pt_ll = np.empty(2,dtype=float)
  deltaTime = 6.; print "Time step is {0} hours: ".format(deltaTime) #since vel in radians/hour
  for tInd in xrange(1,nTimes):
    #gather cells within search region of new vortex
    #advect search vortex
    pt_ll[0] = latCell[vortSite0]; pt_ll[1] = lonCell[vortSite0];
    pt_ll[0] += uv[1]*deltaTime; pt_ll[1] += uv[0]*deltaTime;
    if (pt_ll[0]>np.pi/2.):
      pt_ll[0] = np.pi/2.
    cellsInRegion = conn.gatherCells_radius(pt_ll, radiusSearchRegion, vortSite0, cellsOnCell, nEdgesOnCell, latCell, lonCell)
    
    #get a set of candidate lows or highs within region
    #all sites within region
    cell2Site1 = cell2SiteHistory[tInd]
    lhSites = np.unique(cell2Site1[cellsInRegion])
    #sites that are just lows or highs
    lSites = np.array([iCell for iCell in lhSites if isExtr[tInd,iCell]>0])
    
    #if more than one candidate, choose the "most similar" one
    nCandidates = len(lSites); print "{0} candidates for feature at time {1}: {2}".format(nCandidates, tInd, lSites)
    if (nCandidates<1):
      print "Stopping tracking since no candidates found"
      break
    if (nCandidates==1):
      vortSite1 = lSites[0]
    else:
      #choose vortex with closest theta_2pvu value
      vals1 = data.variables['theta_pv'][time0+tInd,lSites]
      indSite = np.argmin(np.abs(vals1-val0))
      vortSite1 = lSites[indSite]
      
    isVortCell = cell2Site1==vortSite1
    area1 = np.sum(areaCell[isVortCell]/1.e6); #in km^2
    print "Area equivalent radius, km:", np.sqrt(area1/np.pi)  
    
    #update vortex info for searching next time
    vortSite0 = vortSite1
    vortTrack.append(vortSite0)
    val0 = data.variables['theta_pv'][time0+tInd,vortSite0]
    
    uZonal = data.variables['uzonal_200hPa'][time0+tInd,:]
    uMerid = data.variables['umeridional_200hPa'][time0+tInd,:]
    uZonal_vort = np.mean(uZonal[isVortCell]); uMerid_vort = np.mean(uMerid[isVortCell]) #m/s
    #print "velocity d(lon,lat)/dt in m/s: ", uZonal_vort, uMerid_vort
    #convert to degrees/hr, change in long depends on radius(latitude)
    latMean = np.mean(latCell[isVortCell])
    uv = velocityToLatLon(uZonal_vort, uMerid_vort, latMean) #radians/hr
    print "velocity d(lon/lat)/dt in deg/hr: ", uv[0]*180/np.pi, uv[1]*180/np.pi
  
  print "\nFeature tracking summary:\n", "cells:", vortTrack
  print "lat:", latCell[vortTrack]*180./np.pi
  print "lon:", lonCell[vortTrack]*180./np.pi
  return vortTrack

def velocityToLatLon(uZonal, uMerid, lat):
  #convert u,v in m/s to radians/hr, change in long depends on radius(latitude)
  #2pi radians in 1 circumference
  #lenDegree = 2.*np.pi*6371.e3/360. #m/degree
  lenDegree = 6371.e3 #m/radian
  lenLat = lenDegree*np.cos(lat); lenLat = np.max((1.,lenLat)) #in case at pole
  degMerid = uMerid*3600./lenDegree
  degZonal = uZonal*3600./lenLat
  
  return np.array([degZonal, degMerid])
  
def calc_vorticity_theta_2PVU(data, t0):
  #vertical height of sign(lat)*2PVU surface.
  #return height of cell center in column using interpolation from top

  pvuTrop = 2.0

  nCells = len(data.dimensions['nCells'])
  nLevels = len(data.dimensions['nVertLevels'])
  epv = data.variables['ertel_pv'][t0,:,:]
  latCell = data.variables['latCell'][:]

  theta = data.variables['theta'][t0,:,:]
  vortVertex = data.variables['vorticity'][t0,:,:]
  nVertexOnCell = data.variables['nEdgesOnCell'][:]
  vertexOnCell = data.variables['verticesOnCell'][:,:]-1

  vort_trop = np.empty(nCells, dtype='float')
  theta_trop = np.empty(nCells,dtype='float')
  for iCell in xrange(nCells):
    pvuVal = pvuTrop
    if (latCell[iCell]<0):
      pvuVal = -pvuTrop

    #Considering top and bottom boundaries + noise/artifacts, question of
    #which levels we want to count for interpolation
    lev0 = 1; levTop = nLevels-1
    interpLevs = range(lev0,levTop); nInterpLevs = len(interpLevs)
    (l,dl) = output_data.calcIndexOfValue(pvuVal,epv[iCell,interpLevs], nInterpLevs)
    #print "Cell {0} has l,dl = {1},{2}".format(i, l, dl)
    theta_trop[iCell] = output_data.calcValueOfIndex(l,dl,theta[iCell,interpLevs])
    
    if (l<0): #missing val
      vort_trop[iCell] = output_data.missingVal
    else:
      #average vorticity to cell centers and interpolate to trop
      vortl = form_vertVorticityCell(iCell, lev0+l, vortVertex, vertexOnCell, nVertexOnCell)
      vorth = form_vertVorticityCell(iCell, lev0+l+1, vortVertex, vertexOnCell, nVertexOnCell)

      dvort_dl = vorth-vortl;
      vort_trop[iCell] = vortl+dvort_dl*dl
    
  return (vort_trop, theta_trop)


def form_vertVorticityCell(iCell, iLevel, vortVertex, vertexOnCell, nVertexOnCell):
  #approx way of getting vorticity at cell center from voronoi cell vertices

  nVerts = nVertexOnCell[iCell]
  verts = vertexOnCell[iCell,0:nVerts]
  vortVerts = vortVertex[verts,iLevel]
  vortMean = np.mean(vortVerts)
  return vortMean

'''
real(kind=RKIND) function calc_verticalVorticity_cell(c0, level, nVerticesOnCell, verticesOnCell, cellsOnVertex, &
                                                         kiteAreasOnVertex, areaCell, vVortVertex)
      !area weighted average of vorticity at vertices (really midpts of faces) to cell center for the specified cell
      !
      implicit none

      real(kind=RKIND), intent(in) :: areaCell
      integer, intent(in) :: c0, level, nVerticesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell, cellsOnVertex
      real(kind=RKIND), dimension(:,:), intent(in) :: kiteAreasOnVertex, vVortVertex

      real(kind=RKIND) :: vVortCell
      integer :: i, iVertex, cellIndOnVertex

      vVortCell = 0.0
      do i = 1,nVerticesOnCell
         iVertex = verticesOnCell(i,c0)
         cellIndOnVertex = elementIndexInArray(c0, cellsOnVertex(:,iVertex), 3)
         vVortCell = vVortCell + kiteAreasOnVertex(cellIndOnVertex, iVertex)*vVortVertex(level, iVertex)/areaCell
      end do

      calc_verticalVorticity_cell = vVortCell
   end function calc_verticalVorticity_cell
'''

if __name__ == '__main__':
  #example_segment()
  ncfname = '/arctic1/nick/cases/vduda/x4.t.output.2006-08-15_00.00.00.nc'
  #ncfname = '/data02/cases/v1.0/x4/longer/longer/x4.kf.output.2006-08-15_00.00.00.nc'
  #ncfname = '/arctic1/nick/cases/v1.0/x4/longer/x4.kf.output.2006-08-15_00.00.00.nc'
  #ncfname = '/arctic1/nick/cases/2011/integrate/x7.kf.o3.output.2011-12-25_00.00.00.nc'
  data = netCDF4.Dataset(ncfname,'r')
  
  tInd = 4
  
  cell2Site = segment_high_low_watershed(data, tInd)
  
  nCells = len(data.dimensions['nCells'])
  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];

  latCell *= 180./np.pi
  lonCell *= 180./np.pi
  
  isSite = cell2Site==range(nCells) #site goes to self
  plotSegment(latCell, lonCell, cell2Site, isSite)
