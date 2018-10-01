import numpy as np
import netCDF4
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import helpers

#watershed has a few options for implementation:
#-for every cell, walk down steepest gradient to the basin

def find_basinBoundaries(cell2Site, cell0, mesh):
  """a basin boundary is a cell with a cell2Site pointing to another basin"""
  
  isBoundary = np.zeros(mesh.nCells, dtype=int)
  for cell in iter(cell0.copy()): #not that cell referring to cell0 changes cell0
    if (not cell.isInRegion()):
      continue
    
    iCell = cell.ind
    nbrs = cell.get_nbrInds()
    
    val0 = cell2Site[iCell]; valNbrs = cell2Site[nbrs]
    if ( np.any(val0 != valNbrs) ): #have nbr with different basin
      isBoundary[iCell] = 1
  
  return isBoundary

def find_minCells_region_flat(vals, cell0, mesh):
  """return array[nCells] with 1 if cell is min and cell in region"""
  
  isMin = np.zeros(mesh.nCells,dtype=int)
  for cell in iter(cell0.copy()): #not that cell referring to cell0 changes cell0
    if (not cell.isInRegion()):
      continue
    
    iCell = cell.ind
    nbrs = cell.get_nbrInds()
    #assume that every cell inRegion has >1 neighbor inRegion
    nbrsInRegion = nbrs[mesh.isIndsInRegion(nbrs)]
    
    val0 = vals[iCell]; valNbrs = vals[nbrsInRegion]
    #print val0, valNbrs
    if (np.all(val0<=valNbrs)): #is site if can't descend from it
      isMin[iCell] = 1

  nMin = np.sum(isMin)
  print "Number of local min: ", nMin
  return isMin

def watershed_region(vals, cellIsMin, cell0, mesh):
  """
  Map every cell to a basin following local steepest gradient. Regional minima map to self.
  Follow steepest-descent gradient until reach a site that is a minimum. Filter basins so have to be a min within specified region (disk of radius dRegion).
  Return map of cell to basin.
  
  Arguments:
  vals - values
  cellIsMin - array >0 if cell is local minimum
  cell0 - cell for iterating through mesh
  mesh - Mesh instance
  """
  
  '''
  #to adapt global watershed to region, make values outside of region huge so 
  don't steepest descend that way. since we pass in minCells, do this before call:
  bigVal = 1.e10
  vals = np.copy(valsIn)
  vals[inRegion<1] = bigVal
  '''
  
  cell2Site = -np.ones(mesh.nCells,dtype=int) #so no cell2Site[iCell]=iCell
  
  #get local steepest path
  dMin = min(1.e-6,mesh.r/mesh.nCells);
  for cell in iter(cell0.copy()):
    if (not cell.isInRegion()):
      continue
      
    iCell = cell.ind
    if (cellIsMin[iCell]>0): #steepest path is to self
      cell2Site[iCell]= iCell
    else:
      nbrs = cell.get_nbrInds()
      nbrs = nbrs[mesh.isIndsInRegion(nbrs)]
      
      val0 = vals[iCell]
      valNbrs = vals[nbrs]
      
      #correspondence is towards minimum gradient.
      lat0, lon0 = mesh.get_latLon_inds(iCell)
      latNbrs, lonNbrs = mesh.get_latLon_inds(nbrs)
      
      dNbrs = helpers.calc_distSphere_multiple(mesh.r, lat0, lon0, latNbrs, lonNbrs)
      dNbrs[dNbrs<dMin]=dMin #avoid divide by 0
      #print valNbrs, dNbrs, val0
      valNbrs = (valNbrs-val0)/dNbrs
      iNbr = np.argmin(valNbrs)
      if (True):
        if (valNbrs[iNbr]>=0):
          print "Uhoh. Steepest descent for cell {0}->{2} is {1}, which isn't negative!".format(iCell, valNbrs[iNbr], nbrs[iNbr])
      cell2Site[iCell] = nbrs[iNbr]
  
  #Filter local extrema by area to limit high (spatial) frequency "noise".
  #An extremum must be an extremum within the filter region
  nRedirect = 0
  for cell in iter(cell0.copy()):
    if (not cell.isInRegion()):
      continue
       
    iCell = cell.ind
    if (cellIsMin[iCell]>0):
      #see if cell is min in region, not just neighbors.
      #if not regional min, update cell2Site so local min goes to another basin
      cellsRegion = cell.get_regionInds()
      valsRegion = vals[cellsRegion]
      minInd = np.argmin(valsRegion)
      minVal = valsRegion[minInd]; minCell = cellsRegion[minInd];
      val0 = vals[iCell]; #print val0, minVal
      if (minVal < val0):
        #print "Redirecting cell {0} to {1}".format(iCell, minCell)
        cellIsMin[iCell] = 0
        cell2Site[iCell] = minCell
        nRedirect = nRedirect+1
      else:
        #cell is min, but not necessarily distinct min (ie strictly less than all other values w/in disk).
        #here, we deal with the case where multiple cells in region all have the exact same value.
        #3 (or nLon) neigboring mins should be 1 tpv, not 3 (physically).
        #we'll redirect to the maximum index within disk (so all cells redirect to accepted min).
         
        #While unlikely for 2 general floats to be equal, this can arise from:
        #-idealized initialization
        #-compressed storage of variables (eg, ERA-I stores as shorts where val=short*scale+offset)
        isDiskMin = valsRegion==minVal
        if (np.sum(isDiskMin)>1):
          indsOfMins = cellsRegion[isDiskMin>0]
          minCell = np.max(indsOfMins)
          if (iCell != minCell):
            cellIsMin[iCell] = 0
            cell2Site[iCell] = minCell
            nRedirect = nRedirect+1
  print "Number of min after redirect: ", np.sum(cellIsMin>0)
  
  #follow local steepest path (and any redirections from, say, regional thresholds) to site
  for cell in iter(cell0.copy()):
    if (not cell.isInRegion()):
      continue
      
    iCell = cell.ind
    nextCell = cell2Site[iCell]
    nCount = 0
    while (not cellIsMin[nextCell]>0):
      nextCell = cell2Site[nextCell]
      #print "Cell {0} going to {1}".format(iCell,nextCell); print vals[iCell], vals[nextCell]
      nCount=nCount+1
      if (nCount>mesh.nCells): #something is probably quite wrong
        #seems to happen if values w/in a region have the exact same value. Storing values as shorts in files makes this more likely.
        print "Uhoh, stuck in while loop for cell {0} with value {1}".format(iCell, vals[iCell])
        nbrs = cell.get_nbrInds(); valNbrs = vals[nbrs]
        print "Neighbor's values are: ", valNbrs
        break
    #end not cellIsMin
    cell2Site[iCell] = nextCell

  return (cell2Site, cellIsMin)

def segment_high_low_watershed_region(theta, vort, cell0, mesh, segRestrictPerc):
  """
  Segment a continuous surface into high and low watershed basins
  Steps: get high and low basin seeds, associate cells to both high and low basins if not extrema.
  to decide whether "really" part of high or low basin, we have options:
  -(anti-)cyclonic for (high) low...is local vorticity noisy?
  -closer theta value to maxima a la color scale grouping...huge min or max value now matters
  -whether steeper gradient is to high or low
  -physical distance
  -concavity of surface a la last closed contour
  
  Arguments:
  theta - potential temperature
  vort - vertical vorticity
  cell0 - cell for iterating through mesh
  mesh - Mesh instance
  segRestrictPerc - Percentile for restricting watershed basins (percentile of amplitudes of cells on boundary of watershed basin)
  """
  
  #mins
  print "Finding minima"
  #to adapt global watershed to region, make values outside of region huge so don't steepest descend that way
  bigVal = 1.e10
  vals = np.copy(theta) #so don't affect variable passed in
  vals[np.logical_not( mesh.get_inRegion1d() )] = bigVal
  
  cellIsMin = find_minCells_region_flat(vals, cell0.copy(), mesh)
  cell2SiteMin, cellIsMin = watershed_region(vals, cellIsMin, cell0.copy(), mesh)
  
  #maxs: perform min on an inverted surface
  print "Finding maxima"
  #adapt global watershed to region
  vals = -np.copy(theta)
  vals[np.logical_not( mesh.get_inRegion1d() )] = bigVal
  
  cellIsMax = find_minCells_region_flat(vals, cell0.copy(), mesh)
  cell2SiteMax, cellIsMax = watershed_region(vals, cellIsMax, cell0.copy(), mesh)
  
  #"voting" procedure for low/high classification ------
  print "Associating to max or min by local vorticity"
  cell2Site = -np.ones(mesh.nCells, dtype=int)
  
  for cell in iter(cell0.copy()):
    if (not cell.isInRegion()):
      continue
      
    iCell = cell.ind
    if (cellIsMin[iCell]>0 or cellIsMax[iCell]>0): #allows for cyclonic max. is that right?
      cell2Site[iCell] = iCell
    else:
      #cyclonic ~ sign(vorticity) depends on hemisphere
      signHem = 1 #sign function is problematic since sign(0)=0
      lat0, lon0 = mesh.get_latLon_inds(iCell)      
      if (lat0<0): #lat=0 gets put in NH
        signHem = -signHem
      
      if (signHem*vort[iCell]<0): #anticyclonic
        cell2Site[iCell] = cell2SiteMax[iCell]
      else: #0 or cyclonic
        cell2Site[iCell] = cell2SiteMin[iCell]
  
  #restrict basins to closed contours -------
  print "Restricting basins to last closed contour w/in watershed"
  
  isBoundary = find_basinBoundaries(cell2Site, cell0.copy(), mesh)
  for cell in iter(cell0.copy()):
    iCell = cell.ind
    if (cell2Site[iCell] == iCell): #loop by basin
      theta0 = theta[iCell]
      boundingBasin = (cell2Site==iCell)*(isBoundary)
      thetaBoundary = theta[boundingBasin>0]
      
      #defining the last closed contour as smallest amplitude on boundary covers min and max.
      #cells outside of contour are set to -1 == background
      #minAmp = np.min( np.absolute(thetaBoundary-theta0) )
      minAmp = np.percentile(np.absolute(thetaBoundary-theta0), segRestrictPerc)
      
      print 'theta0, thetaMinBound, thetaMaxBound, minAmp, site: ', theta0, np.min(thetaBoundary), np.max(thetaBoundary), minAmp, iCell
      #[cell2Site[i]=-1 for i in xrange(mesh.nCells) if (cell2Site[i]==iCell and abs(theta[i]-theta0)>minAmp)] #i don't think we can set values in list comprehension
      '''
      for i in xrange(mesh.nCells):
        if (cell2Site[i]==iCell and abs(theta[i]-theta0)>minAmp):
          cell2Site[i]=-1
      '''
      toRemove = (cell2Site==iCell)*(np.absolute(theta-theta0)>minAmp); print 'Removed cells {0}/{1}'.format(np.sum(toRemove), np.sum(cell2Site==iCell))
      cell2Site[toRemove>0] = -1
          
  return (cell2Site, cellIsMin, cellIsMax)

def segment(theta, vort, cell0, mesh, segRestrictPerc):
  """Wrapper for segmenting surface into high and low watershed basins"""
  cell2Site, cellIsMin, cellIsMax = segment_high_low_watershed_region(theta, vort, cell0, mesh, segRestrictPerc)
  sitesMin = cell2Site[cellIsMin>0]
  sitesMax = cell2Site[cellIsMax>0]
  
  return (cell2Site, sitesMin, sitesMax)

def write_netcdf_header_seg(fName, info, nCells, nSitesMax):
  """Make and write header for segmentation tpvTrack netcdf file"""
  #I don't know how to make ragged arrays, so we'll use 
  #array[nTimes,nMax] and nElts[nTimes]
  
  data = netCDF4.Dataset(fName, 'w', format='NETCDF4')
  data.description = info
  
  # dimensions
  data.createDimension('time', None)
  data.createDimension('nCells', nCells)
  data.createDimension('nMax', nSitesMax)

  # variables
  cell2Site_data = data.createVariable('cell2Site', 'i4', ('time','nCells',))
  sitesMin_data = data.createVariable('sitesMin', 'i4', ('time','nMax',)) #could make unsigned...careful about any arithmetic later though
  nSitesMin_data = data.createVariable('nSitesMin', 'i4', ('time',))
  sitesMax_data = data.createVariable('sitesMax', 'i4', ('time','nMax',))
  nSitesMax_data = data.createVariable('nSitesMax', 'i4', ('time',))
  
  #units and descriptions
  cell2Site_data.description = 'Map[cell]->basin'
  sitesMin_data.description = 'Cell indices of minima'
  sitesMax_data.description = 'Cell indices of maxima'
  nSitesMin_data.description = '# minima = # cyclonic tpvs'
  
  return data
  
def write_netcdf_iTime_seg(data, iTime, cell2Site, sitesMin, sitesMax, nSitesMax):
  """Write 1 time into segmentation tpvTrack netcdf file"""
  # fill file. with time as unlimited, dimension will just keep growing
  
  data.variables['cell2Site'][iTime,:] = cell2Site[:]
  
  nSites = len(sitesMin)
  if (nSites>nSitesMax):
    print "Uhoh. Only storing {0}/{1} sites".format(nSitesMax,nSites)
    nSites = nSitesMax
  data.variables['sitesMin'][iTime,0:nSites] = sitesMin[0:nSites]
  data.variables['nSitesMin'][iTime] = nSites
  
  nSites = len(sitesMax)
  if (nSites>nSitesMax):
    print "Uhoh. Only storing {0}/{1} sites".format(nSitesMax,nSites)
    nSites = nSitesMax
  data.variables['sitesMax'][iTime,0:nSites] = sitesMax[0:nSites]
  data.variables['nSitesMax'][iTime] = nSites
#

def run_segment(fSeg, info, dataMetr, cell0, mesh, nTimes, segRestrictPerc=5.):
  """Run segmentation over multiple times of input data"""
  nSitesMax = (mesh.get_inRegion1d()).sum() #can't have more sites than cells...
  nSitesMax = 3000
  dataSeg = write_netcdf_header_seg(fSeg, info, mesh.nCells, nSitesMax)
  
  for iTime in xrange(nTimes):
    print "Segmenting time index: ", iTime
    
    theta = dataMetr.variables['theta'][iTime,:]
    vort = dataMetr.variables['vort'][iTime,:]
    cell2Site, sitesMin, sitesMax = segment(theta, vort, cell0.copy(), mesh, segRestrictPerc)
    
    write_netcdf_iTime_seg(dataSeg, iTime, cell2Site, sitesMin, sitesMax, nSitesMax)
    
  dataSeg.close()

def plot_basins_save(fNameSave, lat, lon, vals, sitesMin, sitesMax):
  """
  Example of plotting segmentation with basins colored by site.
  Input all as 1d arrays. lat/lon in radians
  """
  
  plt.figure()

  #m = Basemap(projection='ortho',lon_0=100,lat_0=60, resolution='l')
  r2d = 180./np.pi
  m = Basemap(projection='ortho',lon_0=0,lat_0=89.95, resolution='l')
  x,y = m(lon*r2d, lat*r2d)
  #print x.shape, y.shape

  m.drawcoastlines(linewidth=.5)
  #m.drawmapboundary()
  
  #plot nan's with different color
  maskedVals = np.ma.array(vals, mask=np.isnan(vals))
  cmap = matplotlib.cm.RdBu_r #matplotlib.cm.jet
  cmap.set_bad('w',1.)
  pPlot = m.pcolor(x,y,maskedVals,tri=True, shading='flat',edgecolors='none',cmap=cmap, vmin=280, vmax=360)
  
  xMin = x[sitesMin]; yMin = y[sitesMin]; m.scatter(xMin, yMin, c='k', marker="v")
  xMax = x[sitesMax]; yMax = y[sitesMax]; m.scatter(xMax, yMax, c='w', marker="^")

  plt.colorbar(pPlot)
  plt.savefig(fNameSave, bbox_inches='tight'); plt.close()

def run_plotBasins(fDirSave, dataMetr, fSeg, mesh):
  """Plot segmentation"""
  lat, lon = mesh.get_latLon_inds(np.arange(mesh.nCells))
  if (True):
    #latLon cells will have duplicate points, which mucks up the triangulation
    #this is a quick hack but a more robust option may be needed...
    sizeDelta = .05*np.pi/180.
    dll = sizeDelta*np.random.uniform(-1.,1,mesh.nCells)
    lat += dll
    lon += dll
  
  dataSeg = netCDF4.Dataset(fSeg,'r')
  info = dataSeg.description
  nTimes = len(dataSeg.dimensions['time'])
  for iTime in xrange(nTimes):
    fName = 'seg_{0}_{1}.png'.format(iTime, info)
    fSave = fDirSave+fName
    
    cell2Site = dataSeg.variables['cell2Site'][iTime,:]
    sitesMin = dataSeg.variables['sitesMin'][iTime,:];
    nMin = dataSeg.variables['nSitesMin'][iTime]; sitesMin = sitesMin[0:nMin]
    sitesMax = dataSeg.variables['sitesMax'][iTime,:]; 
    nMax = dataSeg.variables['nSitesMax'][iTime]; sitesMax = sitesMax[0:nMax]
    
    theta = dataMetr.variables['theta'][iTime,:]
    vals = theta[cell2Site]; vals[cell2Site<0] = np.nan; #print vals
    
    print "Saving file: "+fSave
    plot_basins_save(fSave, lat, lon, vals, sitesMin, sitesMax)
  
  dataSeg.close()

