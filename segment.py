import numpy as np
import netCDF4

import helpers

#watershed has a few options for implementation:
#-for every cell, walk down steepest gradient to the basin

def find_minCells_region_flat(vals, cell0, mesh):
  #return array[nCells] with 1 if cell is min and cell in region
  
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
  
  cell2Site = -np.ones(mesh.nCells,dtype=int) #so no cell2Site[iCell]=iCell      

  #get local steepest path
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
      dMin = mesh.r/1.e16; dNbrs[dNbrs<dMin]=dMin #avoid divide by 0
      #print valNbrs, dNbrs, val0
      valNbrs = (valNbrs-val0)/dNbrs
      iNbr = np.argmin(valNbrs)
      cell2Site[iCell] = nbrs[iNbr]
  
  #Filter local extrema by area to limit high (spatial) frequency "noise".
  #For multiple close mins, the smallest counts as min for that region.
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
  print "Number of redirects for regional min: ", nRedirect
  
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

def segment_high_low_watershed_region(theta, vort, cell0, mesh):
  #get high and low basin seeds, associate cells to both high and low basins if not extrema.
  #to decide whether "really" part of high or low basin, we have options:
  #-(anti-)cyclonic for (high) low...is local vorticity noisy?
  #-closer theta value to maxima a la color scale grouping...huge min or max value now matters
  #-whether steeper gradient is to high or low
  #-physical distance
  #-concavity of surface a la last closed contour
  
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
          
  return (cell2Site, cellIsMin, cellIsMax)

def segment(theta, vort, cell0, mesh):
  cell2Site, cellIsMin, cellIsMax = segment_high_low_watershed_region(theta, vort, cell0, mesh)
  sitesMin = cell2Site[cellIsMin>0]
  sitesMax = cell2Site[cellIsMax>0]
  
  return (cell2Site, sitesMin, sitesMax)

def write_netcdf_header_seg(fName, info, nCells, nSitesMax):
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
  # fill file. with time as unlimited, dimension will just keep growing
  
  data.variables['cell2Site'][iTime,:] = cell2Site[:]
  
  nSites = len(sitesMin)
  if (nSites>nSitesMax):
    print "Uhoh. Only storing {0}/{1} sites".format(nSitesMax/nSites)
    nSites = sitesMax
  data.variables['sitesMin'][iTime,0:nSites] = sitesMin[0:nSites]
  data.variables['nSitesMin'][iTime] = nSites
  
  nSites = len(sitesMax)
  if (nSites>nSitesMax):
    print "Uhoh. Only storing {0}/{1} sites".format(nSitesMax/nSites)
    nSites = sitesMax
  data.variables['sitesMax'][iTime,0:nSites] = sitesMax[0:nSites]
  data.variables['nSitesMax'][iTime] = nSites
#

def run_segment(fSeg, info, dataMetr, cell0, mesh):
  
  nSitesMax = (mesh.get_inRegion1d()).sum() #can't have more sites than cells...
  dataSeg = write_netcdf_header_seg(fSeg, info, mesh.nCells, nSitesMax)
  
  nTimes = len(dataMetr.dimensions['time'])
  for iTime in xrange(nTimes):
    print "Segmenting time index: ", iTime
    
    theta = dataMetr.variables['theta'][iTime,:]
    vort = dataMetr.variables['vort'][iTime,:]
    cell2Site, sitesMin, sitesMax = segment(theta, vort, cell0.copy(), mesh)
    
    write_netcdf_iTime_seg(dataSeg, iTime, cell2Site, sitesMin, sitesMax, nSitesMax)
    
  dataSeg.close()


