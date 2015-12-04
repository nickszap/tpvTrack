import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os

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
  #Keep the regional min as the min for that area.
  
  #The edge case is when there are >1 regional extremum, ie not distinct min (duplicate values).
  #3 (or nLon) neigboring mins should be 1 tpv, not 3 (physically).
  #we'll redirect to the maximum index within disk (so all cells redirect to accepted min).
  #While unlikely for 2 general floats to be equal, this can arise from:
  #-idealized initialization
  #-compressed storage of variables (eg, ERA-I stores as shorts where val=short*scale+offset)
  
  #The logic is simplified if we just compare each min to all others (rather than looking at cells in region about each extremum). 
  #Plus side is that now we don't even need a function for cell.get_regionInds()
  localMins = np.arange(mesh.nCells,dtype=int)[cellIsMin>0]
  nMins = len(localMins)
  distMatrix = np.empty((nMins,nMins),dtype=float)
  for iMin in xrange(nMins):
    #just calculate distances for upper right triangle of nMins x nMins matrix, since distances are symmetric)
    iCell = localMins[iMin]
    #can't skip if already found not to be min in disk when checking some previous min if want iteration-order independence
    
    #distance to other mins
    nbrs = localMins[iMin+1:]
    lat0, lon0 = mesh.get_latLon_inds(iCell)
    latNbrs, lonNbrs = mesh.get_latLon_inds(nbrs)
    dNbrs = helpers.calc_distSphere_multiple(mesh.r, lat0, lon0, latNbrs, lonNbrs)
    distMatrix[iMin,iMin] = 0.;
    distMatrix[iMin,iMin+1:] = dNbrs[:]
  for iMin in xrange(nMins):
    #make d symmetric
    distMatrix[iMin+1:,iMin] = distMatrix[iMin,iMin+1:]
  
  for iMin in xrange(nMins):
    #take regional min for disk
    iCell = localMins[iMin]
    dNbrs = distMatrix[iMin,:]
    
    cellsRegion = localMins[dNbrs<mesh.rDisk]
    valsRegion = vals[cellsRegion]
    minVal = np.min(valsRegion)
    #for >1 mins, keep the one with largest index.
    #the particular condition (eg, keep largest index) doesn't matter, it just needs to be globally valid.
    dupMins = cellsRegion[valsRegion==minVal]
    minCell = np.max(dupMins)
    #update if min isn't kept as regional min
    if (iCell != minCell):
      cellIsMin[iCell] = 0
      cell2Site[iCell] = minCell
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

def run_segment(fSeg, info, dataMetr, cell0, mesh, iTimeStart, iTimeStop):
  
  #nSitesMax = (mesh.get_inRegion1d()).sum() #can't have more sites than cells...
  nSitesMax = 3000
  dataSeg = write_netcdf_header_seg(fSeg, info, mesh.nCells, nSitesMax)
  
  #for iTime in xrange(nTimes):
  iTimeFile = 0
  for iTime in xrange(iTimeStart,iTimeStop+1):
    print "Segmenting time index: ", iTime
    
    theta = dataMetr.variables['theta'][iTime,:]
    vort = dataMetr.variables['vort'][iTime,:]
    cell2Site, sitesMin, sitesMax = segment(theta, vort, cell0.copy(), mesh)
    
    write_netcdf_iTime_seg(dataSeg, iTimeFile, cell2Site, sitesMin, sitesMax, nSitesMax)
    iTimeFile = iTimeFile+1 #so you don't have 10k times in a file w/ only 100 valid times...
    
  dataSeg.close()

def plot_basins_save(fNameSave, lat, lon, vals, sitesMin, sitesMax):
  #Input all as 1d arrays. lat/lon in radians
  
  plt.figure()

  #m = Basemap(projection='ortho',lon_0=100,lat_0=60, resolution='l')
  r2d = 180./np.pi
  m = Basemap(projection='ortho',lon_0=0,lat_0=89.95, resolution='l')
  x,y = m(lon*r2d, lat*r2d)
  #print x.shape, y.shape

  m.drawcoastlines(linewidth=.5)
  #m.drawmapboundary()
  
  pPlot = m.pcolor(x,y,vals,tri=True, shading='flat',edgecolors='none',cmap=plt.cm.jet, vmin=280, vmax=360)
  
  xMin = x[sitesMin]; yMin = y[sitesMin]; m.scatter(xMin, yMin, c='k', marker="v")
  xMax = x[sitesMax]; yMax = y[sitesMax]; m.scatter(xMax, yMax, c='w', marker="^")

  plt.colorbar(pPlot)
  plt.savefig(fNameSave, bbox_inches='tight'); plt.close()

def run_plotBasins(fDirSave, dataMetr, fSeg, mesh):
  
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
  #for iTime in range(242,246+1)+range(364,368+1)+range(485,nTimes):
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

def combineParallelFiles(fOut, iTimesStart, iTimesEnd, filesIn, keys1d = ['nSitesMin','nSitesMax'], keys2d=['cell2Site','sitesMin','sitesMax']):
  #stitch the contiguous chunks of times from the separate workers together into 1 file
  
  #rather than remake the header of a netcdf file, append to a copy of one of the existing
  cmd = 'cp {0} {1}'.format(filesIn[0], fOut)
  print cmd; os.system(cmd)
  
  dataOut = netCDF4.Dataset(fOut,'a')
  nFilesIn = len(filesIn)
  #I'm not sure if it's better to do the following by reading xor writing contiguously...???
  for iFile in xrange(nFilesIn):
    fIn = filesIn[iFile]; iStart = iTimesStart[iFile]; iEnd = iTimesEnd[iFile]
    dataIn = netCDF4.Dataset(fIn,'r')
    for key in keys1d:
      vals = dataIn.variables[key][:]
      dataOut.variables[key][iStart:iEnd+1] = vals[:]
    for key in keys2d:
      vals = dataIn.variables[key][:,:]
      dataOut.variables[key][iStart:iEnd+1,:] = vals[:,:]
    
    dataIn.close()
  dataOut.close()
  


