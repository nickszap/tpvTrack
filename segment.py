import numpy as np
import netCDF4
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import helpers

#watershed has a few options for implementation:
#-for every cell, walk down steepest gradient to the basin

def find_basinBoundaries(cell2Site, cell0, mesh):
  #a basin boundary is a cell with a cell2Site pointing to another basin
  
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

def get_connectedCells_blob(cell0, isValidVal, mesh):
  connCells = []
  addCells = [cell0.ind]
  while (len(addCells)>0):
    cell0Ind = addCells.pop()
    if (cell0Ind not in connCells):
      connCells.append(cell0Ind)
    
    cell0 = cell0.new(cell0Ind)
    nbrs = cell0.get_nbrInds()
    nbrs = nbrs[mesh.isIndsInRegion(nbrs)*(isValidVal[nbrs]>0)]
    for nbr in nbrs:
      if ((nbr not in addCells) and (nbr not in connCells)):
        addCells.append(nbr)
  return connCells

def connect_blobs(isValidVal, cell0, mesh):
  #ID blobs by associating each blob to maximum cell index w/in blob (-1 is background).
  #np.unique(cell2Site) will give blob inds.
  
  cell2Site = -np.ones(mesh.nCells,dtype=int) #so no cell2Site[iCell]=iCell
  
  for cell in iter(cell0.copy()):
    if (not cell.isInRegion()):
      continue
    iCell = cell.ind
    
    if (isValidVal[iCell]>0 and cell2Site[iCell]<0): #valid and not visited
      connCells = get_connectedCells_blob(cell.copy(), isValidVal, mesh)
      cell2Site[connCells] = iCell
  return cell2Site

def associateBlob_toMax(vals, cell2Site, fillVal=-1.e6):
  
  cell2SiteOut = np.copy(cell2Site)
  cellIsExtr = np.zeros(len(cell2Site),dtype=int)
  
  blobInds = np.unique(cell2Site); blobInds = blobInds[blobInds>-1]
  
  for iBlob in blobInds:
    valsInBlob = np.copy(vals)
    valsInBlob[cell2Site!=iBlob] = fillVal
    
    iCell = np.argmax(valsInBlob)
    cell2SiteOut[cell2Site==iBlob] = iCell
    cellIsExtr[iCell] = 1
  
  return (cell2SiteOut, cellIsExtr)

def segment_vorticityThresh(vort, vortThresh, cell0, mesh, fillVal=0.0):
  #Make connected blobs out of \pm large vorticity
  
  print "Finding cyclones"
  #to adapt global watershed to region, make values outside of region huge so don't steepest descend that way
  vals = np.copy(vort) #so don't affect variable passed in
  vals[np.logical_not( mesh.get_inRegion1d() )] = fillVal
  #mask by terrain height
  
  cell2Site = connect_blobs(vort>vortThresh, cell0, mesh)
  cell2SiteCyclone, isCycloneExtr = associateBlob_toMax(vort, cell2Site)
  
  print "Finding anti-cyclones"
  cell2Site = connect_blobs(vort<-vortThresh, cell0, mesh)
  cell2SiteAntiCyclone, isAntiCycloneExtr = associateBlob_toMax(-vort, cell2Site)
  
  cell2Site = np.maximum(cell2SiteCyclone, cell2SiteAntiCyclone)
          
  return (cell2Site, isCycloneExtr, isAntiCycloneExtr)

def segment(vort, vortThresh, cell0, mesh):
  #cell2Site, cellIsMin, cellIsMax = segment_high_low_watershed_region(theta, vort, cell0, mesh)
  cell2Site, cellIsMin, cellIsMax = segment_vorticityThresh(vort, vortThresh, cell0, mesh)
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

def run_segment(fSeg, info, dataMetr, cell0, mesh, nTimes, vortThresh=1.e-5):
  
  nSitesMax = (mesh.get_inRegion1d()).sum() #can't have more sites than cells...
  nSitesMax = 3000
  dataSeg = write_netcdf_header_seg(fSeg, info, mesh.nCells, nSitesMax)
  
  for iTime in xrange(nTimes):
    print "Segmenting time index: ", iTime
    
    #theta = dataMetr.variables['theta'][iTime,:]
    vort = dataMetr.variables['vort'][iTime,:]
    cell2Site, sitesMin, sitesMax = segment(vort, vortThresh, cell0.copy(), mesh)
    
    write_netcdf_iTime_seg(dataSeg, iTime, cell2Site, sitesMin, sitesMax, nSitesMax)
    
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
  
  #plot nan's with different color
  maskedVals = np.ma.array(vals, mask=np.isnan(vals))
  cmap = matplotlib.cm.jet
  cmap.set_bad('w',1.)
  pPlot = m.pcolor(x,y,maskedVals,tri=True, shading='flat',edgecolors='none',cmap=cmap, vmin=-1.e-4, vmax=1.e-4)
  
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

