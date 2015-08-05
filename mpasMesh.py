import numpy as np
import netCDF4

import helpers

def findOwner_horizNbrs_latLon(pt_ll, cellId, latCell, lonCell, nEdgesOnCell, cellsOnCell):
  #given a guess cell, walk towards input latLon until no cell center is closer.
  #Then, the input point is w/in voronoi region ("owner") of that cell.
  #Our spherical domain is convex so we can walk from one point to any other.
  
  radEarth = 1.0 #we don't need the "actual earth's" distance, just relative (closer and farther)
  
  dCell = helpers.calc_distSphere_multiple(radEarth, pt_ll[0], pt_ll[1], latCell[cellId], lonCell[cellId])

  flag = 1;
  while (flag==1):
    #keep going towards the cell that's closer to point until no nbr is closer.
    flag =0;
    nNbrs = nEdgesOnCell[cellId]
    nbrs = cellsOnCell[cellId,0:nNbrs]
    dNbrs = helpers.calc_distSphere_multiple(radEarth, pt_ll[0], pt_ll[1], latCell[nbrs], lonCell[nbrs])
    iNbr = np.argmin(dNbrs)
    if (dNbrs[iNbr]<dCell): #then use closest cell for next iteration
      dCell = dNbrs[iNbr]
      cellId = nbrs[iNbr]
      flag=1
  return(cellId);

def gatherCells_radius(pt_ll, rDisk, rEarth, c0, cellsOnCell, nEdgesOnCell, latCell, lonCell):
  #using seed cell, gather all points with cell centers within radius of the specified point
  #with a flood fill type approach.
  #return a list of cells w/in specified radius.

  candidates = [c0]
  closeCells = []; farCells=[]
  while (len(candidates)>0):
    c0 = candidates.pop()
    dCell = helpers.calc_distSphere_multiple(rEarth, pt_ll[0], pt_ll[1], latCell[c0], lonCell[c0])

    if (dCell<rDisk):
      closeCells.append(c0)
      nbrs = cellsOnCell[c0,0:nEdgesOnCell[c0]]
      for iCell in nbrs:
        if ((iCell not in closeCells) and (iCell not in farCells) and (iCell not in candidates)):
          candidates.append(iCell)
    else:
      farCells.append(c0)

  return closeCells

def make_nearestRegrid(fMPAS, latsNew, lonsNew):
  #Input lat/lon arrays (in radians)
  #Return corresponding MPAS cell indices
  
  nLat = len(latsNew); nLon = len(lonsNew);
  ll2MPASCell = np.empty((nLat,nLon),dtype=int)
  
  #mpas connectivity info
  data = netCDF4.Dataset(fMPAS, 'r')
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];
  data.close()
  
  #associate to nearest neighbor
  iCell=0; pt_ll = np.empty(2,dtype=float)
  for iLat in xrange(nLat):
    for iLon in xrange(nLon):
      pt_ll[0] = latsNew[iLat]; pt_ll[1] = lonsNew[iLon]
      iCell = findOwner_horizNbrs_latLon(pt_ll, iCell, latCell, lonCell, nEdgesOnCell, cellsOnCell)
      ll2MPASCell[iLat,iLon] = iCell
  
  if (True):
    #print some diagnostics on how close the regrid cells are
    dError = np.empty((nLat,nLon),dtype=float)
    for iLat in xrange(nLat):
      for iLon in xrange(nLon):
        pt_ll[0] = latsNew[iLat]; pt_ll[1] = lonsNew[iLon]
        dError[iLat,iLon] = helpers.calc_distSphere_multiple(6371., pt_ll[0], pt_ll[1], latCell[ll2MPASCell[iLat,iLon]], lonCell[ll2MPASCell[iLat,iLon]])
    #
    print 'Mean,max,min nearest neighbor regrid distance errors (km): ', np.mean(dError), np.max(dError), np.min(dError)
  
  return ll2MPASCell

class Mesh(object):
  def __init__(self,lat,lon, areaCell, cellsOnCell, nCellsOnCell, r, rDisk):
    self.r = r
    self.lat = lat
    self.lon = lon
    nCells = len(lat)
    self.nCells = nCells
    self.areaCell = areaCell
    self.inRegion = np.ones(nCells,dtype=int)
    self.cellsOnCell = cellsOnCell
    self.nCellsOnCell = nCellsOnCell
    self.rDisk = rDisk
    self.info = 'mpas'
  
  def get_closestCell2Pt(self, latPt, lonPt, guessCell=0):
    #using a seed cell is not currently needed for the latLon caller
    pt_ll = np.empty(2,dtype=float);
    pt_ll[0] = latPt; pt_ll[1] = lonPt;
    iCell = findOwner_horizNbrs_latLon(pt_ll, guessCell, self.lat, self.lon, self.nCellsOnCell, self.cellsOnCell)
    
    #print 'closestCell2Pt (latPt,lonPt), (latCell,lonCell) ({0},{1}), ({2},{3})'.format(latPt, lonPt, self.lat[iCell], self.lon[iCell])
    return iCell
  
  def fill_inRegion(self, latThresh):
    notInRegion = self.lat<latThresh
    self.inRegion[notInRegion] = 0
  
  def get_inRegion1d(self):
    return self.inRegion>0
  
  def isIndsInRegion(self, inds):
    return self.inRegion[inds]>0
    
  def get_latLon_inds(self, inds):
    return (self.lat[inds], self.lon[inds])
    
  def get_area_inds(self, inds):
    return self.areaCell[inds]
    
class Cell(object):
  def __init__(self,mesh,ind):
    self.mesh = mesh
    self.ind = ind
  
  def __iter__(self):
    return self

  def next(self): # Python 3: def __next__(self)
    if self.ind < self.mesh.nCells-1:
      self.ind = self.ind + 1
      return self
    else:
      raise StopIteration
  
  def copy(self):
    return Cell(self.mesh, self.ind)
      
  def isInRegion(self):
    return self.mesh.inRegion[self.ind]>0

  def get_nbrInds(self):
    nNbrs = self.mesh.nCellsOnCell[self.ind]
    iNbrs = self.mesh.cellsOnCell[self.ind,0:nNbrs]
    return iNbrs
  
  def get_regionInds(self):
    pt_ll = np.empty(2,dtype=float)
    pt_ll[0] = self.mesh.lat[self.ind]; pt_ll[1] = self.mesh.lon[self.ind]
    iCells = gatherCells_radius(pt_ll, self.mesh.rDisk, self.mesh.r, self.ind, self.mesh.cellsOnCell, self.mesh.nCellsOnCell, self.mesh.lat, self.mesh.lon)
    return np.array(iCells,dtype=int)
  
  def get_areaCell(self):
    return self.mesh.areaCell[self.ind]
