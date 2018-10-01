import numpy as np

import helpers

def findOwner_horizNbrs_latLon(pt_ll, cellId, latCell, lonCell, nEdgesOnCell, cellsOnCell):
  """
  Given a guess cell, walk towards input latLon until no cell center is closer.
  Then, the input point is w/in voronoi region ("owner") of that cell.
  Our spherical domain is convex so we can walk from one point to any other.
  
  Arguments:
  pt_ll - [latitude,longitude] of input point (in radians)
  cellId - initial guess cell
  latCell - latitudes of mesh
  lonCell - longitudes of mesh
  nEdgesOfCell - number of edges for each cell in mesh
  cellsOnCell - indices of neighbors for each cell in mesh
  """
  
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
  """
  Using seed cell, gather all points with cell centers within radius of the specified point
  with a flood fill type approach.
  Return a list of cells w/in specified radius.
  
  Arguments:
  pt_ll - [latitude,longitude] of input point
  rDisk - radius of disk about input point
  rEarth - radius of sphere
  c0 - index of cell containing pt_ll
  nEdgesOfCell - number of edges for each cell in mesh
  cellsOnCell - indices of neighbors for each cell in mesh
  latCell - latitudes of mesh
  lonCell - longitudes of mesh
  """

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

class Mesh(object):
  """ Define topology of domain for MPAS dual of triangulation mesh"""
  def __init__(self,lat,lon, areaCell, cellsOnCell, nCellsOnCell, r, rDisk):
    """
    Initialize Mesh instance
    
    Arguments:
    lat - latitudes of mesh
    lon - longitudes of mesh
    areaCell - area of cells
    cellsOnCell - indices of neighbors of each cell
    nCellsOnCell - number of neighbors for each cell
    r - radius of sphere
    rDisk - radius of neighborhood
    """
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
    """
    Return index of closest cell to an input point
    
    Arguments:
    latPt - latitude of point
    lonPt - longitude of point
    guessCell - (optional) initial guess cell
    """
    pt_ll = np.empty(2,dtype=float);
    pt_ll[0] = latPt; pt_ll[1] = lonPt;
    iCell = findOwner_horizNbrs_latLon(pt_ll, guessCell, self.lat, self.lon, self.nCellsOnCell, self.cellsOnCell)
    
    #print 'closestCell2Pt (latPt,lonPt), (latCell,lonCell) ({0},{1}), ({2},{3})'.format(latPt, lonPt, self.lat[iCell], self.lon[iCell])
    return iCell
  
  def fill_inRegion(self, latThresh):
    """Define the subset of domain that is used """
    notInRegion = self.lat<latThresh
    self.inRegion[notInRegion] = 0
  
  def get_inRegion1d(self):
    """Return boolean of if cell is in domain that is used """
    return self.inRegion>0
  
  def isIndsInRegion(self, inds):
    """Slice inRegion with inds"""
    return self.inRegion[inds]>0
    
  def get_latLon_inds(self, inds):
    """Return the latitudes and longitudes of input cells"""
    return (self.lat[inds], self.lon[inds])
    
  def get_area_inds(self, inds):
    """Return the areas of input cells"""
    return self.areaCell[inds]
    
class Cell(object):
  """Define the location of a cell in a mesh"""
  def __init__(self,mesh,ind):
    """
    Initialize instance
    
    Arguments:
    mesh - Mesh instance
    ind - cell index
    """
    self.mesh = mesh
    self.ind = ind
  
  def __iter__(self):
    return self

  def next(self): # Python 3: def __next__(self)
    """Iterate to next cell in mesh"""
    if self.ind < self.mesh.nCells-1:
      self.ind = self.ind + 1
      return self
    else:
      raise StopIteration
  
  def copy(self):
    """Copy cell"""
    return Cell(self.mesh, self.ind)
      
  def isInRegion(self):
    """Boolean of cell in region that is used"""
    return self.mesh.inRegion[self.ind]>0

  def get_nbrInds(self):
    """Indices of neighbors of cell in mesh"""
    nNbrs = self.mesh.nCellsOnCell[self.ind]
    iNbrs = self.mesh.cellsOnCell[self.ind,0:nNbrs]
    return iNbrs
  
  def get_regionInds(self):
    """Indices of cells in disk neighborhood"""
    pt_ll = np.empty(2,dtype=float)
    pt_ll[0] = self.mesh.lat[self.ind]; pt_ll[1] = self.mesh.lon[self.ind]
    iCells = gatherCells_radius(pt_ll, self.mesh.rDisk, self.mesh.r, self.ind, self.mesh.cellsOnCell, self.mesh.nCellsOnCell, self.mesh.lat, self.mesh.lon)
    return np.array(iCells,dtype=int)
  
  def get_areaCell(self):
    """Return area of cell"""
    return self.mesh.areaCell[self.ind]
