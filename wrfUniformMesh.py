#This is for a LAM WRF mesh that has uniform grid (projected) dx and dy.

import numpy as np
import helpers

def calc_area(dx, dy, mapfac):
  """
  Calculate cell area
  #input computational dx,dy in meters
  #output area in m^2
  """
  
  area = dx*dy/(mapfac*mapfac)
  return area

#It may be faster to use the map projection to find the closest cell to a point...
#maybe someone else wants to figure out how to do that :) !!!
def findOwner_horizNbrs_latLon(mesh, latPt, lonPt, cellId):
  """
  #given a guess cell, walk towards input latLon until no cell center is closer.
  #Then, the input point is w/in voronoi region ("owner") of that cell.
  #If domain is convex, can walk from one point to any other.
  """
  
  radEarth = 1.0 #we don't need the "actual earth's" distance, just relative (closer and farther)
  latNew, lonNew = mesh.get_latLon_inds(cellId)
  dCell = helpers.calc_distSphere_multiple(radEarth, latPt, lonPt, latNew, lonNew)
  
  flag = 1;
  while (flag==1):
    #keep going towards the cell that's closer to point until no nbr is closer.
    flag =0;
    cell0 = Cell(mesh, cellId)
    nbrs = cell0.get_nbrInds()
    latNew, lonNew = mesh.get_latLon_inds(nbrs)
    
    dNbrs = helpers.calc_distSphere_multiple(radEarth, latPt, lonPt, latNew, lonNew)
    iNbr = np.argmin(dNbrs)
    if (dNbrs[iNbr]<dCell): #then use closest cell for next iteration
      dCell = dNbrs[iNbr]
      cellId = nbrs[iNbr]
      flag=1
  return cellId;

def get_closestSeed2Pt(mesh, latPt, lonPt):
    """return closest of predefined sites """
    
    nSeedsx = 6; nSeedsy = 7;
    seedCandidates_x = np.linspace(0, mesh.nx, num = nSeedsx).astype(int)
    seedCandidates_y = np.linspace(0, mesh.ny, num = nSeedsy).astype(int)
    seedCandidates = helpers.index_2dTo1d(seed_candidates_y, seed_candidates_x, mesh.nx)
    
    d = helpers.calc_distSphere_multiple(1.0, latPt, lonPt, mesh.lat[seedCandidates], mesh.lon[seedCandidates])
    ind = np.argmin(d)
    iCell = seedCandidates[ind]
    if (False):
      r2d = 180./np.pi
      print "Closest cell {0},{1} to coordinate {2},{3}".format(self.lat[iCell]*r2d, self.lon[iCell]*r2d, latPt*r2d, lonPt*r2d)
    return iCell

class Mesh(object):
  """Define topology of domain for a uniformly spaced WRF mesh"""
  def __init__(self,lat,lon, dx, dy, r, rDisk):
    self.r = r
    self.dx = dx
    self.dy = dy
    ny, nx = lat.shape; nCells = nx*ny
    self.nx = nx
    self.ny = ny
    self.nCells = nCells
    self.lat = helpers.flatten_2dTo1d(lat, ny, nx)
    self.lon = helpers.flatten_2dTo1d(lon, ny, nx)
    self.areaCell = np.empty(nCells,dtype=float)
    self.inRegion = np.ones(nCells,dtype=int)
    self.rDisk = rDisk
    self.info = 'wrf'
    
  def fill_cellArea(self, mapFac2d):
    mapFac1d = helpers.flatten_2dTo1d(mapFac2d, self.ny, self.nx) 
    #The above takes scalar to 1-element array since it's just calling np.ravel(). If flatten_2dTo1d changes, make sure that still works!!!
    areaCell = calc_area(self.dx, self.dy, mapFac1d)
    self.areaCell[:] = areaCell #works for both cases of RHS is single value or array with same size
  
  def isPointInDomain(self, latPt, lonPt, iCell):
    #given a point and closest cell, return whether that point is actually within the domain.
    #for a temporary solution, we'll say that a point who's closest cell is on the border of the
    #domain is actually outside the domain.
    iy, ix = helpers.index_1dTo2d(iCell, self.nx)
    
    inDomain = True
    if (iy==0 or iy==self.ny-1):
      inDomain = False
    elif (ix==0 or ix==self.nx-1):
      inDomain = False
    #fwiw, False==0 and True==1 are standard in Python 2.6 and 3. It's likely entrenched enough to work with as such 
    
    return inDomain
  
  def get_closestCell2Pt(self, latPt, lonPt, guessCell=0):
    #distance from all points is a brute force solution.
    #there are some more sophisticated things like walking nearest nbrs from, say, the closest of 10 predefined sites.
    
    iCell = findOwner_horizNbrs_latLon(self, latPt, lonPt, guessCell)
    if (False):
      r2d = 180./np.pi
      print "Closest cell {0},{1} to coordinate {2},{3}".format(self.lat[iCell]*r2d, self.lon[iCell]*r2d, latPt*r2d, lonPt*r2d)
    return iCell
  
  def fill_inRegion(self, latThresh):
    self.inRegion[self.lat<latThresh] = 0
  
  def get_inRegion1d(self):
    return self.inRegion[:]>0
  
  def isIndsInRegion(self, inds):
    return self.inRegion[inds]>0
    
  def get_latLon_inds(self, inds):
    return (self.lat[inds], self.lon[inds])
    
  def get_area_inds(self, inds):
    return self.areaCell[inds]
    
class Cell(object):
  def __init__(self,mesh,ind):
    iy, ix = helpers.index_1dTo2d(ind, mesh.nx)
    self.ix = ix
    self.iy = iy
    self.mesh = mesh
    self.ind = ind
  
  def __iter__(self):
    return self

  def next(self): # Python 3: def __next__(self)
    if self.ind < self.mesh.nCells-1:
      self.ind = self.ind + 1
      self.iy,self.ix = helpers.index_1dTo2d(self.ind, self.mesh.nx)
      return self
    else:
      raise StopIteration
  
  def copy(self):
    return Cell(self.mesh, self.ind)
      
  def isInRegion(self):
    return self.mesh.inRegion[self.ind]>0
      
  def nbrInds_2d(self):
    #return list of 8-conn nbr indices.
    #don't have neighbors to the sides outside the domain
    
    ix = self.ix; iy = self.iy
    nx = self.mesh.nx; ny = self.mesh.ny
    
    iLeft =ix-1; iRight = ix+1;
    iDown =iy-1; iUp = iy+1;
    
    haveL = iLeft>-1; haveR = iRight<nx
    haveD = iDown>-1; haveU = iUp<ny
    nbrx = []; nbry = []
    #top and bottom strips
    if (haveU):
      nbrx.append(ix); nbry.append(iUp)
      if (haveL):
        nbrx.append(iLeft); nbry.append(iUp)
      if (haveR):
        nbrx.append(iRight); nbry.append(iUp)
    if (haveD):
      nbrx.append(ix); nbry.append(iDown)
      if (haveL):
        nbrx.append(iLeft); nbry.append(iDown)
      if (haveR):
        nbrx.append(iRight); nbry.append(iDown)
    #just left and right left
    if (haveL):
      nbrx.append(iLeft); nbry.append(iy)
    if (haveR):
      nbrx.append(iRight); nbry.append(iy)
    
    if (False):
      print "Cell (i{0},ix{1},iy{2}) has nbrs: {3}, {4}".format(self.ind, self.ix, self.iy, nbrx, nbry)
      
    return (nbrx, nbry)

  def get_nbrInds(self):
    nbrx, nbry = self.nbrInds_2d()
    nbrx = np.array(nbrx,dtype=int); nbry=np.array(nbry, dtype=int)
    return helpers.index_2dTo1d(nbry,nbrx,self.mesh.nx)
  
  def get_areaCell(self):
    return self.mesh.areaCell[self.ind]

