#neighborhood operations for latLon,mpas,wrf,... meshes

import numpy as np

#print "\nFor latLon: lats[0->1] goes South. lons[0->1] goes east.\n" #it doesn't "really" matter since we only use neighborhoods

def index_2dTo1d(iLat, iLon, nLon):
  return iLat*nLon+iLon
  
def index_1dTo2d(ind, nLon):
  iLon = ind%nLon
  iLat = (ind-iLon)/nLon
  return (iLat, iLon)

def flatten_2dTo1d(vals, nLat, nLon):
  #vals[lat][lon] goes to vals[iLat*nLon+iLon]
  
  valsOut = np.ravel(vals)
  return valsOut

def unflatten_1dTo2d(vals, nLat, nLon):
  #vals[iLat*nLon+iLon] goes to vals[lat][lon]
  
  valsOut = np.reshape(vals, (nLat,nLon))
  return valsOut

def area_latLonCell(latN, latS, dLon, r):
  #input angles in radians of northern bound, southern bound, and width of rectangle in radians
  #r is radius of circle in meters
  
  solidAngle = dLon*(np.sin(latN)-np.sin(latS))
  area = solidAngle*r*r
  return area
    
def calc_areaLatStrips(lat, r):
  #return the areas of the latitude strips centered around the specified points.
  #lat[0] is north and lat[nLat-1] is south.
  #for multiple evenly spaced cells at a given latitude, areaCell = areaLats[iLat]/nLon
  
  #caps of sphere for poles and latitude strips for rest.
  dLon = 2.*np.pi
  nLat = len(lat)
  areaLats = np.empty(nLat, dtype=float)
  nonPole = np.arange(1,nLat-1) #[1,nlat-2]
  latsN = .5*(lat[nonPole]+lat[nonPole-1]) #face is halfway between centers
  latsS = .5*(lat[nonPole]+lat[nonPole+1])
  areaLats[nonPole] = area_latLonCell(latsN, latsS, dLon, r)
  
  #assuming symmetric latitudes for poles, area of cap is just residual.
  #this probably has more roundoff error than just computing area of cap zone, but
  #this doesn't matter too much for weighting since we normalize by sum(weights) anyway.
  nonPoleArea = np.sum(areaLats[nonPole]); totArea = 4.*np.pi*r*r
  #print nonPoleArea, totArea
  areaCap = (totArea-nonPoleArea)/2
  areaLats[0] = areaCap; areaLats[nLat-1] = areaCap
  
  return areaLats

class Mesh(object):
  def __init__(self,lat,lon, r):
    self.r = r
    self.lat = lat
    self.lon = lon
    self.nLat = len(lat)
    self.nLon = len(lon)
    self.areaCell = self.fill_latCellArea()
    self.inRegion = np.ones((nLat,nLon),dtype=int)
  
  def fill_latCellArea(self):
    areaLat = calc_areaLatStrips(self.lat, self.r)
    self.areaCell = areaLat/self.nLon

class Cell(object):
  def __init__(self,mesh,iLat,iLon):
    self.iLat = 0
    self.iLon = 0
    self.mesh = mesh
  
  def isInRegion(self):
    return self.mesh.inRegion[self.iLat,self.iLon]>0
      
  def nbrInds_ll(self):
    #return list of 8-conn nbr indices: [(latNbr1, lonNbr1), (latNbr2, lonNbr2),...]
    #always have east, west neighbors. not north/south at respective pole
    
    iLat = self.iLat; iLon = self.iLon
    nLat = self.mesh.nLat; nLon = self.mesh.nLon
    
    iWest = (iLon-1)%nLon # -1%4=3 so don't worry about negatives
    iEast = (iLon+1)%nLon
    iSouth = iLat+1
    iNorth = iLat-1
    
    haveSouth = iSouth<nLat
    haveNorth = iNorth>-1
    
    nbrLats = [iLat, iLat]; nbrLons = [iWest, iEast]
    if (haveSouth):
      nbrLats.extend([iSouth, iSouth, iSouth])
      nbrLons.extend([iWest, iLon, iEast])
      
    if (haveNorth):
      nbrLats.extend([iNorth, iNorth, iNorth])
      nbrLons.extend([iWest, iLon, iEast])
      
    return (nbrLats, nbrLons)

  def nbrInds_ll_flat(self):
    nbrInds_lat, nbrInds_lon = self.nbrInds_ll()
    nbrInds_lat = np.array(nbrInds_lat); nbrInds_lon = np.array(nbrInds_lon)
    return index_2dTo1d(self.iLat,self.iLon,self.mesh.nLon)
