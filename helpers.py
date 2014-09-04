import numpy as np

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