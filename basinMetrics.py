import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
#from scipy import stats

#import sys; sys.path.append("/home/nickszap/Dropbox/pythonScripts/minAreaBoundingRectangle/")
#import bbox_test
import helpers

'''
#We can characterize the shapes of TPVs in various ways.
Given discrete basins in cells, we have >= info needed.
#For measures of asymmetry:
-eccentricity (fit ellipse, "spindles" to every boundary point,...)
-skewness about extremum
-moment of inertia about extremum
-distance of extremum to boundary 
  (can do spherical distance to boundary or 
  time=distance/speed with sign as vel.distance)

The boundary is cell inside w/ neighbor outside.
-arc length of boundary
-great circle distance vs distance w/in object

"Volume" measures:
-1/A \int(theta-theta0 dA)
'''

#circulation, amplitude, area, moment of inertia
#solidity: area(basin)/area(convexHull) but apparently convex hull of region on a sphere
#isn't straightforward...could do on map projection
r2d = 180./np.pi
  
def calc_circulation(basinInds, vort, mesh):
  #\int gradxu . dA
  
  vortCells = vort[basinInds]
  cellAreas = mesh.get_area_inds(basinInds)
  return np.dot(vortCells, cellAreas)

def calc_secondMomentArea(site0, basinInds, mesh):  
  #moment of inertia about extremum: \int r^2 dA
  #if want rotational moment, need some form of mass for sum(m r^2)
  #return in units of km^4
  
  #moment about extremum
  lat0, lon0 = mesh.get_latLon_inds(site0)
  latBasin, lonBasin = mesh.get_latLon_inds(basinInds)
  
  d = helpers.calc_distSphere_multiple(mesh.r/1.e3, lat0, lon0, latBasin, lonBasin)
  #d /= 1.e3
  cellAreas = mesh.get_area_inds(basinInds)/1.e6
  J = np.dot(d*d, cellAreas)
  return J
  
def calc_amplitude_maxMin(basinInds, thetaFlat):
  #max-min value
  
  valsBasin = thetaFlat[basinInds]
  minVal = np.min(valsBasin); maxVal = np.max(valsBasin)
  return maxVal-minVal

def calc_fieldVolume(basinInds, thetaFlat, valRef, mesh):
  #"volume": \int (theta-thetaRef) dA
  #return in units of K*km^2
  
  #moment about extremum
  vals = thetaFlat[basinInds]
  #valRef = np.max(vals)
  vals = vals-valRef
  cellAreas = mesh.get_area_inds(basinInds)/1.e6
  vol = np.dot(vals, cellAreas)
  return vol
  
def calc_area(basinInds, mesh):
  
  cellAreas = mesh.get_area_inds(basinInds)
  return np.sum(cellAreas)

metricKeys = 'circ vortMean ampMaxMin rEquiv thetaVol ampMean thetaExtr latExtr lonExtr'.split()
def calc_metrics(sites, cell2Site, vort, theta, mesh):
  #calculate input basins' properties and return metrics as dictionary
  
  #initialize
  nSites = len(sites)
  metrics = {}
  for key in metricKeys:
    metrics[key] = np.empty(nSites, dtype=float)
  
  #calc fields
  for iSite in xrange(nSites):
    site0 = sites[iSite]
    inBasin = cell2Site==site0
    basinInds = np.arange(mesh.nCells,dtype=int)[inBasin]
    
    #extreme values --------------------
    theta0 = theta[site0]
    metrics['thetaExtr'][iSite] = theta0 #K
    
    lat0, lon0 = mesh.get_latLon_inds(site0)
    metrics['latExtr'][iSite] = lat0*r2d #degrees
    metrics['lonExtr'][iSite] = lon0*r2d
    
    ampMaxMin = calc_amplitude_maxMin(basinInds, theta)
    metrics['ampMaxMin'][iSite] = ampMaxMin #K
    
    #areal measures -------------------------------
    area = calc_area(basinInds, mesh)
    #metrics['area'][iSite] = area/1.e6 #km^2
    rEquiv = np.sqrt(area/np.pi) #pi r^2 = A
    metrics['rEquiv'][iSite] = rEquiv/1.e3 #km
    
    circ = calc_circulation(basinInds, vort, mesh)
    metrics['circ'][iSite] = circ/1.e6 #km^2/s
    
    vortMean = circ/area
    metrics['vortMean'][iSite] = vortMean #1/s
    
    thetaVol = calc_fieldVolume(basinInds, theta, theta0, mesh)
    metrics['thetaVol'][iSite] = thetaVol #K*km^2
    
    ampArea = thetaVol/(area/1.e6)
    metrics['ampMean'][iSite] = ampArea #K
  
  return metrics
    
def run_metrics(fNameOut, info, mesh, dataMetr, dataSeg, iTimeStart, iTimeEnd):
  
  nSitesMax = len(dataSeg.dimensions['nMax'])
  dataMetrics = write_netcdf_header(fNameOut, info, nSitesMax)
  
  for iTime in xrange(iTimeStart, iTimeEnd+1):
    theta = dataMetr.variables['theta'][iTime,:]
    vort = dataMetr.variables['vort'][iTime,:]
    
    cell2Site = dataSeg.variables['cell2Site'][iTime,:]
    sitesMin = dataSeg.variables['sitesMin'][iTime,:];
    nMin = dataSeg.variables['nSitesMin'][iTime]; sitesMin = sitesMin[0:nMin]
    sitesMax = dataSeg.variables['sitesMax'][iTime,:]; 
    nMax = dataSeg.variables['nSitesMax'][iTime]; sitesMax = sitesMax[0:nMax]
    
    sites = np.concatenate((sitesMin,sitesMax))
    
    metrics = calc_metrics(sites, cell2Site, vort, theta, mesh)
    
    write_netcdf_iTime(dataMetrics, iTime, metrics, sites)
  
  dataMetrics.close()
  
def write_netcdf_header(fName, info, nSitesMax):
  #I don't know how to make ragged arrays, so we'll use 
  #array[nTimes,nMax] and nElts[nTimes]
  
  data = netCDF4.Dataset(fName, 'w', format='NETCDF4')
  data.description = info
  
  # dimensions
  data.createDimension('time', None)
  data.createDimension('nMax', nSitesMax)

  # variables
  data.createVariable('nSites', 'i4', ('time',))
  data.createVariable('sites', 'i4', ('time','nMax',))
  
  for key in metricKeys:
    data.createVariable(key, 'f8', ('time','nMax',))
  return data
  
def write_netcdf_iTime(data, iTime, metrics, sites):
  # fill file. with time as unlimited, dimension will just keep growing
  
  nSites = len(sites)
  data.variables['nSites'][iTime] = nSites
  data.variables['sites'][iTime,0:nSites] = sites[0:nSites]
  
  for key in metrics:
    data.variables[key][iTime,0:nSites] =  metrics[key][0:nSites]
#

def print_metrics(fName):
  data = netCDF4.Dataset(fName,'r')
  
  nTimes = len(data.dimensions['time'])
  for iTime in xrange(nTimes):
    nSites = data.variables['nSites'][iTime]
    for key in metricKeys:
      vals = data.variables[key][iTime,:]; vals = vals[0:nSites]
      print "Time ", iTime, key; print vals
  
  data.close()

def get_metrics_basin(data, iTime, site):
  
  vals = []
  sites = data.variables['sites'][iTime,:]
  iSite = np.where(sites==site)[0][0]
  
  for key in metricKeys:
    vals.append(data.variables[key][iTime,iSite])
    
  return vals

def calc_diff_metricSpace(data, iTime0, site0, iTime1, sites1, r):
  #here, we're calculating something like a metric in a metric space so value >= 0.
  #Input sites1 as list,array,...index-able
  #return "distance" for >=1 basins from a reference basin.
  #measures can include: distance, diffArea, diffIntensity,...
  
  #ansatz: basin metrics are useful for discrimating differences (ie, tpv properties should be semi-persistent)
  #diffKeys = ['thetaExtr', 'vortMean', 'rEquiv', 'latExtr']; nKeys = len(diffKeys)
  #refDiffs = [3.0, 5.e-5, 300.0, 5.0]
  diffKeys = ['thetaExtr', 'latExtr']; nKeys = len(diffKeys)
  refDiffs = [1.0, 2.0]
  
  #get site0 and sites1 measures ---------
  
  #site0
  vals0 = np.empty(nKeys, dtype=float)
  
  iTime = iTime0
  sites = data.variables['sites'][iTime,:]
  iSite = np.where(sites==site0)[0][0]
  for iKey in xrange(nKeys):
    key = diffKeys[iKey]
    val = data.variables[key][iTime,iSite]
    vals0[iKey] = val
    
  #sites1
  nSites1 = len(sites1)
  vals1 = np.empty((nKeys, nSites1),dtype=float) #ordered so columns hold keyValues for all sites (better cache?)
  
  iTime = iTime1
  sites = data.variables['sites'][iTime,:]
  for iKey in xrange(nKeys):
    key = diffKeys[iKey]
    vals = data.variables[key][iTime,:]
    for iSite1 in xrange(nSites1):
      site1 = sites1[iSite1]
      ind = np.where(sites==site1)[0][0]
      val = vals[ind]
      vals1[iKey, iSite1] = val
    #end iSite1  
  #end iKey
  
  #calc difference measures --------------------------
  diffMeasures = np.empty((nKeys, nSites1), dtype=float)
  
  for iKey in xrange(nKeys):
    diffMeasures[iKey,:] = np.absolute(vals1[iKey,:]-vals0[iKey])
  
  #print vals0; print vals1; 
  #print diffMeasures
  
  #return "distances" in metric space -------------------
  #imagine that dTheta~2,5,10,25K, distance~60,100,300,1000km, vortMean~1.e-5...so want a way to aggregate.
  #1) d = sum( diffMeasures_i/mean(diffMeasures_i) ) gives normalized deviations. But, consider if diffMeasures_i
  #is all within 1.e-12. Then, the values are "really" all the same, but normalizing by mean(diffMeasures_i) will amplify
  #the differences.
  #2) have some equivalent conversion value between variables (5K ~ 100km ~...). Deciding on values is arbitrary, but matches
  #may not be too sensitive to the specific choices (hopefully?)
  d = np.zeros(nSites1, dtype=float)
  for iKey in xrange(nKeys):
    d += diffMeasures[iKey,:]/refDiffs[iKey]
  
  if (True): #print out some info for basins near pole
    indLat = diffKeys.index('latExtr')
    if (vals0[indLat]>85.):
      print "-------------------\n"
      print site0, vals0
      print sites1, vals1
      print d
      print "-------------------\n"
   
  return d


