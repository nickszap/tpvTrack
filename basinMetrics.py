import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
import glob
import datetime as dt
import os
#from scipy import stats

import segment_ll_flat
import track_ll as track
import sys; sys.path.append("/home/nickszap/Dropbox/pythonScripts/minAreaBoundingRectangle/")
import bbox_test

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
#isn't straightforward
r = 6370.e3
r2d = 180./np.pi

def calc_aspectRatio(site, cell2Site, lat, lon):
  #project all points onto a plane. Then use properties of a bounding box fit around cells.
  
  if (np.max(lat)>20):
    print "Uhoh. Check that lat, lon input in radians"
  
  nLat = len(lat); nLon = len(lon)
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  latBasin = lat[iLat]; lonBasin = lon[iLon]; nPts = len(latBasin)
  
  lat0 = np.median(latBasin)*r2d; lon0 = np.median(lonBasin)*r2d
  m = Basemap(projection='ortho',lon_0=lon0,lat_0=lat0, resolution='l')
  x, y = m(lonBasin*r2d, latBasin*r2d);
  
  if (False):
    plt.figure()
    m.drawcoastlines()
    m.drawmapboundary()
    m.scatter(x, y, s=1, marker='o')
    plt.show()
  
  xyPoints = np.empty((nPts,2), dtype=float)
  xyPoints[:,0] = x[:]; xyPoints[:,1] = y[:]
  
  maxL, minL, rotAngle = bbox_test.runCase(xyPoints, False)
  print 'min/max radii of bounding box (km): ', minL/2.e3, maxL/2.e3
  return maxL/minL
  
def calc_majorAsymmetry(site, cell2Site, lat, lon, areaLatCell):
  #let the longest dimension of the bounding box be the major axis.
  #then have areas (+ and -) to both sides of extremum.
  #return ratio
  
  if (np.max(lat)>20):
    print "Uhoh. Check that lat, lon input in radians"
  
  nLat = len(lat); nLon = len(lon)
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  latBasin = lat[iLat]; lonBasin = lon[iLon]; nPts = len(latBasin)
  
  lat0 = np.median(latBasin)*r2d; lon0 = np.median(lonBasin)*r2d
  m = Basemap(projection='ortho',lon_0=lon0,lat_0=lat0, resolution='l')
  x, y = m(lonBasin*r2d, latBasin*r2d);
  
  xyPoints = np.empty((nPts,2), dtype=float)
  xyPoints[:,0] = x[:]; xyPoints[:,1] = y[:]
  maxL, minL, rotAngle = bbox_test.runCase(xyPoints, False)
  
  #project areas onto line centered at basin extremum
  iLat0, iLon0 = segment_ll_flat.index_1dTo2d(site, nLon)
  x0, y0 = m(lon[iLon0]*r2d, lat[iLat0]*r2d);
  dxLine = np.cos(rotAngle); dyLine = np.sin(rotAngle)
  #just need sign of dot product
  d = (x-x0)*dxLine+(y-y0)*dyLine
  sgnD = np.sign(d)
  
  if (True):
    plt.figure()
    plt.scatter(x,y)
    meanD = np.max(x-x0); plt.plot([x0,x0+dxLine*meanD], [y0,y0+dyLine*meanD], 'r')
    #m.drawcoastlines()
    #m.drawmapboundary()
    #m.scatter(x, y, s=1, marker='o')
    plt.show()
  
  areaCells = areaLatCell[iLat]/1.e6
  areaNeg = np.sum(areaCells[sgnD<0]); areaPos = np.sum(areaCells[sgnD>0])
  return max(areaNeg,areaPos)/min(areaNeg, areaPos)
  
def calc_circulation(site, cell2Site, vort, areaLatCell, nLat, nLon):
  #\int gradxu . dA
  
  #lat/lon cells in each basin
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  
  vortCells = vort[iLat, iLon]
  cellAreas = areaLatCell[iLat]
  return np.dot(vortCells, cellAreas)

def calc_secondMomentArea(site, cell2Site, areaLatCell, lat, lon):  
  #moment of inertia about extremum: \int r^2 dA
  #if want rotational moment, need some form of mass for sum(m r^2)
  #return in units of km^4
  
  nLat = len(lat); nLon = len(lon)
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  
  #moment about extremum
  iLat0, iLon0 = segment_ll_flat.index_1dTo2d(site, nLon)
  lat0 = lat[iLat0]; lon0 = lon[iLon0]
  
  d = segment_ll_flat.calc_distSphere_multiple(r/1.e3, lat0, lon0, lat[iLat], lon[iLon])
  #d /= 1.e3
  cellAreas = areaLatCell[iLat]/1.e6
  J = np.dot(d*d, cellAreas)
  return J
  
def calc_amplitude(site, cell2Site, thetaFlat):
  #max-min value
  
  inBasin = cell2Site==site;
  if (not inBasin.any()):
    print "Uhoh, site not in cell2Site (cell,value): "+str(site)+' '+str(thetaFlat[site])
    return 0.0
  valsBasin = thetaFlat[inBasin]
  minVal = np.min(valsBasin); maxVal = np.max(valsBasin)
  return maxVal-minVal

def calc_fieldVolume(site, cell2Site, areaLatCell, thetaFlat, valRef, nLat, nLon):  
  #"volume": \int (theta-thetaRef) dA
  #if want rotational moment, need some form of mass for sum(m r^2)
  #return in units of km^4
  
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  
  #moment about extremum
  vals = thetaFlat[inBasin]
  #valRef = np.max(vals)
  vals = valRef-vals
  cellAreas = areaLatCell[iLat]/1.e6
  J = np.dot(vals, cellAreas)
  return J

def calc_amplitudeBoundary(site, cell2Site, nBoundaryFlat, thetaFlat):
  #<boundary>-extremum value
  
  vals = calc_boundaryValues(site, cell2Site, nBoundaryFlat, thetaFlat)
  valb = np.median(vals) #could do area weighted mean, mode, min,...
  #valb = stats.mode(vals)
  val0 = thetaFlat[site]
  return valb-val0
  
def calc_area(site, cell2Site, areaLatCell, nLat, nLon):
  
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  
  return np.sum(areaLatCell[iLat])

def calc_boundaryLength(site, cell2Site, nBoundaryFlat, dLatCell, nLat, nLon):
  
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  inBasin *= nBoundaryFlat>0
  
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  #nBoundaryBasinCell = nBoundaryFlat[inBasin]
  return np.sum(dLatCell[iLat])
  #return np.dot(nBoundaryBasinCell, dLatCell[iLat])

def calc_minMaxDistToBoundary(site, cell2Site, nBoundaryFlat, lat, lon):
  
  nLat = len(lat); nLon = len(lon)
  nCells = nLat*nLon
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  inBasin *= nBoundaryFlat>0
  
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  iLat, iLon = segment_ll_flat.index_1dTo2d(indsBasin, nLon)
  iLat0, iLon0 = segment_ll_flat.index_1dTo2d(site, nLon)
  
  latb = lat[iLat]; lonb = lon[iLon]
  lat0 = lat[iLat0]; lon0 = lon[iLon0]
  
  d = segment_ll_flat.calc_distSphere_multiple(r, lat0, lon0, latb, lonb)
  minD = np.min(d); maxD = np.max(d)
  
  return (minD, maxD)
  
def calc_boundaryValues(site, cell2Site, nBoundaryFlat, thetaFlat):
  
  inBasin = cell2Site==site; #print '# of cells in basin: ', np.sum(inBasin)
  inBasin *= nBoundaryFlat>0
  return thetaFlat[inBasin]

def find_boundaryCells(cell2Site, nLat, nLon, inRegion):
  #return isBoundary[lat,lon] >0  if cell2Site[nbr] != cell2Site[self]
  
  isBoundary = np.zeros((nLat, nLon), dtype=int)
  for iLat0 in xrange(nLat):
    for iLon0 in xrange(nLon):        
      ind0 = segment_ll_flat.index_2dTo1d(iLat0, iLon0, nLon)
      if (inRegion[ind0]<=0):
        continue;
      val0 = cell2Site[ind0]
      
      indNbrs = segment_ll_flat.nbrInds_ll_flat(iLat0, iLon0, nLat, nLon)
      valNbrs = cell2Site[indNbrs]
      
      nDiff = np.sum(valNbrs!=val0)
      isBoundary[iLat0, iLon0] = nDiff
      
  return isBoundary

def find_boundaryCells_basin(site, cell2Site, nLat, nLon):
  #return lat/lon indices of cells in basin with cell2Site[nbr] != cell2Site[self]
  
  nCells = nLat*nLon
  inBasin = cell2Site==site
  indsBasin = np.arange(nCells, dtype=int)[inBasin]
  isBoundary = np.zeros((nLat, nLon), dtype=int)
  for iCell in indsBasin:
    iLat0, iLon0 = segment_ll_flat.index_1dTo2d(iCell, nLon)
    indNbrs = segment_ll_flat.nbrInds_ll_flat(iLat0, iLon0, nLat, nLon)
    valNbrs = cell2Site[indNbrs]
    
    nDiff = np.sum(valNbrs!=site)
    isBoundary[iLat0, iLon0] = nDiff
  
  return isBoundary

def calc_distToSites(site0, allSites, lat, lon):
  
  nLat = len(lat); nLon = len(lon)
  
  iLat0, iLon0 = segment_ll_flat.index_1dTo2d(site0, nLon)
  iLat, iLon = segment_ll_flat.index_1dTo2d(allSites[allSites!=site0], nLon)
  
  latb = lat[iLat]; lonb = lon[iLon]
  lat0 = lat[iLat0]; lon0 = lon[iLon0]
  
  d = segment_ll_flat.calc_distSphere_multiple(r, lat0, lon0, latb, lonb)
  return d

def calc_basinMetrics(fSeg, fMetr):
  #pass in 0 or nSites filesSeg and filesMetr.
  #calculate basin properties and return metrics as dictionary
  
  #gather mesh+metr info
  data = np.load(fMetr)
  lat = data['lat']; lon = data['lon']; nLat = len(lat); nLon = len(lon)
  u = data['u']; v = data['v']; thetaFlat = data['theta']
  data.close()
  
  #gather seg info
  data_seg = np.load(fSeg)
  cell2Site = data_seg['cell2Site'][:]
  cellIsMin = data_seg['cellIsMin'][:]
  data_seg.close()
  sites = cell2Site[cellIsMin>0]; nSites = len(sites)
  
  #initialize
  nSites = len(sites)
  metricKeys = 'circ vortMean amp area rEquiv aspectRatio thetaVol ampArea thetaMin latMin lonMin'.split()
  metrics = {}
  for key in metricKeys:
    metrics[key] = np.empty(nSites, dtype=float)
  
  #start calculating fields
  areaLatCell = track.calc_areaLatStrips(lat, r)/nLon
  vort = segment_ll_flat.calc_vertVorticity_ll(u, v, nLat, nLon, lat, r)
  
  for iSite in xrange(nSites):
    site0 = sites[iSite]
    #inTPV = cell2Site == site0
    
    thetaMin = thetaFlat[site0]
    metrics['thetaMin'][iSite] = thetaMin
    
    iLat, iLon = segment_ll_flat.index_1dTo2d(site0, nLon)
    latMin = lat[iLat]*r2d; lonMin = lon[iLon]*r2d
    metrics['latMin'][iSite] = latMin #degrees
    metrics['lonMin'][iSite] = lonMin
    
    #areal measures -------------------------------
    circ = calc_circulation(site0, cell2Site, vort, areaLatCell, nLat, nLon);
    metrics['circ'][iSite] = circ/1.e6 #km^2/s

    amp = calc_amplitude(site0, cell2Site, thetaFlat)
    metrics['amp'][iSite] = amp #K
    
    area = calc_area(site0, cell2Site, areaLatCell, nLat, nLon)
    metrics['area'][iSite] = area/1.e6 #km^2
    
    rEquiv = np.sqrt(area/np.pi) #pi r^2 = A
    metrics['rEquiv'][iSite] = rEquiv/1.e3 #km
    
    vortMean = circ/area
    metrics['vortMean'][iSite] = vortMean #1/s
    
    if (False):
      aspectRatio = calc_aspectRatio(site0, cell2Site, lat, lon)
    else:
      aspectRatio = -999.  
    metrics['aspectRatio'][iSite] = aspectRatio
    
    #metrics for extremum environment ----------------------------- 
    thetaVol = calc_fieldVolume(site0, cell2Site, areaLatCell, thetaFlat, thetaMin, nLat, nLon)
    metrics['thetaVol'][iSite] = thetaVol #K*m^2
    
    ampArea = thetaVol/(area/1.e6)
    metrics['ampArea'][iSite] = ampArea #K
  
  return metrics
    
def demo_calcTPVMetrics():
  
  dirData = '/data02/cases/summer2006/eraI/'
  filesMetr = sorted(glob.glob(dirData+'/seg/fields_*.npz'), key=os.path.getmtime)
  filesSeg = sorted(glob.glob(dirData+'/seg/seg_*.npz'), key=os.path.getmtime)
  dirSave = dirData+'/seg/'
  
  nFiles = len(filesSeg)
  print "Saving metrics files in "+dirSave
  for iFile in xrange(nFiles): #looping over files == looping over times
    
    fMetr = filesMetr[iFile]
    fSeg = filesSeg[iFile]
    
    metrics = calc_basinMetrics(fSeg, fMetr)
    
    fSave = 'metrics_'+fSeg.split('/')[-1] #+'.npz'
    fSave = dirSave+fSave; print "Saving file: "+fSave
    
    np.savez(fSave,
    circ=metrics['circ'], vortMean=metrics['vortMean'], amp=metrics['amp'], 
    area=metrics['area'], rEquiv=metrics['rEquiv'], aspectRatio=metrics['aspectRatio'],
    thetaVol=metrics['thetaVol'], ampArea=metrics['ampArea'], thetaMin=metrics['thetaMin'])

def plotTimeSeries(dates0, vals0, yInfo):
  
  plt.figure()
  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
  plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=4))
  plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
  
  plt.plot(dates0, vals0)
  
  plt.xlabel('Date')
  plt.ylabel(yInfo)
  
  plt.gcf().autofmt_xdate()

def plot_basinCompare(dates0, vals0, label0, dates1, vals1, label1, yInfo):
  
  plt.figure() #---------------------------
  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
  plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=4))
  plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
  
  plt.plot(dates0, vals0, 'b', label=label0)
  plt.plot(dates1, vals1, 'r', label=label1)
  
  plt.xlabel('Date')
  plt.ylabel(yInfo)
  plt.legend(loc='best')
  
  plt.gcf().autofmt_xdate()
  
  plt.show()
  
if __name__=='__main__':
  #demo_trackBasin()
  #demo_compareTPVs()
  demo_calcTPVMetrics()
  
