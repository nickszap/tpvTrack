import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import helpers
import llMesh

#watershed has a few options for implementation:
#-for every cell, walk down steepest gradient to the basin

def find_minCells_region_flat(vals, nLat, nLon, inRegion):
  #return array[nCells] with 1 if cell is min and cell in region
  
  isMin = np.zeros((nLat*nLon),dtype=int)
  for iLat  in xrange(nLat):
    for iLon in xrange(nLon):
      ind0 = index_2dTo1d(iLat, iLon, nLon)
      if (inRegion[ind0]<1):
        continue
      
      nbrs = nbrInds_ll_flat(iLat, iLon, nLat, nLon)
      val0 = vals[ind0]
      valNbrs = vals[nbrs]
      if (np.all(val0<=valNbrs)): #is site if can't descend from it
        isMin[ind0] = 1

  nMin = np.sum(isMin)
  print "Number of local min in region: ", nMin
  return isMin

def watershed_region(vals, cellIsMin, nLat, nLon, r, dRegion, latCell, lonCell, inRegion):
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
  
  nCells = nLat*nLon
  cell2Site = -np.ones(nCells,dtype=int) #so no cell2Site[iCell]=iCell
  for iCell in xrange(nCells):
    if (inRegion[iCell]<1):
      continue
    if (cellIsMin[iCell]>0):
      cell2Site[iCell]= iCell

  #get local steepest path
  for iCell in xrange(nCells):
    if (inRegion[iCell]<1):
      continue
    if (cellIsMin[iCell]>0):
      continue
    
    iLat0, iLon0 = index_1dTo2d(iCell, nLon)
    nbrInds_lat, nbrInds_lon = nbrInds_ll(iLat0, iLon0, nLat, nLon)
    nbrInds = index_2dTo1d(np.array(nbrInds_lat), np.array(nbrInds_lon), nLon)
    #print iLat0, iLon0, nbrInds_lat; print nbrInds_lon; print nbrInds
    #nNbrs = len(nbrInds_lat)
    
    val0 = vals[iCell]
    valNbrs = vals[nbrInds]

    #correspondence is towards minimum gradient.
    lat0 = latCell[iLat0]; lon0 = lonCell[iLon0]
    latNbrs = latCell[nbrInds_lat]; lonNbrs = lonCell[nbrInds_lon]
    dNbrs = calc_distSphere_multiple(r, lat0, lon0, latNbrs, lonNbrs)
    dMin = r/1.e16; dNbrs[dNbrs<dMin]=dMin #avoid divide by 0
    #print valNbrs, dNbrs, val0
    valNbrs = (valNbrs-val0)/dNbrs
    iNbr = np.argmin(valNbrs)
    cell2Site[iCell] = nbrInds[iNbr]
  
  nRedirect = 0
  
  #Filter local extrema by area to limit high (spatial) frequency "noise".
  #For multiple close mins, the smallest counts as min for that region.
  #dRegion = 300.e3 #radius in meters of disk of filtering region
  nLatIndsLength = calc_latIndicesWithinLength(nLat, r, dRegion)
  nLonIndsLength = calc_lonIndicesWithinLength(latCell, nLon, r, dRegion)
    
  for iLat in xrange(nLat):
    iLonRef = 0
    inDiskLat, inDiskLon_ref = gatherInds_region_latBox_1AtPole(iLat, iLonRef, nLat, nLon, latCell, lonCell,
                                                nLatIndsLength, nLonIndsLength, r, dRegion)
    #
    for iLon in xrange(nLon):
      iCell = index_2dTo1d(iLat, iLon, nLon)
      if (inRegion[iCell]<1):
        continue
      if (cellIsMin[iCell]>0):
        #see if cell is min in region, not just neighbors.
        #if not regional min, update cell2Site so local min goes to another basin
        diffLonInd = iLon-iLonRef
        inDiskLon = (inDiskLon_ref+diffLonInd)%nLon

        cellsRegion = index_2dTo1d(inDiskLat, inDiskLon, nLon)
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
  for iCell in xrange(nCells):
    if (inRegion[iCell]<1):
      #what should these correspond to?
      continue
    nextCell = cell2Site[iCell]
    nCount = 0
    while (not cellIsMin[nextCell]>0):
      nextCell = cell2Site[nextCell]
      #print "Cell {0} going to {1}".format(iCell,nextCell); print vals[iCell], vals[nextCell]
      nCount=nCount+1
      if (nCount>nCells):
        iLat, iLon = index_1dTo2d(iCell, nLon)
        print "Uhoh, stuck in while loop for iLat,iLon ({0},{1}) with value {2}".format(iLat, iLon, vals[iCell])
        break

    cell2Site[iCell] = nextCell

  return (cell2Site, cellIsMin)

def segment_high_low_watershed_region(lat, lon, theta, vort, r, inRegion, dRegion):
  #get high and low basin seeds, associate cells to both high and low basins if not extrema.
  #to decide whether "really" part of high or low basin, we have options:
  #-(anti-)cyclonic for (high) low...is local vorticity noisy?
  #-closer theta value to maxima a la color scale grouping...huge min or max value now matters
  #-whether steeper gradient is to high or low
  #-physical distance
  #-concavity of surface a la last closed contour
  
  #theta, vort, inRegion come in as 1d arrays
  #use inRegion for (1) region of interest (2)ignore missing values
  
  nLat = len(lat); nLon = len(lon); nCells = nLat*nLon
  
  #segment ------------------------------------------
  #dRegion = 300.e3 #radius in same units as r of disk for min basin
  
  #mins
  print "Finding minima"
  #to adapt global watershed to region, make values outside of region huge so don't steepest descend that way
  bigVal = 1.e10
  vals = np.copy(theta) #so don't affect variable passed in
  vals[inRegion<1] = bigVal
  
  cellIsMin = find_minCells_region_flat(vals, nLat, nLon, inRegion)
  cell2SiteMin, cellIsMin = watershed_region(vals, cellIsMin, nLat, nLon, r, dRegion, lat, lon, inRegion)
  
  #maxs: perform min on an inverted surface
  print "Finding maxima"
  #adapt global watershed to region
  vals = -np.copy(theta)
  vals[inRegion<1] = bigVal
  
  cellIsMax = find_minCells_region_flat(vals, nLat, nLon, inRegion)
  cell2SiteMax, cellIsMax = watershed_region(vals, cellIsMax, nLat, nLon, r, dRegion, lat, lon, inRegion)
  
  #"voting" procedure for low/high classification ------
  print "Associating to max or min"
  cell2Site = -np.ones(nCells, dtype=int)
  
  classifyMethod = 0 #see if-checks below for corresponding number for method
  if (classifyMethod==0):
    #local vorticity
    print("Classifying by local vorticity")
    for iCell in xrange(nCells):
      if (inRegion[iCell]<1):
        continue
      #cyclonic depends on hemisphere
      if (cellIsMin[iCell]>0 or cellIsMax[iCell]>0): #allows for cyclonic max. is that right?
        cell2Site[iCell] = iCell
      else:
        signHem = 1 #sign function is problematic since sign(0)=0
        iLat0, iLon0 = index_1dTo2d(iCell, nLon)
        if (lat[iLat0]<0): #lat=0 gets put in NH
          signHem = -signHem
        if (signHem*vort[iCell]<0): #anticyclonic
          cell2Site[iCell] = cell2SiteMax[iCell]
        else: #0 or cyclonic
          cell2Site[iCell] = cell2SiteMin[iCell]
  else:
    print "Error. Classification method not recognized/specified."
          
  return (cell2Site, cellIsMin, cellIsMax)

def demo_segment_steven():
  #Seeing what's happening with Steven's NH file
  #Latitude is reversed
  
  rEarth = 6370.e3; dFilter = 300.e3 #radius in same units as r of disk for min basin
  fDirSave = '/data02/cases/test_segment/stevenCase/'
  
  fDirData = '/home/scavallo/'
  fnames = sorted(glob.glob(fDirData+'erainterim_pv_201309*.nc'), key=os.path.getmtime)
  print fnames
  for iFile in xrange(len(fnames)):
    #gather persistent info like mesh, times,... ----------------------
    fname = fnames[iFile]; 
    data = netCDF4.Dataset(fname,'r')
    
    nTimes = len(data.dimensions['time']); times = range(nTimes)
    d2r = np.pi/180.; 
    lat = data.variables['xlat'][:,0]*d2r; lon = data.variables['xlong'][0,:]*d2r
    nLat = len(lat); nLon = len(lon)
    lat = lat[::-1] # reverse so north pole is lat[0]
    
    iLevel = 3; #for 2 pvu
    print "PVU level: ", data.variables['levels'][iLevel]
    
    #specify region for segmentation, ----------------------------
    #including halo so make sure nbr values are decent
    latThreshHalo = 43.*np.pi/180.
    latThresh = 45.*np.pi/180.
    inRegionHalo = np.zeros((nLat,nLon), dtype=int); inRegionHalo[lat>latThreshHalo,:] = 1
    inRegion = np.zeros((nLat,nLon), dtype=int); inRegion[lat>latThresh,:] = 1
    
    #loop over individual times ------------------------------
    for iTime in xrange(nTimes):
      theta = data.variables['th'][iTime,iLevel,:,:]
      u = data.variables['u'][iTime,iLevel,:,:]; v = data.variables['v'][iTime,iLevel,:,:]
      
      #reorient latitudes
      theta = theta[::-1,:]; u = u[::-1,:]; v = v[::-1,:]
      
      #fill missing values
      isMissing = np.isnan(theta)
      u = fill_missingVals_region(u, nLat, nLon, isMissing, inRegionHalo)
      v = fill_missingVals_region(v, nLat, nLon, isMissing, inRegionHalo)
      theta = fill_missingVals_region(theta, nLat, nLon, isMissing, inRegionHalo)
      
      #calc vorticity
      vort = calc_vertVorticity_ll(u, v, nLat, nLon, lat, rEarth)
    
      #segment ------------------------------------------
      theta = flatten_2dTo1d(theta, nLat, nLon)
      vort = flatten_2dTo1d(vort, nLat, nLon)
      inRegion = flatten_2dTo1d(inRegion, nLat, nLon)
      cell2Site, cellIsMin, cellIsMax = segment_high_low_watershed_region(lat, lon, theta, vort, rEarth, inRegion, dFilter)
    
      #save segmented result
      t0 = int(times[iTime])
      fInfo = eraTimeToCalendarTime(t0)
      fNameSave = 'seg_'+fInfo+'.npz'
      f = fDirSave+fNameSave
      print "Saving cell2Site to file "+f
      np.savez(f, cell2Site=cell2Site, cellIsMin=cellIsMin, cellIsMax=cellIsMax)
      fNameSave = 'fields_'+fInfo+'.npz'
      f = fDirSave+fNameSave
      print "Saving cell2Site data to file "+f
      np.savez(f, lat=lat, lon=lon, u=u, v=v, theta=theta, inRegion=inRegion)
    
      #plot segmented result ----------
      if (True):
      #if (False):
        print "Plotting segmentation"
        thetaFlat = theta; #flatten_2dTo1d(theta, nLat, nLon)
        #print thetaFlat[cellIsMin>0]; print thetaFlat[cellIsMax>0]
        isMin = unflatten_1dTo2d(cellIsMin, nLat, nLon)
        isMax = unflatten_1dTo2d(cellIsMax, nLat, nLon)
        vals = thetaFlat[cell2Site];
        vals = unflatten_1dTo2d(vals, nLat, nLon)
        #plot_segment(lat, lon, vals, isMin, isMax)
        fNameSave = fDirSave+fNameSave+'.png'
        print "Saving segmentation plot to: "+fNameSave
        plot_segment_save(fNameSave, lat, lon, vals, isMin, isMax)
    #looped over all times in this file
    data.close()

if __name__ == '__main__':
  #example_segment()
  #demo_plot_segInds()
  #demo_plot_segSites()
  #demo_segment()
  #demo_segment_era()
  demo_segment_steven()
