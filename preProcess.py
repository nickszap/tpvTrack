'''
Take input data and write mesh+metr info to our format.
mesh: lat,lon,areaCell
metr: u,v,theta,verticalVorticity,inRegion on DT
'''

import numpy as np
import netCDF4

def get_missingCells_file(data):
  #if search vertical column down from top, there can be no 2pvu value if:
  #-entire column is above 2pvu: DT below sfc
  #-entire column is below 2pvu: wrong hemisphere or low pv column (tropics, anticyclonic,...)
  
  pvInd = 1
  isMissing = data.variables['TMP_P0_L109_GLL0'][pvInd,:,:].mask
  return isMissing
  
def fill_missingVals_region(valsIn, nLat, nLon, isMissing, inRegion):
  #fill value is average of non-missing neighbors
  
  vals = np.copy(valsIn)
  needFill = isMissing*inRegion;
  nNeedFill = np.sum(needFill); print "Filling {0} values".format(nNeedFill)
  
  while (np.sum(needFill)>0):
    for iLat in xrange(nLat):
      for iLon in xrange(nLon):
        if (needFill[iLat,iLon]>0): #True>0 is True, False==0 is True
          nbrInds_lat, nbrInds_lon = nbrInds_ll(iLat, iLon, nLat, nLon)
          nbrsNeedFill = needFill[nbrInds_lat, nbrInds_lon]
          if (False in nbrsNeedFill): #have neighbor with value
            #fill value is average of valid nbrs
            validNbrs = nbrsNeedFill==False;
            valsNbrs = vals[nbrInds_lat, nbrInds_lon]
            vals[iLat, iLon] = np.mean(valsNbrs[validNbrs])+1.e-10 #so don't have same value
            needFill[iLat,iLon]=False
  
  return vals

def demo_eraI(fMesh, fData):
  fDirData = '/data02/cases/2006/eraI/pv/'
  fnames = sorted(glob.glob(fDirData+'eraI_theta-u-v_2pvu_2006-07-20*.nc'), key=os.path.getmtime)
  nFiles = len(fData)
  for iFile in xrange(len(fnames)):
    #gather persistent info like mesh, times,... ----------------------
    fname = fnames[iFile]; 
    data = netCDF4.Dataset(fname,'r')
    
    times = data.variables['time'][:]; nTimes = len(times)
    d2r = np.pi/180.; 
    lat = data.variables['latitude'][:]*d2r; lon = data.variables['longitude'][:]*d2r
    nLat = len(lat); nLon = len(lon)
    
    #specify region for segmentation, ----------------------------
    #including halo so make sure nbr values are decent
    latThreshHalo = 43.*np.pi/180.
    latThresh = 45.*np.pi/180.
    inRegionHalo = np.zeros((nLat,nLon), dtype=int); inRegionHalo[lat>latThreshHalo,:] = 1
    inRegion = np.zeros((nLat,nLon), dtype=int); inRegion[lat>latThresh,:] = 1
    
    #loop over individual times ------------------------------
    for iTime in xrange(nTimes):
      theta = data.variables['pt'][iTime,:,:]
      u = data.variables['u'][iTime,:,:]; v = data.variables['v'][iTime,:,:]
      vort = calc_vertVorticity_ll(u, v, nLat, nLon, lat, rEarth)
    
      #segment ------------------------------------------
      theta = flatten_2dTo1d(theta, nLat, nLon)
      vort = flatten_2dTo1d(vort, nLat, nLon)
      inRegion = flatten_2dTo1d(inRegion, nLat, nLon)

def writeNetcdf_mesh(fNameOut, info, lat,lon, inRegion, areaCell):

  data = netCDF4.Dataset(fNameOut, 'w', format='NETCDF4')
  data.description = 'Mesh information for: '+info

  # dimensions
  nCells = len(areaCell)
  data.createDimension('nCells', nCells)

  # variables
  dataLat = data.createVariable('lat', 'f4', ('nCells',))
  dataLon = data.createVariable('lon', 'f4', ('nCells',))
  dataInRegion = data.createVariable('inRegion', 'i2', ('nCells',))
  dataAreaCell = data.createVariable('areaCell', 'f4', ('nCells',))
  
  #units and descriptions
  dataLat.units = 'radians N'
  dataLon.units = 'radians E'
  dataInRegion.units = 'flag'
  dataAreaCell.units = 'm^2'
  
  # fill data
  dataLat[:] = lat[:]
  dataLon[:] = lon[:]
  dataInRegion[:] = inRegion[:]
  dataAreaCell[:] = areaCell[:]

  data.close()

