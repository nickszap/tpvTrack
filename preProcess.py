'''
Take input data and write mesh+metr info to our format.
mesh: lat,lon,areaCell
metr: u,v,theta,verticalVorticity,inRegion on DT
'''

import numpy as np
import netCDF4

import helpers
import llMesh
import mpasMesh

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
  
def get_segmentVars_file(data):
  #return data of variables needed from file, no time index.
  #return SI units
  
  #fname = '/data02/cases/2014/gfs_4_20140101_0000_123.nc'
  #data = netCDF4.Dataset(fname,'r')
  lat = data.variables['lat_0'][:] * np.pi/180. #deg 2 radians
  lon = data.variables['lon_0'][:] * np.pi/180.
  
  #for pvu levels, lv_PVL4 = [-2, 2]*E-6
  pvInd = 1
  tmp = data.variables['TMP_P0_L109_GLL0'][pvInd,:,:].data #K
  press = data.variables['PRES_P0_L109_GLL0'][pvInd,:,:].data #Pa
  
  u = data.variables['UGRD_P0_L109_GLL0'][pvInd,:,:].data #m/s
  v = data.variables['VGRD_P0_L109_GLL0'][pvInd,:,:].data
  
  return (lat, lon, u, v, tmp, press)

def calc_vertVorticity_ll(u, v, nLat, nLon, lat, r):
  '''
  Pulled from: http://www.ncl.ucar.edu/Document/Functions/Built-in/uv2vr_cfd.shtml :
  According to H.B. Bluestein [Synoptic-Dynamic Meteorology in Midlatitudes, 1992, 
  Oxford Univ. Press p113-114], 
  let D represent the partial derivative, a the radius of the earth, 
  phi the latitude and dx2/dy2 the appropriate longitudinal and latitudinal spacing, 
  respectively. Then, letting j be the latitude y-subscript, and i be the longitude x-subscript:

    rv = Dv/Dx - Du/Dy + (u/a)*tan(phi)


    rv(j,i) = (v(j,i+1)-v(j,i-1))/dx2(j)
              - (u(j+1,i)-u(j-1,i))/dy2(j)
              + (u(j,i)/a)*tan(phi(j)) #since meridians aren't parallel

  The last terms accounts for the convergence of the meridians on a sphere. 
  '''
  
  #we'll do the above centered finite differencing for the non-pole latitudes.
  #for the poles, we'll do a finite volume \int gradxu dA = \int u.n dS since 
  #trying to finite difference it confuses me. remember that the poles are really
  #nLon copies of the same point
  
  vort = np.empty((nLat, nLon), dtype=float)
  
  dRadLat = np.pi/(nLat-1) #[-pi/2, pi/2], ie with values at both poles
  dRadLon = 2.*np.pi/nLon #[0,2pi)
  dy = r*dRadLat; dy2 = 2.*dy #arc length on a sphere
  
  #calc values for non poles
  for iLat in xrange(1,nLat-1):
    tanphi = np.tan(lat[iLat])/r
    dx = r*np.cos(lat[iLat])*dRadLon; dx2 = 2.*dx
    for iLon in xrange(nLon):
      iWest = (iLon-1)%nLon # -1%4=3 so don't worry about negatives
      iEast = (iLon+1)%nLon
      iSouth = iLat+1; iNorth = iLat-1
      dv_dx = (v[iLat, iEast]-v[iLat, iWest])/dx2
      du_dy = (u[iNorth, iLon]-u[iSouth, iLon])/dy2
      meridTerm = u[iLat, iLon]*tanphi
      vort[iLat, iLon] = dv_dx-du_dy+meridTerm
      
  #calc values for north and south poles with finite volume approach
  #around next latitude equatorward of pole
  iLat = nLat-1;
  dx = r*np.cos(lat[iLat-1])*dRadLon #for evenly spaced lats, same dx for south and north
  #for area of the cap, http://mathworld.wolfram.com/SphericalCap.html
  a = dx/dRadLon
  h = r-np.sqrt(r*r-a*a)
  areaCap = 2.*np.pi*r*h
  #around south pole, remember integrate with domain on left
  undS = np.sum(u[iLat-1,:])*-dx #since +dx has domain on right
  vort[iLat,:] = undS/areaCap
  iLat = 0
  undS = np.sum(u[iLat+1,:])*dx
  vort[iLat,:] = undS/areaCap
  
  return vort

def calc_potentialTemperature(tmp, press):
    
  Cp = 1004.5; Rd = 287.04;
  Rd_cp = Rd/Cp; p0 = 1.e5
  theta = tmp*((p0/press)**Rd_cp)
  
  return theta

def eraTimeToCalendarTime(hrs):
  #in the ERA file, time is stored as "hours since 1900-01-01 00:00:0.0"
  #we'll convert that to a datetime object and return a nice looking string
  
  tBase = dt.datetime(1900, 1, 1, 0)
  #note that TypeError: unsupported type for timedelta hours component: numpy.int32
  tNew = tBase + dt.timedelta(hours=hrs)
  tTuple = dt.datetime.timetuple(tNew);
  s = time.strftime('%Y-%m-%d_%H', tTuple)
  return s  

def demo_eraI(fMesh, filesDataIn, fNameOut, r, dRegion, latThresh, iTimeStart_fData, iTimeEnd_fData, info='eraI case'):
  #mesh ---------------------
  data = netCDF4.Dataset(fMesh,'r')
  d2r = np.pi/180.; 
  lat = data.variables['latitude'][:]*d2r; lon = data.variables['longitude'][:]*d2r
  #want latitudes to be in [-pi/2, pi/2] and longitudes in [0, 2pi)
  lon = lon%(2.*np.pi)
  data.close()
  
  mesh = llMesh.Mesh(lat,lon, r, dRegion)
  #mesh.fill_latCellArea()
  mesh.fill_inDisk()
  mesh.fill_inRegion(latThresh)
  cell0 = llMesh.Cell(mesh,0)
  
  #metr fields -----------------
  nFiles = len(filesDataIn)
  if (nFiles<1):
    return mesh, cell0
  
  dataOut = write_netcdf_header_metr(fNameOut, info, mesh.nCells)
  iTimeGlobal = 0
  for iFile in xrange(nFiles):
    fPath = filesDataIn[iFile]
    data = netCDF4.Dataset(fPath,'r')
    
    #loop over individual times ------------------------------
    #times = data.variables['time'][:]; nTimes = len(times); nTimes = 20
    #for iTime in xrange(nTimes):
    iTimeStart = iTimeStart_fData[iFile]; iTimeEnd = iTimeEnd_fData[iFile]
    if (iTimeEnd<0): #use all times in file
      times = data.variables['time'][:]; nTimes = len(times);
      iTimeEnd = nTimes-1
    for iTime in xrange(iTimeStart,iTimeEnd+1):
      #read from file
      theta = data.variables['pt'][iTime,:,:]
      u = data.variables['u'][iTime,:,:]; v = data.variables['v'][iTime,:,:]
      
      #fill in missing values w/in region
      #ERA-I doesn't appear to have any missing values...I don't know how their interpolation works.
      #Some old documentation described PP2DINT that extrapolates using constant values.
      #This rando site: https://badc.nerc.ac.uk/data/ecmwf-op/levels.html
      #says 
      #"The ECMWF Operational and ERA-40 datasets also provide data on a "PV=+/-2" surface on which the potential vorticity takes the value 2PVU in the northern hemisphere and -2PVU in the southern hemisphere (1PVU = 10-6 m2 s-1 K kg-1), provided such a surface can be found searching downwards from the Model level close to 96hPa. Values at this model level are used where the search is unsuccessful."
      
      #compute additional fields
      vort = calc_vertVorticity_ll(u, v, mesh.nLat, mesh.nLon, mesh.lat, r)
    
      #write to file
      u = helpers.flatten_2dTo1d(u, mesh.nLat, mesh.nLon)
      v = helpers.flatten_2dTo1d(v, mesh.nLat, mesh.nLon)
      theta = helpers.flatten_2dTo1d(theta, mesh.nLat, mesh.nLon)
      vort = helpers.flatten_2dTo1d(vort, mesh.nLat, mesh.nLon)
      
      write_netcdf_iTime_metr(dataOut, iTimeGlobal, u,v,theta,vort)
      iTimeGlobal = iTimeGlobal+1
    #end iTime
  #end iFile
  dataOut.close()
  
  return mesh, cell0

def demo_mpas(fMesh, filesDataIn, fNameOut, r, dRegion, latThresh, iTimeStart_fData, iTimeEnd_fData, info='mpas case'):
  #mesh ---------------------
  data = netCDF4.Dataset(fMesh,'r')
  lat = data.variables['latCell'][:]; lon = data.variables['lonCell'][:]
  lon = lon%(2.*np.pi) #want latitudes to be in [-pi/2, pi/2] and longitudes in [0, 2pi)
  
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  areaCell = data.variables['areaCell'][:]
  data.close()
  
  mesh = mpasMesh.Mesh(lat,lon, areaCell, cellsOnCell, nEdgesOnCell, r, dRegion)
  mesh.fill_inRegion(latThresh)
  cell0 = mpasMesh.Cell(mesh,0)
  
  #metr fields -----------------
  nFiles = len(filesDataIn)
  if (nFiles<1):
    return mesh, cell0
  
  dataOut = write_netcdf_header_metr(fNameOut, info, mesh.nCells)
  iTimeGlobal = 0
  for iFile in xrange(nFiles):
    fPath = filesDataIn[iFile]
    data = netCDF4.Dataset(fPath,'r')
    
    #loop over individual times ------------------------------
    #times = data.variables['time'][:]; nTimes = len(times); nTimes = 20
    #for iTime in xrange(nTimes):
    iTimeStart = iTimeStart_fData[iFile]; iTimeEnd = iTimeEnd_fData[iFile]
    if (iTimeEnd<0): #use all times in file
      nTimes = len(data.dimensions['Time'])
      iTimeEnd = nTimes-1
    for iTime in xrange(iTimeStart,iTimeEnd+1):
      #read from file
      theta = data.variables['theta_pv'][iTime,:]
      u = data.variables['u_pv'][iTime,:]; v = data.variables['v_pv'][iTime,:]
      vort = data.variables['vort_pv'][iTime,:]
      
      #fill in missing values w/in region
      #MPAS will have surface values if whole column is above 2pvu.
      #Unclear what to do if whole column is below 2pvu (like near equator)
      
      #compute additional fields
    
      #write to file      
      write_netcdf_iTime_metr(dataOut, iTimeGlobal, u,v,theta,vort)
      iTimeGlobal = iTimeGlobal+1
    #end iTime
  #end iFile
  dataOut.close()
  
  return mesh, cell0
  
def write_netcdf_header_metr(fName, info, nCells):
  
  data = netCDF4.Dataset(fName, 'w', format='NETCDF4')
  data.description = info
  
  # dimensions
  data.createDimension('time', None)
  data.createDimension('nCells', nCells)

  # variables
  u_data = data.createVariable('u', 'f8', ('time','nCells',))
  v_data = data.createVariable('v', 'f8', ('time','nCells',))
  theta_data = data.createVariable('theta', 'f8', ('time','nCells',))
  vort_data = data.createVariable('vort', 'f8', ('time','nCells',))
  
  #units and descriptions
  u_data.units = 'm/s'
  v_data.units = 'm/s'
  theta_data.units = 'K'
  vort_data.units = '1/s'
  
  return data
  
def write_netcdf_iTime_metr(data, iTime, u,v,theta,vort):
  # fill file. with time as unlimited, dimension will just keep growing
  
  data.variables['u'][iTime,:] = u[:]
  data.variables['v'][iTime,:] = v[:]
  data.variables['theta'][iTime,:] = theta[:]
  data.variables['vort'][iTime,:] = vort[:]
#

