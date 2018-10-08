'''
Take input data and write mesh+metr info to our format.
mesh: lat,lon,areaCell
metr: u,v,theta,verticalVorticity,inRegion on DT
'''

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import helpers
import llMesh
import mpasMesh
import wrfUniformMesh as wrfMesh

def get_missingCells_file(data):
  """Return mask of missing values, coded for GFS netCDF4 object"""
  #if search vertical column down from top, there can be no 2pvu value if:
  #-entire column is above 2pvu: DT below sfc
  #-entire column is below 2pvu: wrong hemisphere or low pv column (tropics, anticyclonic,...)
  
  pvInd = 1
  isMissing = data.variables['TMP_P0_L109_GLL0'][pvInd,:,:].mask
  return isMissing
  
def fill_missingVals_region(valsIn, nLat, nLon, isMissing, inRegion):
  """For missing values in region that is used, fill value is average of non-missing neighbors"""
  
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
  """Return data of variables needed from file, with no time index, in SI units. Coded for GFS"""
  
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
  Calculate vertical vorticity on a latitude/longitude mesh
  
  Arguments:
  u - zonal wind
  v - meridional wind
  nLat - number of latitude points
  nLon - number of longitude points
  lat - latitudes
  r - radius of sphere
  
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
   
def calc_vorticity_wrfTrop_uniform(u, v, dx, dy, mapFac=1.0):
  """Calculate vertical vorticity on a uniformly spaced WRF domain"""
  #Steven's files have variables already processed to cell centers.
  #u,v come in ordered [south_north=y, west_east=x]
  #we'll use numpy.gradient for the finite difference,
  #e.g., http://stackoverflow.com/questions/17901363/gradient-calculation-with-python
  
  du_dy, du_dx = np.gradient(u, dy, dx)
  dv_dy, dv_dx = np.gradient(v, dy, dx)
  
  #In physical space, dx and dy are not constants.
  #If we act like we were finite differencing across the faces of each cell,
  #grid stretching gets lumped into O(approx) and we just scale the d/dDirection
  dv_dx *= mapFac # d/dxEarth = d/(dxGrid/mapFac)
  du_dy *= mapFac
  
  return dv_dx-du_dy

def calc_potentialTemperature(tmp, press):
  """
  Return potential temperature
  
  Arguments:
  tmp - temperature (K)
  press - pressure (Pa)
  """
  Cp = 1004.5; Rd = 287.04;
  Rd_cp = Rd/Cp; p0 = 1.e5
  theta = tmp*((p0/press)**Rd_cp)
  
  return theta

def eraTimeToCalendarTime(hrs):
  """
  In the ERA file, time is stored as "hours since 1900-01-01 00:00:0.0"
  we'll convert that to a datetime object and return a nice looking string
  """
  
  tBase = dt.datetime(1900, 1, 1, 0)
  #note that TypeError: unsupported type for timedelta hours component: numpy.int32
  tNew = tBase + dt.timedelta(hours=hrs)
  tTuple = dt.datetime.timetuple(tNew);
  s = time.strftime('%Y-%m-%d_%H', tTuple)
  return s  

def demo_eraI(fMesh, filesDataIn, fNameOut, r, dRegion, latThresh, iTimeStart_fData, iTimeEnd_fData, info='eraI case'):
  """
  Pre-process ERA-Interim data into tpvTrack format
  
  Arguments:
  fMesh - path to file with mesh inforamtion
  filesDataIn - Input ERA-I filepaths
  fNameOut - filepath for output file
  r - radius of sphere
  dRegion - radius of neighborhood
  latThresh - latitude cutoff for subset of domain used for segmentation, tracking,...
  iTimeStart_fData - integer index for start time of each file
  iTimeEnd_fData - integer index for end time of each file (can use -1 for last time in file)
  """
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
  cell0 = llMesh.Cell(mesh,-1)
  
  #metr fields -----------------
  nFiles = len(filesDataIn)
  if (nFiles<1):
    return mesh, cell0
  
  dataOut = write_netcdf_header_metr(fNameOut, info, mesh)
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
  """
  Pre-process MPAS data into tpvTrack format
  
  Arguments:
  fMesh - path to file with mesh inforamtion
  filesDataIn - Input MPAS output filepaths
  fNameOut - filepath for output file
  r - radius of sphere
  dRegion - radius of neighborhood
  latThresh - latitude cutoff for subset of domain used for segmentation, tracking,...
  iTimeStart_fData - integer index for start time of each file
  iTimeEnd_fData - integer index for end time of each file (can use -1 for last time in file)
  """
  
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
  cell0 = mpasMesh.Cell(mesh,-1)
  
  #metr fields -----------------
  nFiles = len(filesDataIn)
  if (nFiles<1):
    return mesh, cell0
  
  dataOut = write_netcdf_header_metr(fNameOut, info, mesh)
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

def demo_wrf_trop(fMesh, filesDataIn, fNameOut, r, dRegion, latThresh, iTimeStart_fData, iTimeEnd_fData, fMapProj, info='wrf case', pvIndex=3):
  """
  Pre-process WRF diagnosed tropopause data into tpvTrack format
  
  Arguments:
  fMesh - path to file with mesh inforamtion
  filesDataIn - Input ERA-I filepaths
  fNameOut - filepath for output file
  r - radius of sphere
  dRegion - radius of neighborhood
  latThresh - latitude cutoff for subset of domain used for segmentation, tracking,...
  iTimeStart_fData - integer index for start time of each file
  iTimeEnd_fData - integer index for end time of each file (can use -1 for last time in file)
  fMapProj - filepath with information about domain's map projection
  """
  
  #For Steven's wrfout_trop files that have already been processed in a particular way.
  #I think the grid is oriented such that u,v are both grid and global zonal,meridional velocities.
  #If this isn't true, there's some figuring out to do.
  
  #mesh ---------------------
  data = netCDF4.Dataset(fMesh,'r')
  dataProj = netCDF4.Dataset(fMapProj,'r')
  d2r = np.pi/180.; 
  lat = data.variables['XLAT'][0,:,:]*d2r; lon = data.variables['XLONG'][0,:,:]*d2r
  dx = data.DX
  dy = data.DY
  #want latitudes to be in [-pi/2, pi/2] and longitudes in [0, 2pi)
  lon = lon%(2.*np.pi)
  data.close()
  
  mesh = wrfMesh.Mesh(lat,lon, dx, dy, r, dRegion)
  mesh.fill_inRegion(latThresh)
  
  mapFac = dataProj.variables['MAPFAC_M'][0,:,:]
  #mapFac = 1.0
  mesh.fill_cellArea(mapFac)
  cell0 = wrfMesh.Cell(mesh,-1)
  
  #metr fields -----------------
  nFiles = len(filesDataIn)
  if (nFiles<1):
    return mesh, cell0
  
  cosalpha = dataProj.variables['COSALPHA'][0,:,:]
  sinalpha = dataProj.variables['SINALPHA'][0,:,:]
  dataProj.close()
  
  dataOut = write_netcdf_header_metr(fNameOut, info, mesh)
  iTimeGlobal = 0
  for iFile in xrange(nFiles):
    fPath = filesDataIn[iFile]
    data = netCDF4.Dataset(fPath,'r')
    
    #loop over individual times ------------------------------
    #times = data.variables['time'][:]; nTimes = len(times); nTimes = 20
    #for iTime in xrange(nTimes):
    iTimeStart = iTimeStart_fData[iFile]; iTimeEnd = iTimeEnd_fData[iFile]
    if (iTimeEnd<0): #use all times in file
      nTimes = len(data.dimensions['time']);
      iTimeEnd = nTimes-1
    for iTime in xrange(iTimeStart,iTimeEnd+1):
      #read from file
      theta = data.variables['THETA'][iTime,pvIndex,:,:]
      u = data.variables['U'][iTime,pvIndex,:,:]
      v = data.variables['V'][iTime,pvIndex,:,:]
      
      #fill in missing values w/in region
      #apparently netCDF4 doesn't return a masked array if all masks are False
      try:
        print "Time {0}, number of missing values: {1}".format(iTimeGlobal, np.sum(theta.mask))
        for var in (theta, u, v):
          #we'll just replace missing values with global mean
          meanVal = np.mean(var); print "Replacing values with: ",meanVal
          var[var.mask==True] = meanVal
        #print "Missing values still in theta? ", True in theta.mask
        theta = theta.data;
        u = u.data
        v = v.data
      except AttributeError:
        print "Time {0}, number of missing values: {1}".format(iTimeGlobal,0)
      
      #compute additional fields
      vort = calc_vorticity_wrfTrop_uniform(u, v, dx, dy, mapFac=mapFac)
      #rotate grid-relative wind to global (apparently the stored rotation is for earth->grid)
      uGlobal = u*cosalpha-v*sinalpha
      vGlobal = u*sinalpha+v*cosalpha
      u = uGlobal; v = vGlobal
      
      #write to file
      u = helpers.flatten_2dTo1d(u, mesh.ny, mesh.nx)
      v = helpers.flatten_2dTo1d(v, mesh.ny, mesh.nx)
      theta = helpers.flatten_2dTo1d(theta, mesh.ny, mesh.nx)
      vort = helpers.flatten_2dTo1d(vort, mesh.ny, mesh.nx)
      
      write_netcdf_iTime_metr(dataOut, iTimeGlobal, u,v,theta,vort)
      iTimeGlobal = iTimeGlobal+1
    #end iTime
  #end iFile
  dataOut.close()
  
  return mesh, cell0  

def write_netcdf_header_metr(fName, info, mesh):
  """Create file and write header for tpvTrack preprocess netcdf file"""
  data = netCDF4.Dataset(fName, 'w', format='NETCDF4')
  data.description = info
  
  # dimensions
  nCells = mesh.nCells
  data.createDimension('time', None)
  data.createDimension('nCells', nCells)

  # variables
  latCell_data = data.createVariable('latCell', 'f8', ('nCells',))
  lonCell_data = data.createVariable('lonCell', 'f8', ('nCells',))
  u_data = data.createVariable('u', 'f8', ('time','nCells',))
  v_data = data.createVariable('v', 'f8', ('time','nCells',))
  theta_data = data.createVariable('theta', 'f8', ('time','nCells',))
  vort_data = data.createVariable('vort', 'f8', ('time','nCells',))
  
  #units and descriptions
  latCell_data.units = 'radians'; latCell_data.long_name='Latitude'
  lonCell_data.units = 'radians'; lonCell_data.long_name='Longitude'
  u_data.units = 'm/s'; u_data.long_name='Zonal velocity'
  v_data.units = 'm/s'; v_data.long_name='Meridional velocity'
  theta_data.units = 'K'; theta_data.long_name='Potential temperature'
  vort_data.units = '1/s'; vort_data.long_name='Vertical vorticity'
  
  #fill lat/lon
  allCells = np.arange(nCells)
  lat, lon = mesh.get_latLon_inds(allCells)
  data.variables['latCell'][:] = lat[:]
  data.variables['lonCell'][:] = lon[:]
  
  return data
  
def write_netcdf_iTime_metr(data, iTime, u,v,theta,vort):
  """Write one time into tpvTrack preprocess netcdf file"""
  # fill file. with time as unlimited, dimension will just keep growing
  
  data.variables['u'][iTime,:] = u[:]
  data.variables['v'][iTime,:] = v[:]
  data.variables['theta'][iTime,:] = theta[:]
  data.variables['vort'][iTime,:] = vort[:]
#

def plot_metr(fMetr):
  """Example of plotting preprocess variables on map"""
  data = netCDF4.Dataset(fMetr,'r')
  nTimes = len(data.dimensions['time'])
  lat = data.variables['latCell'][:]*180./np.pi
  lon = data.variables['lonCell'][:]*180./np.pi
  
  m = Basemap(projection='ortho',lon_0=0,lat_0=90., resolution='l')
  x,y = m(lon, lat)
  
  keys = ['u','v','theta','vort']
  bounds = [[-50.,50.], [-50.,50.], [270.,360.], [-1.e-4,1.e-4]]
  for iTime in xrange(nTimes):
    for iKey in xrange(len(keys)):
      key = keys[iKey]
      keyRange = bounds[iKey]; keyMin = keyRange[0]; keyMax = keyRange[1]
      plt.figure()
      vals = data.variables[key][iTime,:]
      m.drawcoastlines()
      m.pcolor(x,y,vals,tri=True, shading='flat',edgecolors='none',cmap=plt.cm.RdBu_r, vmin=keyMin, vmax=keyMax)
      plt.colorbar()
      plt.title(key)
      
      s = '{0}_t{1}.png'.format(key, iTime)
      plt.savefig(s, bbox_inches='tight'); plt.close()

if __name__ == '__main__':
  fMetr = '/data01/tracks/wrf/algo/fields_debug.nc'
  plot_metr(fMetr)
      
# ------------------------- Untested code --------------------------------

def calc_vertVorticity_wrf(U,V,MSFU,MSFV,MSFT,DX,DY,NX,NY):
  print "Uhoh. calc_vertVorticity_wrf function is untested!!!"
  '''
  Adapted from DCOMPUTEABSVORT(AV,U,V,MSFU,MSFV,MSFT,COR,DX,DY,NX,NY,
     +                           NZ,NXP1,NYP1)
  from NCL source code (https://github.com/yyr/ncl/blob/master/ni/src/lib/nfpfort/wrf_pvo.f)
  
  u,v: unstaggered grid-relative winds
  msf{u,v,t}: map-scale factors
  d{x,y}: computational grid spacing
  N{x,y}: # of cells in each direction
  2D arrays are indexed on the grid (S_N=v, W_E=u)
  
  For a uniform grid, could just do:
  du_dy, du_dx = np.gradient(u, dy, dx)
  dv_dy, dv_dx = np.gradient(v, dy, dx)
  return dv_dx-du_dy
  '''
  
  #for interior cells, do finite difference between neighboring cells, e.g., (valNorth-valSouth)/2dy.
  #to get u at midpt of horizontal cell, average left and right boundaries.
  vort = np.empty((NX,NY),dtype=float)
  for J in xrange(NY):
    JP1 = min(J+1,NY-1)
    JM1 = max(J-1,0)
    for I in xrange(NX):
      IP1 = min(I+1,NX-1)
      IM1 = max(I-1,0)
      
      DSX = (IP1-IM1)*DX
      DSY = (JP1-JM1)*DY
      MM = MSFT[j,i]*MSFT[j,i]
      
      du_dy = .5* (U[JP1,I]/MSFU[JP1,I]+ U[JP1,I+1]/MSFU[JP1,I+1] -
                   U[JM1,I]/MSFU[JM1,I]- U[JM1,I+1]/MSFU[JM1,I+1])/(DSY/MM)
      #
      dv_dx = .5*(V[J,IP1]/MSFV[J,IP1]+ V[J+1,IP1]/MSFV[J+1,IP1] -
                  V[J,IM1]/MSFV[J,IM1]+ V[J+1,IM1]/MSFV[J+1,IM1])/(DSX/MM)
      #
      vort[I,J] = dv_dx-du_dy
      
  return vort
  
def demo_wrfUntested(fMesh, filesDataIn, fNameOut, r, dRegion, latThresh, iTimeStart_fData, iTimeEnd_fData, info='wrf case'):
  #mesh ---------------------
  data = netCDF4.Dataset(fMesh,'r')
  d2r = np.pi/180.; 
  lat = data.variables['XLAT'][0,:,:]*d2r; lon = data.variables['XLONG'][0,:,:]*d2r
  dx = data.variables['RDX'][0]; dx = 1./dx;
  dy = data.variables['RDY'][0]; dy = 1./dy;
  #want latitudes to be in [-pi/2, pi/2] and longitudes in [0, 2pi)
  lon = lon%(2.*np.pi)
  data.close()
  
  mesh = wrfMesh.Mesh(lat,lon, r, dRegion)
  mesh.fill_inRegion(latThresh)
  cell0 = wrfMesh.Cell(mesh,-1)
  
  #metr fields -----------------
  nFiles = len(filesDataIn)
  if (nFiles<1):
    return mesh, cell0
  
  dataOut = write_netcdf_header_metr(fNameOut, info, mesh)
  iTimeGlobal = 0
  for iFile in xrange(nFiles):
    fPath = filesDataIn[iFile]
    data = netCDF4.Dataset(fPath,'r')
    
    #loop over individual times ------------------------------
    #times = data.variables['time'][:]; nTimes = len(times); nTimes = 20
    #for iTime in xrange(nTimes):
    iTimeStart = iTimeStart_fData[iFile]; iTimeEnd = iTimeEnd_fData[iFile]
    if (iTimeEnd<0): #use all times in file
      nTimes = len(data.dimensions['Time']);
      iTimeEnd = nTimes-1
    for iTime in xrange(iTimeStart,iTimeEnd+1):
      #read from file
      theta = data.variables['pt'][iTime,:,:]
      ug = data.variables['U'][iTime,:,:]; #south-north, west-east-stagger
      vg = data.variables['V'][iTime,:,:]; #south-north-stagger, west-east
      
      mapfac = data.variables['MAPFAC_M'][iTime,:,:]
      sinalpha = data.variables['SINALPHA'][iTime,:,:]
      cosalpha = data.variables['COSALPHA'][iTime,:,:]
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

