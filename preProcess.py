'''
Take input data and write mesh+metr info to our format.
mesh: lat,lon,areaCell
metr: u,v,theta,verticalVorticity,inRegion on DT
'''

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

import helpers
import llMesh
#import mpasMesh
#import wrfUniformMesh as wrfMesh

def get_missingCells_file(inVar):
  #if search vertical column down from top, there can be no 2pvu value if:
  #-entire column is above 2pvu: DT below sfc
  #-entire column is below 2pvu: wrong hemisphere or low pv column (tropics, anticyclonic,...)
  
  isMissing = inVar.mask
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
   
def calc_vorticity_wrfTrop_uniform(u, v, dx, dy, mapFac=1.0):
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
  lat = data.variables['lat'][:]
  lat = lat*d2r

  lon = data.variables['lon'][:] 
  lon = lon*d2r

  #want latitudes to be in [-pi/2, pi/2] and longitudes in [0, 2pi)
  lon = lon%(2.*np.pi)  #% is a modulus operator used to find the remainder of lon divided by 2*pi

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
  
  data_g = netCDF4.Dataset(filesDataIn[0],'r')
  data_u = netCDF4.Dataset(filesDataIn[1],'r')
  data_v = netCDF4.Dataset(filesDataIn[2],'r') 
    
  levs = np.array(data_u.variables['lev'])
  lev_1000 = np.where(levs==1000)[0][0]     
  lev_500 = np.where(levs==500)[0][0] 
  lev_700 = np.where(levs==700)[0][0]     
  
  dataOut = write_netcdf_header_metr(fNameOut, info, mesh)
      
  iTimeGlobal = 0

  iTimeStart = iTimeStart_fData; iTimeEnd = iTimeEnd_fData
  
  if (iTimeEnd<0): #use all times in file
    times = data_g.variables['time'][:]; nTimes = len(times);
    iTimeEnd = nTimes-1  
    
  for iTime in xrange(iTimeStart,iTimeEnd+1): 
  
    g1000 = data_g.variables['g'][iTime,lev_1000,:,:] # 1000-hPa geopotential height
    g500 = data_g.variables['g'][iTime,lev_500,:,:] # 500-hPa geopotential height

    theta = g500-g1000 # This is really 1000-500-hPa thickness, but in order to reduce the 
                       # need to alter the original TPV code, I will just refer to thickness as theta
    
    u1000 = data_u.variables['u'][iTime,lev_1000,:,:] # 1000-hPa u-wind
    v1000 = data_v.variables['v'][iTime,lev_1000,:,:] # 1000-hPa v-wind

    u500 = data_u.variables['u'][iTime,lev_500,:,:] # 500-hPa u-wind
    v500 = data_v.variables['v'][iTime,lev_500,:,:] # 500-hPa v-wind

    u700 = data_u.variables['u'][iTime,lev_700,:,:] # 700-hPa u-wind
    v700 = data_v.variables['v'][iTime,lev_700,:,:] # 700-hPa v-wind

    vort500 = calc_vertVorticity_ll(u500, v500, mesh.nLat, mesh.nLon, mesh.lat, r) # 500-hPa relative vorticity
    vort1000 = calc_vertVorticity_ll(u1000, v1000, mesh.nLat, mesh.nLon, mesh.lat, r) #1000-hPa relative vorticity
    
    vort = vort500-vort1000  # 500-1000-hPa thermal vorticity

    #write to file
    u700 = helpers.flatten_2dTo1d(u700, mesh.nLat, mesh.nLon)
    v700 = helpers.flatten_2dTo1d(v700, mesh.nLat, mesh.nLon)
    theta = helpers.flatten_2dTo1d(theta, mesh.nLat, mesh.nLon)
    vort = helpers.flatten_2dTo1d(vort, mesh.nLat, mesh.nLon)
        
    write_netcdf_iTime_metr(dataOut,iTimeGlobal,u700,v700,theta,vort)
    iTimeGlobal = iTimeGlobal+1
    print(iTime)
    print("another time done preprocessing...")
    #end iTime
  dataOut.close()
  
  return mesh, cell0    



def write_netcdf_header_metr(fName, info, mesh):
  
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
  latCell_data.units = 'radians'
  lonCell_data.units = 'radians'
  u_data.units = 'm/s'
  v_data.units = 'm/s'
  theta_data.units = 'dam'
  vort_data.units = '1/s'
  #ws_data.units = 'm/s'

  
  
  #fill lat/lon
  allCells = np.arange(nCells)
  lat, lon = mesh.get_latLon_inds(allCells)
  data.variables['latCell'][:] = lat[:]
  data.variables['lonCell'][:] = lon[:]
  
  return data
  
def write_netcdf_iTime_metr(data, iTime,u,v,theta,vort):
  # fill file. with time as unlimited, dimension will just keep growing
  
  data.variables['u'][iTime,:] = u[:]
  data.variables['v'][iTime,:] = v[:]
  data.variables['theta'][iTime,:] = theta[:]
  data.variables['vort'][iTime,:] = vort[:]
  
#

def plot_metr(fMetr):
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
      
