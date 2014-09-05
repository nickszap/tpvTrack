
import netCDF4
import matplotlib.pyplot as plt

import my_settings
import preProcess
import llMesh
import segment

def demo():
  #setup -----------------
  info = my_settings.info
  filesData = my_settings.filesData
  fMesh = filesData[0]  
  fMetr = my_settings.fDirSave+'fields_debug.nc'
  fSeg = my_settings.fDirSave+'seg_debug.nc'
  
  rEarth = my_settings.rEarth
  dRegion = my_settings.dFilter
  latThresh = my_settings.latThresh
  
  #pre-process ------------------------
  #mesh = preProcess.demo_eraI(fMesh, filesData, fMetr, my_settings.rEarth, dRegion, latThresh, info=info)
  mesh = preProcess.demo_eraI(fMesh, [], fMetr, my_settings.rEarth, dRegion, latThresh) #if already processed input data
  
  if (False):
    print mesh.lat; print mesh.nLat
    print mesh.lon; print mesh.nLon
  
  #segment --------------------------
  cell0 = llMesh.Cell(mesh,0)
  if (True):
    print 'index: ', cell0.ind, 'nbrs: ', cell0.get_nbrInds()
  
  dataMetr = netCDF4.Dataset(fMetr,'r')
  #segment.run_segment(fSeg, info, dataMetr, cell0, mesh)
  
  segment.run_plotBasins(my_settings.fDirSave, dataMetr, fSeg, mesh)
  dataMetr.close()
  
  #spatial metrics ------------------------
  
  #time correspondence -----------------
  
  #time tracks -------------------------
  
  #time metrics ----------------------

if __name__=='__main__':
  demo()
