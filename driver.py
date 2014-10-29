import netCDF4
import matplotlib.pyplot as plt

import my_settings
import preProcess
import llMesh
import segment
import basinMetrics
import correspond
import tracks

def demo():
  #setup -----------------
  info = my_settings.info
  filesData = my_settings.filesData
  fMesh = my_settings.fMesh  
  fMetr = my_settings.fMetr
  fSeg = my_settings.fSeg
  fCorr = my_settings.fCorr
  fTrack = my_settings.fTrack
  fMetrics = my_settings.fMetrics
  
  rEarth = my_settings.rEarth
  dRegion = my_settings.dFilter
  latThresh = my_settings.latThresh
  
  #pre-process ------------------------
  if (my_settings.doPreProc):
    mesh = preProcess.demo_eraI(fMesh, filesData, fMetr, 
                                my_settings.rEarth, dRegion, latThresh, 
                                my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData, info=info)
  else:
    #if already processed input data
    mesh = preProcess.demo_eraI(fMesh, [], fMetr, my_settings.rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData)
  
  cell0 = llMesh.Cell(mesh,0)
  if (True):
    print 'index: ', cell0.ind, 'nbrs: ', cell0.get_nbrInds()
  
  #segment --------------------------
  dataMetr = netCDF4.Dataset(fMetr,'r'); nTimes = len(dataMetr.dimensions['time'])
  if (my_settings.doSeg):
    segment.run_segment(fSeg, info, dataMetr, cell0, mesh)
    segment.run_plotBasins(my_settings.fDirSave, dataMetr, fSeg, mesh)
  
  dataMetr.close()
  
  #spatial metrics ------------------------
  dataMetr = netCDF4.Dataset(fMetr,'r')
  dataSeg = netCDF4.Dataset(fSeg,'r')
  if (my_settings.doMetrics):
    basinMetrics.run_metrics(fMetrics, info, mesh, dataMetr, dataSeg, 0, nTimes-1)
  
  dataMetr.close()
  dataSeg.close()
  
  #basinMetrics.print_metrics(fMetrics)
  
  #time correspondence -----------------
  dataMetr = netCDF4.Dataset(fMetr,'r')
  dataSeg = netCDF4.Dataset(fSeg,'r')
  dataMetrics = netCDF4.Dataset(fMetrics, 'r')
  if (my_settings.doCorr):
    correspond.run_correspond(fCorr, dataMetr, dataSeg, mesh, my_settings.deltaT, my_settings.trackMinMaxBoth, my_settings.areaOverlap, 0, nTimes-1, dataMetrics)
    correspond.plot_correspondences(my_settings.fDirSave, fCorr, nTimes-1, mesh)
  
  dataMetrics.close()
  dataSeg.close()
  dataMetr.close()
  
  #time tracks -------------------------
  if (my_settings.doTracks):
    #since appending to fTrack over time, wipe file before starting (if it exists)
    my_settings.silentremove(fTrack)
    
    tracks.run_tracks_timeInterval(fTrack, fCorr, 0, nTimes-2, fMetrics=fMetrics, trackOnlyMajor=True)
    #tracks.plot_tracks_cells(fTrack, mesh, my_settings.fDirSave)
  
  #time metrics ----------------------

if __name__=='__main__':
  demo()
