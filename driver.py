import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import datetime as dt

import my_settings
import preProcess
import llMesh
import segment
import basinMetrics
import correspond
import tracks

printTiming = True
from datetime import datetime

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
  
  myRank = my_settings.myRank
  nRanks = my_settings.nRanks
  commWorld = my_settings.commWorld
  
  #datetime objects
  timeStartGlobal = my_settings.timeStart
  deltaTGlobal = my_settings.timeDelta
  
  #pre-process ------------------------
  if (my_settings.doPreProc):
    if (printTiming):
      tStart = datetime.now()
      
    if (my_settings.inputType=='eraI'):
      mesh, cell0 = preProcess.demo_eraI(fMesh, filesData, fMetr, rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData, info=info)
    elif (my_settings.inputType=='mpas'):
      mesh, cell0 = preProcess.demo_mpas(fMesh, filesData, fMetr, rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData, info=info)
    elif (my_settings.inputType=='wrf_trop'):
      mesh, cell0 = preProcess.demo_wrf_trop(fMesh, filesData, fMetr, rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData, 
                                             my_settings.fileMap, info=info, pvIndex=3)
    else:
      print "Unrecognized input type in my_settings: ",my_settings.inputType
      
    if (printTiming):
      tEnd = datetime.now()
      print "Time doPreProc: ", tEnd-tStart
  else:
    #if already processed input data
    if (my_settings.inputType=='eraI'):
      mesh, cell0 = preProcess.demo_eraI(fMesh, [], fMetr, rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData)
    elif (my_settings.inputType=='mpas'):
      mesh, cell0 = preProcess.demo_mpas(fMesh, [], fMetr, rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData)
    elif (my_settings.inputType=='wrf_trop'):
      mesh, cell0 = preProcess.demo_wrf_trop(fMesh, [], fMetr, rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData, my_settings.fileMap)
    else:
      print "Unrecognized input type in my_settings: ",my_settings.inputType
  
  #segment --------------------------
  dataMetr = netCDF4.Dataset(fMetr,'r'); 
  nTimes = len(dataMetr.dimensions['time']); #nTimes = 5
  iStartGlobal=0; iEndGlobal = nTimes-1
  if (my_settings.doSeg):
    if (printTiming):
      tStart = datetime.now()
    
    #each worker segments over time independently
    iStartRank,iEndRank = my_settings.getLimits_startStop(iStartGlobal, iEndGlobal, iWork=myRank, nWork=nRanks)  
    #segment.run_segment(fSeg, info, dataMetr, cell0.copy(), mesh, iStartRank, iEndRank)
    
    #combine each workers files into 1 consistent one
    commWorld.Barrier() #all workers need to finish their seg to combine into 1 file
    if (myRank == 0):
      iTimesStart_ranks = []; iTimesEnd_ranks = []; filesIn_ranks = []
      for iRank in xrange(nRanks):
        iStartRank,iEndRank = my_settings.getLimits_startStop(iStartGlobal, iEndGlobal, iWork=iRank, nWork=nRanks)
        iTimesStart_ranks.append(iStartRank); iTimesEnd_ranks.append(iEndRank)
        filesIn_ranks.append( my_settings.fSegFmt.format(iRank) )
      #
      print 'Combining parallel seg files: ', myRank, iTimesStart_ranks, iTimesEnd_ranks, filesIn_ranks  
      segment.combineParallelFiles(my_settings.fSegFinal, iTimesStart_ranks, iTimesEnd_ranks, filesIn_ranks)
    commWorld.Barrier() #I/O process needs to finish with serial file before we move on
    
    if (False):
      segment.run_plotBasins(my_settings.fDirSave, dataMetr, fSeg, mesh)
      
    if (printTiming):
      tEnd = datetime.now()
      print "Time doSeg: ", tEnd-tStart
  dataMetr.close()
  
  #spatial metrics ------------------------
  dataMetr = netCDF4.Dataset(fMetr,'r')
  dataSeg = netCDF4.Dataset(fSeg,'r')
  if (my_settings.doMetrics):
    if (printTiming):
      tStart = datetime.now()
      
    basinMetrics.run_metrics(fMetrics, info, mesh, dataMetr, dataSeg, 0, nTimes-1)
    
    if (printTiming):
      tEnd = datetime.now()
      print "Time doMetrics: ", tEnd-tStart
  dataMetr.close()
  dataSeg.close()
  
  #basinMetrics.print_metrics(fMetrics)
  
  #time correspondence -----------------
  dataMetr = netCDF4.Dataset(fMetr,'r')
  dataSeg = netCDF4.Dataset(fSeg,'r')
  dataMetrics = netCDF4.Dataset(fMetrics, 'r')
  if (my_settings.doCorr):
    if (printTiming):
      tStart = datetime.now()
    
    correspond.run_correspond(fCorr, dataMetr, dataSeg, mesh, my_settings.deltaT, my_settings.trackMinMaxBoth, my_settings.areaOverlap, 0, nTimes-1, dataMetrics)
    if (False):
      correspond.plot_correspondences(my_settings.fDirSave, fCorr, nTimes-1, mesh)
      
    if (printTiming):
      tEnd = datetime.now()
      print "Time doCorr: ", tEnd-tStart
  dataMetrics.close()
  dataSeg.close()
  dataMetr.close()
  
  #time tracks -------------------------
  if (my_settings.doTracks):
    if (printTiming):
      tStart = datetime.now()
      
    #since appending to fTrack over time, wipe file before starting (if it exists)
    my_settings.silentremove(fTrack)
    tracks.run_tracks_timeInterval(fTrack, fCorr, 0, nTimes-1, timeStartGlobal, deltaTGlobal, fMetrics=fMetrics, trackOnlyMajor=True)
    if (False):
      tracks.plot_tracks_metrics(fTrack, my_settings.fDirSave+'test_tracks.png')
      #tracks.plot_tracks_cells(fTrack, mesh, my_settings.fDirSave)
      
    if (printTiming):
      tEnd = datetime.now()
      print "Time doTracks: ", tEnd-tStart
  #time metrics ----------------------
  
  commWorld.Free()

def demo_plotTracks():
  fTrack = my_settings.fTrack
  tracks.plot_tracks_metrics(fTrack, 'test_tracks.png')

def demo_algo_plots():
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
  #if already processed input data
  if (my_settings.inputType=='eraI'):
    mesh, cell0 = preProcess.demo_eraI(fMesh, [], fMetr, my_settings.rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData)
  elif (my_settings.inputType=='mpas'):
    mesh, cell0 = preProcess.demo_mpas(fMesh, [], fMetr, my_settings.rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData)
  elif (my_settings.inputType=='wrf_trop'):
      mesh, cell0 = preProcess.demo_wrf_trop(fMesh, [], fMetr, rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData, my_settings.fileMap)
  else:
    print "Unrecognized input type in my_settings: ",my_settings.inputType
  
  #segment --------------------------
  if (True):
    dataMetr = netCDF4.Dataset(fMetr,'r');
    segment.run_plotBasins(my_settings.fDirSave, dataMetr, fSeg, mesh)
    
    dataMetr.close()
  
  #time correspondence -----------------
  if (False):
    dataMetr = netCDF4.Dataset(fMetr,'r')
    dataSeg = netCDF4.Dataset(fSeg,'r')
    dataMetrics = netCDF4.Dataset(fMetrics, 'r')
    
    nTimes = len(dataSeg.dimensions['time'])
    correspond.plot_correspondences(my_settings.fDirSave, fCorr, nTimes-1, mesh)
    
    dataMetrics.close()
    dataSeg.close()
    dataMetr.close()
  
  #time tracks -------------------------
  if (True):
    tracks.plot_tracks_metrics(fTrack, my_settings.fDirSave+'test_tracks.png')
  
  #time metrics ----------------------

def debug_helper():
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
  #if already processed input data
  if (my_settings.inputType=='eraI'):
    mesh, cell0 = preProcess.demo_eraI(fMesh, [], fMetr, my_settings.rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData)
  elif (my_settings.inputType=='mpas'):
    mesh, cell0 = preProcess.demo_mpas(fMesh, [], fMetr, my_settings.rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData)
  elif (my_settings.inputType=='wrf_trop'):
      mesh, cell0 = preProcess.demo_wrf_trop(fMesh, [], fMetr, rEarth, dRegion, latThresh, my_settings.iTimeStart_fData, my_settings.iTimeEnd_fData, None)
  else:
    print "Unrecognized input type in my_settings: ",my_settings.inputType
    
  cells = np.arange(1000,1010)
  print mesh.get_latLon_inds(cells)
  print mesh.get_area_inds(cells)
  print mesh.isIndsInRegion(cells)

if __name__=='__main__':
  demo()
  #demo_testMods()
  #debug_helper()
  #demo_algo_plots()
  #tracks.plot_tracks_metrics(my_settings.fTrack, my_settings.fDirSave+'test_tracks.png')
  #tracks.demo_plotMetrics('/data02/cases/test_segment/testUnified/summer2006/tracks_debug.txt')
  #tracks.demo_plotLifetimes('/data02/cases/test_segment/testUnified/summer2006/tracks_debug.txt')
  #tracks.demo_compareMetrics('/data02/cases/test_segment/testUnified/200608/tracks_debug.txt')
