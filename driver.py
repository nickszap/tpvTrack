
import matplotlib.pyplot as plt

import my_settings
import preProcess

def demo():

  filesData = my_settings.filesData
  fMesh = filesData[0]  
  fMetr = my_settings.fDirSave+'fields.nc'
  
  rEarth = my_settings.rEarth
  dRegion = my_settings.dFilter
  latThresh = my_settings.latThresh
  
  #mesh = preProcess.demo_eraI(fMesh, filesData, fMetr, my_settings.rEarth, dRegion, latThresh)
  mesh = preProcess.demo_eraI(fMesh, [], fMetr, my_settings.rEarth, dRegion, latThresh)
  
  if (True):
    print "inRegion: ", mesh.inRegion;
    #plt.figure(); plt.plot(mesh.areaCell); plt.show()
    print "region: ", mesh.inDiskLat[0:6]


if __name__=='__main__':
  demo()
