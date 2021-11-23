import os
import numpy as np
from datetime import date, datetime, timedelta
from GPRIBatchProcessFunctions import *
import pickle
import time
import pylab as plt

## In Python3 Shell: exec(open('Main.py').read())

def main():
    ## All Folderst with "/" at the end
    Root = '../proc2/'
    #RootDirectoryWithSlcFolders = f'/esd/esd02/data/radar_data_vol1/Switzerland/Lauterbrunne/201207/20120730LB02_Final/slc/20120730/'
    RootDirectoryWithSlcFolders = f'/esd/esd02/data/radar_data_vol1/Switzerland/Lauterbrunnen/201202_201207/'
    IntProcFolder = f'{Root}/ints/'
    GammaConfigFileFullPath = 'GAMMA_config4Python.py'

    ## Load Processing Paramters from GammaConfigFile
    GetProcessingParametersFromFile(GammaConfigFileFullPath)

    # Get SLC Structure with lower antenna removed
    SlcStructure = GetSlcStructure(RootDirectoryWithSlcFolders)
    SlcStructure = RemoveLowerAntenna(SlcStructure)

    ##Print contents of SLC structure to get a feel for it
    for i, x in enumerate(SlcStructure["SlcDate"]):
        print(f'This is the SLC date: {i} -- {x}')
    #
    #Get Master-rectified SLC for output tifs as background:
    GetMasterRecMli(SlcStructure["SlcFullPath"][1])
    #
    # ## Get and Int Structure:
    # ##-------------------------------------------------------
    IntStructure = SetupIntStructure(IntProcFolder)
    ## We don't care about temporal Baselin here. Choose your Index1 and Index2 wisely
    IntStructure = AddToIntStructure(IntStructure,SlcStructure,0,1)

    GetInterferogramm(IntStructure,0,1,0,0)





main()
