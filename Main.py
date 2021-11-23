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
    Root = '../proc/'
    #RootDirectoryWithSlcFolders = f'/esd/esd02/data/radar_data_vol1/Switzerland/Lauterbrunne/201207/20120730LB02_Final/slc/20120730/'
    RootDirectoryWithSlcFolders = f'/esd/esd02/data/radar_data_vol1/Switzerland/Lauterbrunnen/201207/20120730LB02_Final/slc/20120730/'
    IntProcFolder = f'{Root}/ints/'
    GammaConfigFileFullPath = 'GAMMA_config4Python.py'

    ## Load Processing Paramters from GammaConfigFile
    GetProcessingParametersFromFile(GammaConfigFileFullPath)

    # Get SLC Structure with lower antenna removed
    SlcStructure = GetSlcStructure(RootDirectoryWithSlcFolders)
    SlcStructure = RemoveLowerAntenna(SlcStructure)
    #
    #Get Master-rectified SLC for output tifs as background:
    GetMasterRecMli(SlcStructure["SlcFullPath"][10])
    #
    # ## Get and Int Structure:
    # ##-------------------------------------------------------
    IntStructure = SetupIntStructure(IntProcFolder)
    ## Check different intererograms with common temporal baseline
    TemporalBaselineInHours = 12.0
    IntStructure = tStructure(IntStructure, SlcStructure,TemporalBaselineInHours)
    
    for Index,x in enumerate(IntStructure["IntID"]):
        print(x)
        print(IntStructure["CenterTime"][Index])
    ### Calculate some random ints (this should be the ones to be stacked.)
    ### ------------------------------------------------------------------------
    # GetInterferogramm(IntStructure,5,1,0,0)
    # GetInterferogramm(IntStructure,6,1,0,0)
    # GetInterferogramm(IntStructure,7,1,0,0)
    # GetInterferogramm(IntStructure,8,1,0,0)
    ### Isolate the interferograms to be stacked from IntStructure
    #IntStructure4Stacking = GetSubsetFromIntStructure([5, 6, 7, 8],IntStructure)
    ### Stack all ints in IntStructure4Stacking. Outputfilename can also be automatically generated (e.g., with times..)
    #StackInts(IntStructure4Stacking,'../proc/stacks/StackedInt.int')

    ### Query a box in a Float32 (e.g. a float64 after cpx_to_real)
    #Array = ArrayFromFloat32(FileName,Width)
    ### Make sure that RowBox1<Robox2<Width and ColBox1<Colbox2<Height
    #BoxFromArray = Array[RowBox1:RowBox2,ColBox1:ColBox2]
    #f_handle = open('BoxOut.txt', 'w')
    #np.savetxt(f_handle, BoxFromArray, fmt="%.3f")



main()
