from datetime import date,datetime,timedelta
import numpy as np
import subprocess
from colorama import Fore, Style
import os
import imp
import glob
import pickle
import sys
import time
import pylab  as plt
from PIL import Image

PathGamma = '/usr/share/modules/Modules/3.2.10/GAMMA_SOFTWARE-20190613/'
#PathGamma = '/usr/share/modules/Modules/3.2.10/GAMMA_SOFTWARE-20200728/'
sys.path.insert(0,PathGamma)

try:
    import py_gamma as pg
    print(f'{Fore.RED} Seeting path to: {PathGamma}. Make sure this is correct.')
except:
    print('{Fore.RED} Could Not Find GAMMA Python module. {Style.RESET_ALL}')

def GetSlcStructure(RootDirectoryWithSlcFolders):
    #Loops through all slc directories and collects all unique SLCs and their Attributes
    SlcStructure = {"SlcName" : [], "SlcID" : [], "SlcRxAntenna": [], "SlcParFullPath": [],"SlcFullPath": [],"SlcDate": []}
    Count = 0;CountFalse=0;
    for root, dirs, files in os.walk(RootDirectoryWithSlcFolders):
        for name in sorted(files):
            #Only take the SLC Files and nothing else
            if (name[-4::]==".slc"):
                if(FileDoesExist(f'{root}/{name[:-4]}.slc.par')):
                    #Remove double entries
                    if (name not in SlcStructure["SlcName"]):
                        SlcStructure["SlcName"].append(name)
                        SlcStructure["SlcID"].append(name[:-5])
                        SlcStructure["SlcRxAntenna"].append(name[-5:-4])
                        SlcStructure["SlcFullPath"].append(root+'/'+name)
                        SlcStructure["SlcParFullPath"].append(root+'/'+name+'.par')
                        SlcStructure["SlcDate"].append(datetime.strptime(name.split('.')[0][:-1], '%Y%m%d_%H%M%S'))
                else:
                    CountFalse = CountFalse + 1


    NumberOfEntries = len(SlcStructure["SlcID"])
    print(f"Added {NumberOfEntries} unique SLCs to the list")
    print(f"Skipped {CountFalse} entries because of missing slc.par files.")
    return SlcStructure

def GetMasterRecMli(SlcFullPath,Scale=1.0):
    #Get rectified version of SLC. The output for Master is defined in GammaConfigFile.
    MasterRecFullPath=params.MasterRecMliFullPath
    if os.path.isfile(MasterRecFullPath):
        print(f'{Fore.GREEN} We already have a Master MLI. Do nothing.{Style.RESET_ALL}')
    else:
        MasterRecSlcFolder, SlcName = os.path.split(MasterRecFullPath)
        SlcID = os.path.splitext(SlcName)
        try:
            os.stat(MasterRecSlcFolder)
            print(f'{Fore.GREEN} Master rectified SLC: {MasterRecSlcFolder} exists. {Style.RESET_ALL}')
        except:
            os.makedirs(MasterRecSlcFolder)
            print(f'{Fore.GREEN} Created {MasterRecSlcFolder} {Style.RESET_ALL}')

        cmd_Slc2Mli = f'multi_look {SlcFullPath} {SlcFullPath}.par {MasterRecSlcFolder}/Master.mli {MasterRecSlcFolder}/Master.mli.par 1 1 >> GammaProcLogFiles.log'
        subprocess.check_call(cmd_Slc2Mli, shell=True)
        cmd_MasterRec = f'pol2rec  {MasterRecSlcFolder}/Master.mli {MasterRecSlcFolder}/Master.mli.par {MasterRecFullPath} {MasterRecFullPath}.par {params.pix_size_p2r} 0 0 - - - - - >> GammaProcLogFiles.log'
        subprocess.check_call(cmd_MasterRec, shell=True)
        MliWidth = GetKeywordFromTxtFile(f'{MasterRecSlcFolder}/Master.mli.par','range_samples')
        if (Scale != "1.0"):
            cmd_Scale = f'resamp_image {MasterRecFullPath} {MliWidth} {Scale} {Scale} {MasterRecSlcFolder}/Master.rec.scaled - - >> GammaProcLogFiles.log'
            subprocess.check_call(cmd_Scale, shell=True)
        MliWidth = GetKeywordFromTxtFile(f'{MasterRecFullPath}.par','range_samples')
        cmd_MasterRecTif = f'raspwr {MasterRecFullPath} {MliWidth} - - - - - - -  {MasterRecFullPath}.tif >> GammaProcLogFiles.log'
        subprocess.check_call(cmd_MasterRecTif, shell=True)
        print(f"{Fore.GREEN} Master rectified SLC created {MasterRecFullPath} {Style.RESET_ALL}")

def GetSubsetFromIntStructure(SubsetIndices,IntStructure):
    NumberOfEntriesInitial = len(IntStructure["IntName"])
    for key in IntStructure:
        if (not(key=="IntFolderTifs" or key == "IntFolder")):
            #Removing Elements from List starting with the largest index.
            #for i in sorted(SubsetIndices, reverse=True):
            for i in sorted(np.arange(0,NumberOfEntriesInitial), reverse=True):
                if i not in SubsetIndices:
                    del IntStructure[key][i]
    NumberOfEntries = len(IntStructure["IntName"])
    print(f"Got subset of {NumberOfEntries} ints from {NumberOfEntriesInitial} total entries.")
    return IntStructure

def RemoveFromSlcStructure(FilterIndices,SlcStructure):
    #Remove Elements from all Lists in SlcStructure based on Indices
    #All Lists must be of same size.
    NumberOfEntries = len(SlcStructure["SlcName"])
    print(f"Old SlcStructure had {NumberOfEntries} entries.")
    for key in SlcStructure:
        #Removing Elements from List starting with the largest index.
        for i in sorted(FilterIndices, reverse=True):
            #print(i)
            del SlcStructure[key][i]
    NumberOfEntries = len(SlcStructure[key])
    print(f"New SlcStructure has {NumberOfEntries} entries.")
    return SlcStructure

def RemoveFromIntStructure(FilterIndices,IntStructure):
    #Remove Elements from all Lists in IntStructure based on Indices
    #All Lists must be of same size.
    NumberOfEntries = len(IntStructure["IntName"])
    print(f"Old IntStructure had {NumberOfEntries} entries.")
    for key in IntStructure:
        if (not(key=="IntFolderTifs" or key == "IntFolder")):
            #Removing Elements from List starting with the largest index.
            for i in sorted(FilterIndices, reverse=True):
                del IntStructure[key][i]
    NumberOfEntries = len(IntStructure["IntName"])
    print(f"New IntStructure has {NumberOfEntries} entries.")
    return IntStructure

def RemoveFromStackStructure(FilterIndices,StackStructure):
    #Remove Elements from all Lists in SlcStructure based on Indices
    #All Lists must be of same size.
    NumberOfEntries = len(StackStructure["StackID"])
    print(f"Old IntStructure had {NumberOfEntries} entries.")
    for k in FilterIndices:
        print(f'Removing {StackStructure["StackID"][k]}')
    for key in StackStructure:
        if (not(key=="IntFolderTifs" or key == "IntFolder")):
            #Removing Elements from List starting with the largest index.
            for i in sorted(FilterIndices, reverse=True):
                del StackStructure[key][i]
    NumberOfEntries = len(StackStructure["StackID"])
    print(f"New IntStructure has {NumberOfEntries} entries.")
    return StackStructure

def RemoveLowerAntenna(SlcStructure):
    #Remove All Slcs from from the lower Antenna
    FilterIndices = [int(i) for i, x in enumerate(SlcStructure["SlcRxAntenna"]) if x == "l"]
    SlcStructure=RemoveFromSlcStructure(FilterIndices,SlcStructure)
    print(f'{Fore.GREEN} RemoveLowerAntenna: Removed {len(FilterIndices)} SLCs from SlcStructure. {Style.RESET_ALL}')
    return SlcStructure

def RemoveSlcsInTimeInterval(SlcStructure,TimeStart,TimeEnd):
    #Remove all SLCs in certain time interval (rotated antenna)
    FilterIndices = [int(i) for i, x in enumerate(SlcStructure["SlcDate"]) if (((TimeStart-x).total_seconds()<=0) and (TimeEnd-x).total_seconds()>=0)]
    SlcStructure=RemoveFromSlcStructure(FilterIndices,SlcStructure)
    print(f'{Fore.GREEN} RemoveSlcsInTimeInterval: Removed {len(FilterIndices)} SLCs from SlcStructure in time interval: {TimeStart} and {TimeEnd} {Style.RESET_ALL}')
    return SlcStructure

def FindClosestTimeInSlcStructure(TargetTime,SlcStructure):
    #Finds closest SLC to TargetTime in SlcStructure
    TmpDiffSeconds = 1.0e12;
    for i, x in enumerate(SlcStructure["SlcDate"]):
        if (np.abs((x - TargetTime).total_seconds()) < TmpDiffSeconds):
            TmpDiffSeconds = np.abs((x - TargetTime).total_seconds())
            FitIndex = i; FitDatetime = x; FitDifferenceSeconds = TmpDiffSeconds
    return FitIndex,FitDatetime,FitDifferenceSeconds

def GetMultiTemporalIntStructure(SlcStructure,IntStructure,CenterTime,TimeIntervalInHours):
    ### Get multitemporal Interferogramstructure
    from itertools import combinations
    T1 = CenterTime - timedelta(hours=TimeIntervalInHours/2.0)
    T2 = CenterTime + timedelta(hours=TimeIntervalInHours/2.0)
    IndexRangeSlc=[]
    for i, x in enumerate(SlcStructure["SlcDate"]):
            if (T1 < x < T2): IndexRangeSlc.append(i)
    print(f'{Fore.GREEN} GetMultiTemporalIntStructure: Found {len(IndexRangeSlc)} SLCs in requested time interval. {Style.RESET_ALL}')
    IndexCombinations = list(combinations(IndexRangeSlc,2))
    print(f'{Fore.GREEN} GetMultiTemporalIntStructure: This results in {len(IndexCombinations)} possible multi-temporal interferograms. {Style.RESET_ALL}')
    [AddToIntStructure(IntStructure,SlcStructure,Inds[0],Inds[1]) for Inds in IndexCombinations]
    return IntStructure

def GetStackFileId(MultiTemporalIntStructure):
    ## Find Start and End Dates of Stack for systematic naming.
    for k,x in enumerate(MultiTemporalIntStructure["IntFullPath"]):
        if k==0:
            T1=MultiTemporalIntStructure["StartTime"][k]
            T2=MultiTemporalIntStructure["EndTime"][k]
        else:
            if MultiTemporalIntStructure["StartTime"][k] < T1:
                T1 = MultiTemporalIntStructure["StartTime"][k]
            if MultiTemporalIntStructure["EndTime"][k] > T2:
                T2 = MultiTemporalIntStructure["EndTime"][k]
        StackFileId = f'Stack_{T1.year}{T1.month}{T1.day}_{T1.hour}{T1.minute}{T1.second}_{T2.year}{T2.month}{T2.day}_{T2.hour}{T2.minute}{T2.second}'
    return StackFileId

def SetupStackStructure():
    StackStructure = {"StackFolder" : [],  \
                     "StackID" : [],      \
                     "StackFullPathInt": [], \
                     "StackOffFileFullPath": [], \
                     "StackSlcFileFullPath": [], \
                     "StackFullPathIntW": [], \
                     "StackFullPathUnw": [], \
                     "StackFullPathStacked": [], \
                     "StackFullPathUnwW": [], \
                     "StackFullPathDispmap": [], \
                     "CenterTime": [],    \
                     "StartTime": [],     \
                     "NumberOfStacks": [],\
                     "EndTime": []}
    return StackStructure

def AddToStackStructure(StackStructure,MultiTemporalIntStructure,StackFolder):
    for k,x in enumerate(MultiTemporalIntStructure["IntFullPath"]):
        if k==0:
            T1=MultiTemporalIntStructure["StartTime"][k]
            T2=MultiTemporalIntStructure["EndTime"][k]
        else:
            if MultiTemporalIntStructure["StartTime"][k] < T1:
                T1 = MultiTemporalIntStructure["StartTime"][k]
            if MultiTemporalIntStructure["EndTime"][k] > T2:
                T2 = MultiTemporalIntStructure["EndTime"][k]
    T1day = "{:0>2}".format(T1.day);T2day = "{:0>2}".format(T2.day)
    T1hour = "{:0>2}".format(T1.hour);T2hour = "{:0>2}".format(T2.hour)
    T1min = "{:0>2}".format(T1.minute);T2min = "{:0>2}".format(T2.minute)
    T1sec = "{:0>2}".format(T1.second);T2sec = "{:0>2}".format(T2.second)

    StackID = f'Stack_{T1.year}{T1.month}{T1day}_{T1hour}{T1min}{T1sec}_{T2.year}{T2.month}{T2day}_{T2hour}{T2min}{T2sec}'
    StackStructure["StackFolder"].append(StackFolder)
    StackStructure["StackID"].append(StackID)
    StackStructure["StartTime"].append(T1)
    StackStructure["EndTime"].append(T2)
    StackStructure["StackFullPathInt"].append(f'{StackFolder}{StackID}.int')
    StackStructure["StackFullPathIntW"].append(f'{StackFolder}{StackID}W.int')
    StackStructure["StackFullPathUnw"].append(f'{StackFolder}{StackID}.unw')
    StackStructure["StackFullPathUnwW"].append(f'{StackFolder}{StackID}W.unw')
    StackStructure["StackFullPathStacked"].append(f'{StackFolder}{StackID}.stacked.unw')
    StackStructure["StackFullPathDispmap"].append(f'{StackFolder}{StackID}.dispmap')
    StackStructure["NumberOfStacks"].append(f'{k}')
    StackStructure["StackOffFileFullPath"].append(MultiTemporalIntStructure["IntOffFullPath"][0])
    StackStructure["StackSlcFileFullPath"].append(MultiTemporalIntStructure["SlcPar1FullPath"][0])
    StackStructure["CenterTime"].append(T1 + (T2-T1)/2)
    return StackStructure

def AddToIntStructure(IntStructure,SlcStructure,Index1,Index2):
    IntID = SlcStructure["SlcID"][Index1]+'_'+SlcStructure["SlcID"][Index2]
    IntName = IntID+'.int'
    UnwName = IntID+'.unw'
    IntProcFolder = IntStructure["IntFolder"][0]
    IntProcFolderTifs = IntStructure["IntFolderTifs"][0]
    IntStructure["IntID"].append(IntID)
    IntStructure["IntName"].append(IntName)
    IntStructure["IntOffFullPath"].append(IntProcFolder+IntID+'.off')
    IntStructure["IntAdfFullPath"].append(IntProcFolder+IntID+'.sm.int')
    IntStructure["IntFullPath"].append(IntProcFolder+IntName)
    IntStructure["UnwFullPath"].append(IntProcFolder+UnwName)
    IntStructure["Slc1FullPath"].append(SlcStructure["SlcFullPath"][Index1])
    IntStructure["Slc2FullPath"].append(SlcStructure["SlcFullPath"][Index2])
    IntStructure["SlcPar1FullPath"].append(SlcStructure["SlcParFullPath"][Index1])
    IntStructure["SlcPar2FullPath"].append(SlcStructure["SlcParFullPath"][Index2])
    TBaselineInHours = (SlcStructure["SlcDate"][Index2]-SlcStructure["SlcDate"][Index1]).total_seconds()/3600.0
    IntStructure["TemporalBaselineInHours"].append(TBaselineInHours)
    IntStructure["CenterTime"].append(SlcStructure["SlcDate"][Index1]+timedelta(hours=TBaselineInHours/2.0))
    IntStructure["StartTime"].append(SlcStructure["SlcDate"][Index1])
    IntStructure["EndTime"].append(SlcStructure["SlcDate"][Index2])
    IntStructure["IntRecFullPath"].append(IntProcFolder+IntName+'.rec')
    IntStructure["IntRecParFullPath"].append(IntProcFolder+IntName+'.rec.par')
    return IntStructure

def GetDecimalDayOfYear(DTime):
    DayOfYear = DTime.timetuple().tm_yday
    hour = DTime.hour;minute = DTime.minute; secs = DTime.second
    DecimalDay = (hour+minute/60.0 + secs/3600.0)/24.0
    DecimalDayOfYear=DayOfYear+DecimalDay
    return DecimalDayOfYear

def DeleteIntsInStructure(IntStructure):
    for i,x in enumerate(["IntFullPath"]):
        try:
            print(f'Deleting Interferogram: {x}')
            os.remove(f'{x}')
            print(f'Deleting Interferogram: {IntStructure["IntOffFullPath"][i]}')
            os.remove(f'{IntStructure["IntOffFullPath"][i]}')
        except:
            print(f'Could not delete: {x}')

def SetupIntStructure(ProcFolder):
    IntStructure = {"IntName" : [],        \
                    "IntFolder" : [],        \
                    "IntFolderTifs" : [],        \
                    "IntID" : [],          \
                    "IntFullPath": [],     \
                    "UnwFullPath": [],     \
                    "IntOffFullPath" : [], \
                    "IntAdfFullPath" : [], \
                    "IntRecFullPath": [],     \
                    "IntRecParFullPath": [],     \
                    "Slc1FullPath" : [],   \
                    "SlcPar1FullPath": [], \
                    "Slc2FullPath" : [],   \
                    "SlcPar2FullPath": [],  \
                    "CenterTime": [],  \
                    "StartTime": [],  \
                    "EndTime": [],  \
                    "TemporalBaselineInHours": []}
    IntStructure["IntFolder"].append(ProcFolder)
    IntStructure["IntFolderTifs"].append(f'{ProcFolder[0:-1]}Tifs')
    MakeFolder(ProcFolder)
    MakeFolder(f'{ProcFolder[:-1]}Tifs')
    return IntStructure

def FindClosestDoyStartTime(IntStructure,TargetDoy):
    MinDiff = 1.0e30
    MinIndex = 0
    for i,x in enumerate(IntStructure["StartTime"]):
        Doy = GetDecimalDayOfYear(IntStructure["StartTime"][i])
        if np.abs(Doy-TargetDoy) < MinDiff:
            MinDiff = np.abs(Doy-TargetDoy)
            MinIndex = i
    print(f'Closest Doy is {MinDiff} away from TargetDoy.')
    return MinDiff,MinIndex

def GetInterferogramm(IntStructure,Index,WithPol2Rec=1,WithMatching=0,WithFileCheck=1):
    #Generate an Interferogram from an IntStructure with or without matching.
    #Generate Quicklookf if WithPol2Rec=1
    print(f'')
    if (WithFileCheck == 1 and os.path.isfile(IntStructure["IntFullPath"][Index])):
        print(f'{Fore.RED} File {IntStructure["IntFullPath"][Index]} exists. Do nothing. {Style.RESET_ALL}')
    else:
        cmd_CreateOffset = f'printf \'\\n\\n\\n\\n\\n\\n\\n\' |create_offset {IntStructure["SlcPar1FullPath"][Index]} {IntStructure["SlcPar2FullPath"][Index]} {IntStructure["IntOffFullPath"][Index]} {params.algorithm} >> GammaProcLogFiles.log 2>&1'
        subprocess.check_call(cmd_CreateOffset, shell=True)
        Slc1= IntStructure["Slc1FullPath"][Index]
        Slc1Par= IntStructure["SlcPar1FullPath"][Index]
        Slc2= IntStructure["Slc2FullPath"][Index]
        Slc2Par= IntStructure["SlcPar2FullPath"][Index]
        OffPar=IntStructure["IntOffFullPath"][Index]
        Int=IntStructure["IntFullPath"][Index]
        if (WithMatching == 1):
            print(f'{Fore.GREEN} Doing Int {Index} with matching: {Int} {Style.RESET_ALL}')
            pg.init_offset(Slc1,Slc2,Slc1Par,Slc2Par,OffPar,1,1,'-','-','-','-','-','-','-','-',stdout_flag=False)
            pg.offset_pwr(Slc1,Slc2,Slc1Par,Slc2Par,OffPar,'Offs.tmp','Ccp.tmp',stdout_flag=False)
            pg.offset_fit('Offs.tmp','Ccp.tmp',OffPar,'Coffs.tmp','Coffsets.tmp',stdout_flag=False)
            pg.SLC_interp(Slc2,Slc1Par,Slc2Par,OffPar,'RSlc.slc','RSlc.slc.par',stdout_flag=False)
            pg.SLC_intf(Slc1,'RSlc.slc',Slc1Par,'RSlc.slc.par',OffPar,Int,params.rlks_intf,params.azlks_intf,params.loff_intf,'-',params.sps_flg,params.azf_flg,params.rp1_flg,params.rp2_flg,stdout_flag=False)
            subprocess.check_call('rm Offs.tmp Ccp.tmp Coffs.tmp Coffsets.tmp RSlc.slc RSlc.slc.par',shell=True)
        else:
            print(f'{Fore.GREEN} Doing Int {Index} without matching: {Int} {Style.RESET_ALL}')
            pg.SLC_intf(Slc1,Slc2,Slc1Par,Slc2Par,OffPar,Int,params.rlks_intf,params.azlks_intf,params.loff_intf,'-',params.sps_flg,params.azf_flg,params.rp1_flg,params.rp2_flg,stdout_flag=False)
        if WithPol2Rec == 1:
            #print(f'{Fore.GREEN} Getting a rectangular version.. {Style.RESET_ALL}')
            cmd_Pol2Rec = f'pol2rec {IntStructure["IntFullPath"][Index]} {IntStructure["SlcPar1FullPath"][Index]} {IntStructure["IntRecFullPath"][Index]} {IntStructure["IntRecParFullPath"][Index]} {params.pix_size_p2r} 1 0 - - - - - >>GammaProcLogFiles.log'
            subprocess.check_call(cmd_Pol2Rec, shell=True)
            RecWidth = GetKeywordFromTxtFile(IntStructure["IntRecParFullPath"][Index],'range_samples')
            Tmp, OutTifName = os.path.split(IntStructure["IntRecFullPath"][Index])
            OutTifNameFullPath = f'{IntStructure["IntFolderTifs"][0]}/{OutTifName}.tif'
            if os.path.isfile(params.MasterRecMliFullPath):
                cmd_RasMphPwrInt = f'rasmph_pwr {IntStructure["IntRecFullPath"][Index]} {params.MasterRecMliFullPath}  {RecWidth} - - - 1 1 - - - {OutTifNameFullPath} >> GammaProcLogFiles.log '
                subprocess.check_call(cmd_RasMphPwrInt, shell=True)
            else:
                cmd_RasMphInt = f'rasmph {IntStructure["IntRecFullPath"][Index]} {RecWidth} - -  1 1 - - - {OutTifNameFullPath} >> GammaProcLogFiles.log'
                subprocess.check_call(cmd_RasMphInt, shell=True)
            print(f'{Fore.GREEN} Done with Int {Index} Pol2Rec {IntStructure["IntName"][Index]} {Style.RESET_ALL}')
            print(f'{Fore.GREEN} Check with:  {OutTifNameFullPath} {Style.RESET_ALL}')
        else:
            print(f'{Fore.GREEN} Done with Int {Index} {IntStructure["IntName"][Index]} and not Pol2Rec. {Style.RESET_ALL}')

def UnwrapInterf(IntIn,UnwOut,Width=22434,CC_Threshold=0.7):
    if os.path.isfile(UnwOut):
        print(f'{Fore.RED} Unw file: {UnwOut} exists. Do nothing. {Style.RESET_ALL}')
    else:
        cmd_adf = f'adf {IntIn} {IntIn[:-4]}.sm.int {IntIn[:-4]}.cc {Width} {params.alpha} {params.nfft} {params.cc_win} - - - {params.wfrac} >> GammaProcLogFiles.log'
        try:
            cmd_unw = f'mcf {IntIn[:-4]}.sm.int {IntIn[:-4]}.cc {params.UnwMask} {UnwOut} {Width} {params.tri_mode_mcf} - - - - {params.npat_r_mcf} {params.npat_az_mcf} -  {params.xCoordPhUnw} {params.yCoordPhUnw} 1 >> GammaProcLogFiles.log '
            print(f'{Fore.GREEN} Unwrapping with Mask: {params.UnwMask} {Style.RESET_ALL}')
            subprocess.check_call(cmd_adf,shell=True)
            subprocess.check_call(cmd_unw,shell=True)
        except:
            cmd_unw = f'mcf {IntIn[:-4]}.sm.int {IntIn[:-4]}.cc lmask.ras {UnwOut} {Width} {params.tri_mode_mcf} - - - - {params.npat_r_mcf} {params.npat_az_mcf} - {params.xCoordPhUnw} {params.yCoordPhUnw} 1 >> GammaProcLogFiles.log'
            print(f'{Fore.GREEN} Adf fringe filtering..  {Style.RESET_ALL}')
            subprocess.check_call(cmd_adf,shell=True)
            print(f'{Fore.GREEN} Making a mask.. {Style.RESET_ALL}')
            pg.rascc_mask(f'{IntIn[:-4]}.cc','-',f'{Width}','-','-','-','-','-',f'{CC_Threshold}','-','-','-','-','-','-','lmask.ras')
            print(f'{Fore.GREEN} Unwrapping: {IntIn}....{Style.RESET_ALL}')
            subprocess.check_call(cmd_unw,shell=True)
        print(f'Done.')

def Unw2Hgt(UnwIn,HgtOut,RefHeight=0.0):
    print(f'{Fore.RED} Converting UNW to HGT with python2 function. Using Parameters from Master-Rec File. {Style.RESET_ALL}')
    print(f'Using: AuxFunctions/PythonFunctions/GPRI_hgtmap/gpri2_hgt.py')
    if os.path.isfile(params.MasterRecMliFullPath):
        MasterRecFolder, MasterRec = os.path.split(params.MasterRecMliFullPath)
        cmd_Unw2Hgt = f'AuxFunctions/PythonFunctions/GPRI_hgtmap/gpri2_hgt2.py {UnwIn} {MasterRecFolder}/Master.mli.par {params.xCoordPhUnw} {params.yCoordPhUnw} {RefHeight} {HgtOut} -m 2 >> GammaProcLogFiles.log'
        subprocess.check_call(cmd_Unw2Hgt,shell=True)
    else:
        print(f'{Fore.RED} Could not find Master-MLI File. Get it with GetMasterRecMli(). {Style.RESET_ALL}')

def GetDifferentialInterfStructure(IntStructure, SlcStructure, TempBaselineHours, MaxTemporalMisfitSecs=300, PickleFile=f'PickleFiles/DifferntialIntStructure.pickle'):
    #Get Interferogramm Pairs with Common Master and Temporal baselin (backward/forward)
    if FileDoesExist(PickleFile):
        print(f'{Fore.RED} Loading Pickle of Differential Int Structure from previous run. {Style.RESET_ALL}')
        IntStructure = LoadPickle(f'{PickleFile}')
    else:
        for i, x in enumerate(SlcStructure["SlcDate"]):
            if (i%1==0):
                TargetTime1 = x + timedelta(hours=TempBaselineHours)
                TargetTime2 = x + timedelta(hours=2*TempBaselineHours)
                [Index1, FitDateTime,FitDifferenceSeconds] = FindClosestTimeInSlcStructure(TargetTime1,SlcStructure)
                [Index2, FitDateTime2,FitDifferenceSeconds2] = FindClosestTimeInSlcStructure(TargetTime2,SlcStructure)
                if (FitDifferenceSeconds < MaxTemporalMisfitSecs and FitDifferenceSeconds2 < MaxTemporalMisfitSecs):
                    IntStructure = AddToIntStructure(IntStructure,SlcStructure,i,Index1)
                    IntStructure = AddToIntStructure(IntStructure,SlcStructure,Index1,Index2)
        SavePickle(IntStructure,PickleFile)
    return IntStructure

def GetDifferentialInterfStructures(IntStructure1, IntStructure2, SlcStructure, TempBaselineHours, MaxTemporalMisfitSecs=300, PickleFile=f'PickleFiles/DifferntialIntStructure.pickle'):
    #Get Interferogramm Pairs with Common Master and Temporal baselin (backward/forward) and store them in two IntStructures
    if FileDoesExist(PickleFile):
        print(f'{Fore.RED} Loading Pickle of Differential Int Structure from previous run. {Style.RESET_ALL}')
        IntStructure1, IntStructure2 = LoadPickle(f'{PickleFile}')
    else:
        for i, x in enumerate(SlcStructure["SlcDate"]):
            if (i%1==0):
                TargetTime1 = x + timedelta(hours=TempBaselineHours)
                TargetTime2 = x + timedelta(hours=2*TempBaselineHours)
                [Index1, FitDateTime,FitDifferenceSeconds] = FindClosestTimeInSlcStructure(TargetTime1,SlcStructure)
                [Index2, FitDateTime2,FitDifferenceSeconds2] = FindClosestTimeInSlcStructure(TargetTime2,SlcStructure)
                if (FitDifferenceSeconds < MaxTemporalMisfitSecs and FitDifferenceSeconds2 < MaxTemporalMisfitSecs):
                    IntStructure1 = AddToIntStructure(IntStructure,SlcStructure,i,Index1)
                    IntStructure2 = AddToIntStructure(IntStructure,SlcStructure,Index1,Index2)
        SavePickle([IntStructure1, IntStructure2],PickleFile)
    return IntStructure

def tStructure(IntStructure,SlcStructure,TempBaselineHours, MaxTemporalMisfitSecs=300):
    #Get Interferogramm Pairs with fixed temporal baseline
    for i, x in enumerate(SlcStructure["SlcDate"]):
        if (i%1==0):
            TargetTime1 = x + timedelta(hours=TempBaselineHours)
            [Index1, FitDateTime,FitDifferenceSeconds] = FindClosestTimeInSlcStructure(TargetTime1,SlcStructure)
            if (FitDifferenceSeconds < MaxTemporalMisfitSecs):
                IntStructure = AddToIntStructure(IntStructure,SlcStructure,i,Index1)
    return IntStructure

def GetTopographicIntStructure(IntStructure,SlcStructure,TargetTime,TimeIntervalInHours):
    TargetTime1 = TargetTime - timedelta(hours=TimeIntervalInHours)/2
    TargetTime2 = TargetTime + timedelta(hours=TimeIntervalInHours)/2
    print(TargetTime1)
    print(TargetTime2)
    for i, x in enumerate(SlcStructure["SlcDate"]):
        if (TargetTime1 <= x <= TargetTime2):
            upper=f'{SlcStructure["SlcName"][i][:-5]}u.slc'
            lower=f'{SlcStructure["SlcName"][i][:-5]}l.slc'
            if ((upper in SlcStructure["SlcName"]) and (lower in SlcStructure["SlcName"])):
                Index1 = SlcStructure["SlcName"].index(upper)
                Index2 = SlcStructure["SlcName"].index(lower)
                IntID=f'{upper[:-5]}_{lower[:-5]}'
                #print(f'RD {IntID}')
                if (IntID not in IntStructure["IntID"]):
                    IntStructure = AddToIntStructure(IntStructure,SlcStructure,Index1,Index2)
    print(f'{Fore.GREEN} Added {len(IntStructure["IntName"])} topographic interferograms to the list. {Style.RESET_ALL}')

    return IntStructure

def CalculateDifferentialInterfs(IntStructure,TripleDiffFolder):
    ## Create Output Folder
    MakeFolder(TripleDiffFolder)
    ## Calculate all Interferograms in Instructure
    for i,x in enumerate(IntStructure["IntName"]):
        if os.path.isfile(IntStructure["IntFullPath"][i]):
            print('Int-File already exists. Do nothing.')
        else:
            GetInterferogramm(IntStructure,i,1,1)
    ## Difference successive interferograms
    for i,x in enumerate(IntStructure["IntName"]):
        if (i%2 ==  0):
            LDiffName = f'{TripleDiffFolder}Diff_{IntStructure["IntName"][i][0:-4]}_{IntStructure["IntName"][i+1]}'
            if os.path.isfile(LDiffName):
                print('Diff file already exists. Do nothing.')
            else:
                DifferenceInterferogramsPhase(IntStructure["IntFullPath"][i],IntStructure["IntFullPath"][i+1],LDiffName,1,IntStructure["SlcPar1FullPath"][0])

def CalculateDifferentialInterf(IntStructure,TripleDiffFolder,Index):
    ## Create Output Folder
    MakeFolder(TripleDiffFolder)
    ## Calculate all Interferograms in Instructure
    if os.path.isfile(IntStructure["IntFullPath"][Index]):
        print('Int-File already exists. Do nothing.')
    else:
        GetInterferogramm(IntStructure,Index,0,1)
    if os.path.isfile(IntStructure["IntFullPath"][Index+1]):
        print('Int-File already exists. Do nothing.')
    else:
        GetInterferogramm(IntStructure,Index+1,0,1)
    ## Difference successive interferograms
    LDiffName = f'{TripleDiffFolder}Diff_{IntStructure["IntName"][Index][0:-4]}_{IntStructure["IntName"][Index+1]}'
    LDiffNameGc = f'{LDiffName[:-4]}.gc.int'
    if os.path.isfile(LDiffName):
        print('Diff file already exists. Do nothing.')
    else:
        DifferenceInterferogramsPhase(IntStructure["IntFullPath"][Index],IntStructure["IntFullPath"][Index+1],LDiffName,1,IntStructure["SlcPar1FullPath"][0])
        pg.geocode_back(f'{LDiffName}','22434','AuxFiles/Geocode/Lut',f'{LDiffNameGc}','4024','4746')
        pg.data2geotiff('AuxFiles/Geocode/Seg.dem_par',f'{LDiffNameGc}',2,f'{LDiffNameGc[:-4]}.tif')
        os.system(f'gdal_translate -tr 10 10 {LDiffNameGc[:-4]}.tif {LDiffNameGc[:-4]}R.tif')
        os.system(f'gdalwarp -srcnodata "0" -dstnodata \"NaN\"  -s_srs EPSG:3031 -t_srs EPSG:4326 {LDiffNameGc[:-4]}R.tif {LDiffNameGc[:-4]}.geo.tif')

def ArrayFromFloat32(FileName,Width):
    arr = np.fromfile(FileName,'>f4')
    print(f'Length of Array: {len(arr)}')
    print(f'Width of Array: {Width}')
    print(f'Height of Array: {len(arr)/Width}')
    arr = np.reshape(arr,(np.int(len(arr)/Width),np.int(Width)))
    return arr

def DifferenceInterferogramsPhase(Int1, Int2, DiffInt, WithPol2Rec=0,IntSlc=None,ByteSwap=0):
    print(f' ')
    print(f'{Fore.GREEN} Differencing Phase of {Int1} and {Int2} {Style.RESET_ALL}')
    cmd_tmp1 = f'cpx_to_real {Int1} Int1Real.int 22434 4 >> GammaProcLogFiles.log'
    cmd_tmp2 = f'cpx_to_real {Int2} Int2Real.int 22434 4 >> GammaProcLogFiles.log'
    subprocess.check_call(cmd_tmp1, shell=True)
    subprocess.check_call(cmd_tmp2, shell=True)

    arr1 = np.fromfile('Int1Real.int',dtype='>f4')
    arr2 = np.fromfile('Int2Real.int',dtype='>f4')
    arr3 = np.mod(arr1-arr2,2*np.pi).astype('>f4')
    arr3.tofile(DiffInt,  format='>f4')
    #subprocess.check_call('rm -f Int1Real.int Int2Real.int',shell=True)
    if WithPol2Rec == 1:
        #print(f'{Fore.GREEN} Getting a rectangular version.. {Style.RESET_ALL}')
        cmd_Pol2Rec = f'pol2rec {DiffInt} {IntSlc} {DiffInt[0:-4]}.rec  {DiffInt[0:-4]}.rec.par  {params.pix_size_p2r} 0 0 - - - - - >>Main_Int.log'
        subprocess.check_call(cmd_Pol2Rec, shell=True)
        RecWidth = GetKeywordFromTxtFile(f'{DiffInt[0:-4]}.rec.par','range_samples')
        if os.path.isfile(params.MasterRecMliFullPath):
            cmd_RasRmgInt = f'rasrmg {DiffInt[0:-4]}.rec {params.MasterRecMliFullPath} {RecWidth} - - -  1 1 1  - - - - {DiffInt[0:-4]}.rec.tif - - - >> GammaProcLogFiles.log'
        else:
            cmd_RasRmgInt = f'rasrmg {DiffInt[0:-4]}.rec - {RecWidth} - - -  1 1 - - - - - {DiffInt[0:-4]}.rec.tif - - - >> GammaProcLogFiles.log'
        subprocess.check_call(cmd_RasRmgInt, shell=True)

    print(f"Done writing {DiffInt}")

def DifferenceDems(PathToDem1,PathToDem2,PathToOutFile):
    arr1 = np.fromfile(PathToDem1,dtype='>f4')
    arr2 = np.fromfile(PathToDem2,dtype='>f4')
    arr3 = (arr1-arr2).astype('>f4')
    arr3.tofile(PathToOutFile,  format='>f4')
    MakePol2RecTifFromHgt(PathToOutFile,'../proc_data/Diff.tif',50.1)

def GetDem(TopoIntStructure,TargetTime,TimeIntervalInHours,OutputFolder):
    DemIdFullPath=f"{OutputFolder}/Dem_{TargetTime.strftime('%Y%m%d_%H%M%S')}_TimeInterval_{TimeIntervalInHours}h"
    if (FileDoesNotExist(f'{DemIdFullPath}.hgt')):
        for i,x in enumerate(TopoIntStructure["IntName"]):
            terferogramm(TopoIntStructure,i)
        print(f'{Fore.GREEN} Start stacking... {Style.RESET_ALL}')
        StackInts(TopoIntStructure, f'{DemIdFullPath}.int',f'{DemIdFullPath}.rec.tif')
        print(f'{Fore.GREEN} Unwrapping... {Style.RESET_ALL}')
        UnwrapInterf(f'{DemIdFullPath}.int',f'{DemIdFullPath}.unw','22434')
        print(f'{Fore.GREEN} Done Unwrapping. {Style.RESET_ALL}')
        Unw2Hgt(f'{DemIdFullPath}.unw',f'{DemIdFullPath}.hgt','800.0')
        MakePol2RecTifFromHgt(f'{DemIdFullPath}.hgt',f'{DemIdFullPath}.hgt.rec.tif')

def StackUnws(MultiTemporalIntStructure,OutputStackFullPath,Width=22434):
    if FileDoesExist('StackTab.txt'): os.remove('StackTab.txt')
    fd = open('StackTab.txt','w')
    for i,x in enumerate(MultiTemporalIntStructure["UnwFullPath"]):
        BaselineInDays = MultiTemporalIntStructure["TemporalBaselineInHours"][i]/(24.0);
        fd.write(f'{x} {BaselineInDays}\n')
    fd.close()
    pg.stacking('StackTab.txt',f'{Width}',f'{OutputStackFullPath}',f'{OutputStackFullPath}_std',f'{OutputStackFullPath}_std_res',f'{params.xCoordPhUnw}',f'{params.yCoordPhUnw}','-','-','1')

def StackInts(IntStructure,StackedFileFullPath,StackedFileRecTifFullPath="None"):
    ## Stacking all Float64 interferogramms in structure using stack_cpx.py
    print(f'Prepare Stacking..')
    if (FileDoesNotExist(StackedFileFullPath)):
        ## Prepare Output
        StackFolder, StackName = os.path.split(StackedFileFullPath)
        MakeFolder(StackFolder)
        if os.path.isfile(f"{StackFolder}/TmpStackingList.txt"):
            subprocess.check_call(f"rm -f {StackFolder}/TmpStackingList.txt",shell="true")
        f = open(f"{StackFolder}/TmpStackingList.txt", "w")
        ## Get Diff_Par (same for all)
        off1 = f'{IntStructure["IntFullPath"][0][0:-4]}.off'
        off2 = f'{IntStructure["IntFullPath"][1][0:-4]}.off'
        pg.create_diff_par(off1,off2,f'{StackFolder}/TmpDiffpar.par',0,0,stdout_flag=False)
        ## Get Tab delimited list for stack_cpx.py
        [ f.write(f'{x} {StackFolder}/TmpDiffpar.par\n') for i, x in enumerate(IntStructure["IntFullPath"])]
        f.close()
        [print(f'{Fore.GREEN}Stack: {x}\n{Style.RESET_ALL}') for i, x in enumerate(IntStructure["IntFullPath"])]
        ## Stack it with hardwired window size for averaging phase offset.
        cmd_stack = f'stack_cpx.py {StackFolder}/TmpStackingList.txt {StackedFileFullPath} tmp_cct -r "{params.xCoordPhUnw} {params.yCoordPhUnw}" -w "50 15" >> GammaProcLogFiles.log'
        subprocess.check_call(cmd_stack,shell="true")
        ## Get Rectified version if you want.
        if (StackedFileRecTifFullPath != "None"):
            MakePol2RecTifFromInt(StackedFileFullPath,StackedFileRecTifFullPath)
        print(f'{Fore.GREEN} Done stacking.  {Style.RESET_ALL}')
    else:
        print(f'{Fore.RED} Stacked outuput file already there. Do nothing. {Style.RESET_ALL}')

def FileDoesNotExist(FileName):
    if os.path.isfile(FileName):
        #print(f'{Fore.RED} {FileName} exists. Do nothing. {Style.RESET_ALL}')
        return False
    else:
        return True

def FileDoesExist(FileName):
    if os.path.isfile(FileName):
        return True
    else:
        return False

def MakePol2RecTifFromInt(IntFullPath,TifFullPath):
    ## Make a Pol2Rec Quicklook using the MasterRec MLI and an incoming FCOMPLEX Float
    if os.path.isfile(params.MasterRecMliFullPath):
        MasterRecFolder, MasterRec = os.path.split(params.MasterRecMliFullPath)
        TifFolder, TifName = os.path.split(TifFullPath)
        MakeFolder(TifFolder)
        cmd_Pol2Rec = f'pol2rec {IntFullPath} {MasterRecFolder}/Master.mli.par tmp.rec  tmp.rec.par  {params.pix_size_p2r} 1 0 - - - - - >>GammaProcLogFiles.log'
        subprocess.check_call(cmd_Pol2Rec, shell=True)
        RecWidth = GetKeywordFromTxtFile(f'{params.MasterRecMliFullPath}.par','range_samples')
        cmd_RasRmgInt = f'rasmph_pwr tmp.rec {params.MasterRecMliFullPath} {RecWidth} - - - 1 1 1  - - {TifFullPath}  >> GammaProcLogFiles.log'
        subprocess.check_call(cmd_RasRmgInt, shell=True)
        subprocess.check_call('rm -f tmp.rec tmp.rec.par',shell=True)
        print(f'{Fore.GREEN} Created quicklook image: {TifFullPath}. {Style.RESET_ALL}')
    else:
        print(f'{Fore.RED} You want a Pol2RecTif but there is no Master rectified MLI. Get it first with: GetMasterRecMli(SlcStructure["SlcFullPath"][0] {Style.RESET_ALL}')

def MakePol2RecFromReal(RealFullPath,RealRecFullPath):
  if os.path.isfile(params.MasterRecMliFullPath):
        MasterRecFolder, MasterRec = os.path.split(params.MasterRecMliFullPath)
        RealFolder, RealName = os.path.split(RealRecFullPath)
        MakeFolder(RealFolder)
        cmd_Pol2Rec = f'pol2rec {RealFullPath} {MasterRecFolder}/Master.mli.par {RealRecFullPath}  {RealRecFullPath}.par  20 0  0 0 - - - - - >>GammaProcLogFiles.log'
        #print(cmd_Pol2Rec)
        subprocess.check_call(cmd_Pol2Rec, shell=True)
        print(f'{Fore.GREEN} Created rec image: {RealRecFullPath}. {Style.RESET_ALL}')
  else:
        print(f'{Fore.RED} You want a Pol2RecTif but there is no Master rectified MLI. Get it first with: GetMasterRecMli(SlcStructure["SlcFullPath"][0] {Style.RESET_ALL}')


def MakePol2RecTifFromReal(IntFullPath,TifFullPath,PhaseScale=10.0):
    ## Make a Pol2Rec Quicklook using the MasterRec MLI and an incoming Float32
    if os.path.isfile(params.MasterRecMliFullPath):
        MasterRecFolder, MasterRec = os.path.split(params.MasterRecMliFullPath)
        TifFolder, TifName = os.path.split(TifFullPath)
        MakeFolder(TifFolder)
        cmd_Pol2Rec = f'pol2rec {IntFullPath} {MasterRecFolder}/Master.mli.par tmp.rec  tmp.rec.par  {params.pix_size_p2r} 0 0 - - - - - >>GammaProcLogFiles.log'
        subprocess.check_call(cmd_Pol2Rec, shell=True)
        RecWidth = GetKeywordFromTxtFile(f'{params.MasterRecMliFullPath}.par','range_samples')
        cmd_RasRmgInt = f'rasrmg tmp.rec {params.MasterRecMliFullPath} {RecWidth} - - - 1 1 {PhaseScale} - - - - {TifFullPath}  >> GammaProcLogFiles.log'
        subprocess.check_call(cmd_RasRmgInt, shell=True)
        subprocess.check_call('rm -f tmp.rec tmp.rec.par',shell=True)
        print(f'{Fore.GREEN} Created quicklook image: {TifFullPath}. {Style.RESET_ALL}')
    else:
        print(f'{Fore.RED} You want a Pol2RecTif but there is no Master rectified MLI. Get it first with: GetMasterRecMli(SlcStructure["SlcFullPath"][0] {Style.RESET_ALL}')

def MakePol2RecTifFromHgt(IntFullPath,TifFullPath,MetersPerColorCycle=400):
    ## Make a Pol2Rec Quicklook using the MasterRec MLI and an incoming Float32
    if os.path.isfile(params.MasterRecMliFullPath):
        MasterRecFolder, MasterRec = os.path.split(params.MasterRecMliFullPath)
        TifFolder, TifName = os.path.split(TifFullPath)
        MakeFolder(TifFolder)
        cmd_Pol2Rec = f'pol2rec {IntFullPath} {MasterRecFolder}/Master.mli.par tmp.rec  tmp.rec.par  {params.pix_size_p2r} 0 0 - - - - - >>GammaProcLogFiles.log'
        subprocess.check_call(cmd_Pol2Rec, shell=True)
        RecWidth = GetKeywordFromTxtFile(f'{params.MasterRecMliFullPath}.par','range_samples')
        cmd_RasRmgInt = f'rashgt tmp.rec {params.MasterRecMliFullPath} {RecWidth} - - - 1 1 {MetersPerColorCycle} - - - {TifFullPath}  >> GammaProcLogFiles.log'
        cmd_RasRmgInt = f'visdt_pwr.py tmp.rec {params.MasterRecMliFullPath} {int(RecWidth)} -80.0 80.0 -p {TifFullPath}  >> GammaProcLogFiles.log'
        print(cmd_RasRmgInt)
        subprocess.check_call(cmd_RasRmgInt, shell=True)

        #subprocess.check_call('rm -f tmp.rec tmp.rec.par',shell=True)
        print(f'{Fore.GREEN} Created quicklook image: {TifFullPath}. {Style.RESET_ALL}')
    else:
        print(f'{Fore.RED} You want a Pol2RecTif but there is no Master rectified MLI. Get it first with: GetMasterRecMli(SlcStructure["SlcFullPath"][0] {Style.RESET_ALL}')

def CreateFloat32Array(SrcDirectory,Suffix,Width,TrgtDirectory,Scale=1.0):
    #Assemble Float32s in matrix and scale them while doing so.
    #Todo: Function too long. Outsource resampling.
    NewWidth=int(np.round(Width*Scale)+1)

    if os.path.isfile(f'{TrgtDirectory}/Array.pickle'):
        print(f'{Fore.GREEN} Loading Data from Previous Scaling: {TrgtDirectory}/Array.pickle {Style.RESET_ALL}')
        Array = LoadPickle(f'{TrgtDirectory}/Array.pickle')
    else:
        print(f'{Fore.GREEN} Parsing Directory: {SrcDirectory} for {Suffix} Files. {Style.RESET_ALL}')
        FileList = glob.glob(f'{SrcDirectory}/*{Suffix}')
        print(f'{Fore.GREEN} Found {len(FileList)} Files.{Style.RESET_ALL}')
        ScaledFolder=f'{TrgtDirectory}/ScaledFiles/'
        MakeFolder(ScaledFolder)
        for i,x in enumerate(FileList):
            RealFolder, RealName = os.path.split(x)
            ResampledFileNameFullPath = f'{TrgtDirectory}/ScaledFiles/{RealName[:-len(Suffix)]}.scaled'
            if os.path.isfile(ResampledFileNameFullPath):
                print(f'{Fore.GREEN} Already resampled. Do nothing. {Style.RESET_ALL}')
            else:
                cmd_Resample=f'resamp_image {x} {Width} {Scale} {Scale} {ResampledFileNameFullPath} - - >>GammaProcLogFiles.log'
                subprocess.check_call(cmd_Resample,shell=True)
        FileList = glob.glob(f'{TrgtDirectory}/ScaledFiles/*scaled')
        print(f'{Fore.RED} Done scaling files. New width probably {NewWidth}, but better double check in output logs.{Style.RESET_ALL}')
        Arr1 = np.fromfile(FileList[0],dtype='>f4')
        NewHeight = int(len(Arr1)/NewWidth)
        Array = np.zeros((NewWidth*NewHeight,int(len(FileList)))).astype('>f4')
        for i,k in enumerate(FileList):
            Array[:,i] = np.fromfile(FileList[i],dtype='>f4').transpose()
        SavePickle(Array,f'{TrgtDirectory}/Array.pickle')

    return Array, NewWidth

def PCA_RD(Array,TrgtDirectoryOutput):
    #Do a PCA of the incoming array.
    ## Todo, Function too long. Outsource plotting of images.
    PickleNameEVs = f'{TrgtDirectoryOutput}/EigenValues.pickle'
    PickleNamePCs = f'{TrgtDirectoryOutput}/PrincipalComponents.pickle'
    if os.path.isfile(PickleNameEVs):
        print(f'{Fore.GREEN} Loading results from previous run.{Style.RESET_ALL}')
        EigenValues = LoadPickle(PickleNameEVs)
        PrincipalComponents= LoadPickle(PickleNamePCs)
    else:
        print(f'{Fore.RED} PCA without standardizing ongoing...{Style.RESET_ALL}')
        CovMatrix = np.cov(np.transpose(Array))
        EigenValues,EigenVectors = np.linalg.eig(CovMatrix)
        PrincipalComponents = np.dot(Array, EigenVectors)
        ## EigenValues are returned already sorted.
        SavePickle(EigenValues,PickleNameEVs)
        SavePickle(PrincipalComponents,PickleNamePCs)
    return EigenValues,PrincipalComponents

def PlotPCA(PrincipalComponents,TrgtDirectoryOutput,NewWidth,InvScale=1.0,PrintPCs=4):
    MasterRecSlcFolder, SlcName = os.path.split(params.MasterRecMliFullPath)
    for k in np.arange(0,PrintPCs):
        Arr = PrincipalComponents[:,k].astype('>f4')
        Arr.tofile('pc.int',format='>f4')
        Width = GetKeywordFromTxtFile(f'{MasterRecSlcFolder}/Master.mli.par','range_samples')
        Height = GetKeywordFromTxtFile(f'{MasterRecSlcFolder}/Master.mli.par','azimuth_lines')
        cmd_Resample = f'resamp_image pc.int {NewWidth} {InvScale} {InvScale} pc.scaled.int {Width} {Height} >> GammaProcLogFiles.log'
        subprocess.check_call(cmd_Resample,shell=True)
        MakePol2RecTifFromReal('pc.scaled.int',f'{TrgtDirectoryOutput}/PC{k}.rec.tif')
        subprocess.check_call('rm -f pc.int pc.scaled.int',shell=True)

def SavePickle(Array,FileNameFullPath):
    PickleOut = open(FileNameFullPath,'wb')
    pickle.dump(Array,PickleOut)
    PickleOut.close()

def LoadPickle(FileNameFullPath):
    PickleIn = open(FileNameFullPath,'rb')
    ArrayName=pickle.load(PickleIn)
    PickleIn.close()
    return ArrayName

def GetKeywordFromTxtFile(filename,keyword):
    #Assuming that parameter is seperated with : after keyword
    with open(filename) as f:
        for line in f:
            if line.startswith(keyword):
                try:
                    return float(line.split(':')[1])
                except ValueError:
                    print('{Something awful has happened while reading {keyword} from {filename}')

def GetProcessingParametersFromFile(filename):
    #from importlib.machinery import SourceFileLoader
    if os.path.isfile(filename):
        global params
        f = open(filename)
        params = imp.load_source('params',filename, f)
        #params = SourcefileLoader(f)
        print(f'{Fore.GREEN} Loaded parameters from {filename}. They are now available via params.XX. {Style.RESET_ALL}')
        f.close()
    else:
        print(f'{Fore.RED} I cannot find: {filename}! {Style.RESET_ALL}')
def convert_bytes(num):
    """
    this function will convert bytes to MB.... GB... etc
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0

def file_size(file_path):
    """
    this function will return the file size
    """
    if os.path.isfile(file_path):
        file_info = os.stat(file_path)
        return convert_bytes(file_info.st_size)

def CheckIfFileIsDone(FilePath,FileSize):
    if FileDoesNotExist(FilePath):
        return False
    else:
        sn = int(''.join(filter(str.isdigit, file_size(FilePath))))
        if (sn < FileSize):
            Done = False
        else:
            Done = True
    return Done

def AuxStringInFile(String,TxtFile):
    if (FileDoesExist(TxtFile)):
        if String in open(TxtFile).read():
            print('Found String on Blacklist.')
            return True
        else:
            return False
        #with open(TxtFile) as f:
        #    cstring = f.read()
        #    cstring = cstring.rstrip()
        #    print(f'{cstring} {String}')
        #    if (String == cstring):
        #        print('Hallo!!!!!')
        #        return True
        #    else:
        #        return False
    else:
        return False

def runningMeanFast(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]

def MakeFolder(Folder):
    try:
        os.stat(Folder)
        print(f'Exists: {Folder}.')
    except:
        os.makedirs(Folder)
        print(f'Created: {Folder}.')

def linearf(B, x):
    '''Linear function y = m*x + b'''
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x + B[1]

def GetDecimalDayOfYear(DTime):
    DayOfYear = DTime.timetuple().tm_yday
    hour = DTime.hour;minute = DTime.minute; secs = DTime.second
    DecimalDay = (hour+minute/60.0 + secs/3600.0)/24.0
    DecimalDayOfYear=DayOfYear+DecimalDay
    return DecimalDayOfYear
