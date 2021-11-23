################# General Config ##########################

############### phase unwrapping stationary target coordinates #################################

xCoordPhUnw=10694
yCoordPhUnw=141
UnwMask='AuxFiles/Mask/lmask.bmp'
HgtMask='AuxFunctions/Masks/MaskHgt.bmp'
## for pretty background in rectified images.
#MasterRecMliFullPath="/esd/esd02/data/rdrews/data_large/Antarctica_K050_2018/PROC/GPRI/12HTripleDifferences/MasterRecMli/MasterRec.mli"
#MasterRecMliFullPath="/esd/esd02/data/rdrews/data_large/Antarctica_K050_2018/PROC/GPRI/12HTripleDifferences/MasterRecMli/MasterRec.mli"
MasterRecMliFullPath="/esd/esd01/docs/rdrews/current/GPRI_Lauterbrunnen/proc2/MasterRec/Master.rec.mli"

# options for the gamma functions used
############### gamma_int_semiauto ########################

#multi_look
rlks=1     # number of range looks
azlks=1    # number of azimuth looks
loff=0     # offset to starting line (default: 0)
nlines="DEFAULT"   # number of SLC lines to process (enter - for default: entire file)
scale=1    # scale factor for output MLI (default: 1.0)
exp=1      # exponent for the output MLI (default: 1.0)

# rasmph_pwr parameters
start_cpx=1  # starting line of cpx (enter - for default: 1)
start_pwr=1  # starting line of pwr (enter - for default: 1)
nlines=0     # number of lines to display (enter - or 0 for default: to end of file)
pixavr=5     # number of pixels to average in range (enter - for default: 1)
pixavaz=1    # number of pixels to average in azimuth (enter - for default: 1)
scale=.5     # pwr display scale factor (enter - for default: 1.0)
exp=.35      # pwr display exponent (enter - for default: 0.35)
LR=1         # left/right flipping flag (enter - for default: 1: normal (default), -1: mirror image)

# adf parameters 0.4 32 3.
alpha=0.4          # exponent for non-linear filtering (enter - for default: 0.40)
nfft=32            # filtering FFT window size, 2**N, 8 --> 512, (enter - for default: 32)
cc_win=3.          # coherence parameter estimation window size odd, max: 15 (enter - for default: 5)
step="DEFAULT"        # processing step (enter - for default: nfft/8)
loff=0             # offset to starting line to process (enter - for default: 0)
nlines_adf=0       # number of lines to process (enter - for default: to end of file)
wfrac=0.700        # minimum fraction of points required to be non-zero in the filter window (enter - for default: 0.200)

# SLC_intf parameters
rlks_intf=1      # number of range looks
azlks_intf=1     # number of azimuth looks
loff_intf=0      # offset to starting line relative to SLC-1 for interferogram (default=0)
nlines_intf="DEFAULT"    # number of SLC lines to process (enter - for default: to end of file)
sps_flg=0        # range spectral shift flag:
                 #   1: apply range spectral shift filter (default)
                 #   0: do not apply range spectral shift filter
azf_flg=0        # azimuth common band filter flag:
                 #   1: apply azimuth common-band filter (default)
                 #   0: do not apply azimuth common band filter
rp1_flg=1        # SLC-1 range phase mode
                 #   0: nearest approach (zero-Doppler) phase
                 #   1: ref. function center (Doppler centroid) (default)
rp2_flg=1        # SLC-2 range phase mode
                 #   0: nearest approach (zero-Doppler) phase
                 #   1: ref. function center (Doppler centroid) (default)


# create_offset parameters
algorithm=1   # algorithm  offset estimation algorithm
              #              1: intensity cross-correlation (default)
              #              2: fringe visibility
rlks_co=1     # rlks       number of interferogram range looks (enter -  for default: 1)
azlks_co=1    # azlks      number of interferogram azimuth looks (enter - for default: 1)
iflg=0        # iflg       interactive mode flag (enter -  for default)
              #              0: non-interactive
              #              1: interactive (default)

# offset_fit parameters
coffs=7.            # (output) culled range and azimuth offset estimates (fcomplex, enter - for none)
coffsets=1          # (output) culled offset estimates and cross-correlation values (text format, enter - for none)
thresh="DEFAULT"			# cross-correlation threshold (enter - for default from OFF_par)
npoly=4             # number of model polynomial parameters (enter - for default, 1, 3, 4, 6, default: 4)
interact_flag="off" # interactive culling of input data:
                    #    0: off (default)
                    #    1: on

# Polar to rectangular conversion for GPRI SLC and MLI image data
# usage: pol2rec <data1> <SLC_par1> <data2> <SLC_par2> <pix_size> <type> [mode] [xmin] [nx] [ymin] [ny] [rmax]
# input parameters:
pix_size_p2r=2 # (output) output pixel size (meters)
type_fcomplex_p2r_1=1    # input data type
type_float_p2r_2=0
mode_p2r_1=1   	 # interpolation algorithm
mode_p2r_2=1
xmin_p2r_1="DEFAULT"   # starting x coordinate (enter - for default: calculated from image)
xmin_p2r_2="DEFAULT"   # starting x coordinate (enter - for default: calculated from image)
nx_p2r_1="DEFAULT"     # number of x samples in the output image (enter - for default: calculated from image)
nx_p2r_2="DEFAULT"
ymin_p2r_1="DEFAULT"   # starting y coordinate (enter - for default: calculated from image)
ymin_p2r_2="DEFAULT"
ny_p2r_1="DEFAULT"     # number of y samples in the output image (enter - for default: calculated from image)
ny_p2r_2="DEFAULT"
rmax_p2r_1="DEFAULT"   # maximum slant range in the GPRI image to resample (enter - for default: maximum slant range of the input image)
rmax_p2r_2="DEFAULT"

# rasmph_pwr parameters (currently defaults) for rectified image
start_cpx_p2r=1  # starting line of cpx (enter - for default: 1)
start_pwr_p2r=1  # starting line of pwr (enter - for default: 1)
nlines_p2r=0     # number of lines to display (enter - or 0 for default: to end of file)
pixavr_p2r=5     # number of pixels to average in range (enter - for default: 1)
pixavaz_p2r=1    # number of pixels to average in azimuth (enter - for default: 1)
scale_p2r=.5     # pwr display scale factor (enter - for default: 1.0)
exp_p2r=.35      # pwr display exponent (enter - for default: 0.35)
LR_p2r=1         # left/right flipping flag (enter - for default: 1: normal (default), -1: mirror image)

####################### gamma_baseline.sh #########################
base_mflag=3	# mflag     baseline estimation method flag (enter - for default)
				# mflag        b_para    b_perp    input
					# 0:         orbits    orbits    p1,p2  (default)
					# 1:         offsets   offsets   p1,p2,off
					# 2:         orbits    fft       p1,p2,off,int
					# 3:         offsets   fft       p1,p2,off,int
					# 4:         fft       fft       p1,off,int

# one of these options segfaults the base_init program
#base_nrfft=-           # size of range FFT   (512, 1024,...) (enter - for default determined from image width)
#base_nazfft=-          # size of azimuth FFT (512, 1024,...) (enter - for default determined from image azimuth lines)
#base_r_samp=center     # range pixel offset to center of the FFT window (enter - for default, default: range center)
#base_az_line=center    # line offset from start of the interf. for the  FFT window (enter - for default, default: azimuth center)

base_time_rev=1 # ?

  # NOTE: Not all  input data files are required for the different methods

        # enter - for files that are not provided


####################### gamma_cpx2real.sh #########################
# the type of information to extract from the complex file
# cpx_to_real

# output data type
#         0: real part
#         1: imaginary part
#         2: intensity (re*re + im*im)
#         3: magnitude (sqrt(re*re + im*im))
#         4: phase (atan2(im, re))

setType_ctr=4
setType_rtc=1

###################### gamma_pol2rec.sh ######################

# pol2rec
rec_pix_size=2. # (output) output pixel size (meters)
rec_dtype=0    	# input data type
rec_mode=1   	# interpolation algorithm
rec_xmin="DEFAULT"   	# starting x coordinate (enter - for default: calculated from image)
rec_nx="DEFAULT"     	# number of x samples in the output image (enter - for default: calculated from image)
rec_ymin="DEFAULT"   	# starting y coordinate (enter "DEFAULT" for default: calculated from image)
rec_ny="DEFAULT"     	# number of y samples in the output image (enter "DEFAULT" for default: calculated from image)i
rec_rmax="DEFAULT"   	# maximum slant range in the GPRI image to resample (enter "DEFAULT" for default: maximum slant range of the input image)
				# NOTE: center image line of the scan defines the direction of the X axis

# rasmph_pwr parameters (currently defaults)
start_cpx_pol2r=1    # starting line of cpx (enter - for default: 1)
start_pwr_pol2r=1    # starting line of pwr (enter - for default: 1)
nlines_pol2r=0       # number of lines to display (enter - or 0 for default: to end of file)
pixavr_pol2r=1       # number of pixels to average in range (enter - for default: 1)
pixavaz_pol2r=1      # number of pixels to average in azimuth (enter - for default: 1)
scale_pol2r=.5       # pwr display scale factor (enter - for default: 1.0)
exp_pol2r=.35        # pwr display exponent (enter - for default: 0.35)
LR_pol2r=1           # left/right flipping flag (enter - for default: 1: normal (default), -1: mirror image)

# rasrmg
start_unw_pol2r=1 	#  starting line of unwrapped phase file (default: 1)
start_pwr_pol2r=1 	#  starting line of intensity file (default: 1)
nlines_pol2r=0    	#  number of lines to display (default 0: to end of file)
pixavr_pol2r=5    	#  number of pixels to average in range (default: 1)
pixavaz_pol2r=1   	#  number of pixels to average in azimuth (default: 1)
ph_scale_pol2r=0.05 #  phase display scale factor (default:  0.33333) NOTE: one color cycle: 2PI/ph_scale
scale_pol2r=1       #  pwr display scale factor (default: 1.0)
exp_pol2r=.35       #  pwr display exponent (default: .35)
ph_offset_pol2r=0.0 #  phase offset in radians subtracted from unw (default: 0.0)
LR_pol2r=1          #  left/right mirror image flag, (1: normal (default), -1: mirror image)
					#  (output) image filename, extension determines the format, enter - for default: *.ras
                    #    *.bmp BMP format; *.ras Sun raster format; *.tif TIFF format

#cc
start_cc_pol2r=1 	#  starting line of threshold data file (default=1)
cc_min_pol2r=.2   	#  pixels with cc values below cc_min are displayed using greyscale (default=.2)

#rashgt
start_hgt_pol2r=1 	# starting line of hgt (enter - for default: 1)
start_pwr_pol2r=1 	# starting line of pwr (enter - for default: 1)
nlines_pol2r=0 		# number of lines to display (enter - or 0 for default: to end of file)
pixavr_pol2r=1 		# number of pixels to average in range (enter - for default: 1)
pixavaz_pol2r=1 	# number of pixels to average in azimuth (enter - for default: 1)
mcycle_pol2r=0.3 	# meters per color cycle (enter - for default: 160.0)
scale_pol2r=1 		# display scale factor (enter - for default: 1.0)
exp_pol2r=0.35 		# display exponent (enter - for default: 0.35)

################## gamma_phunw.sh ######################
## ph_slope_base
#int_type_phunw_1=1
#inverse_phunw_1=0
#int_type_phunw_2=0
#inverse_phunw_2=1

## adf
alpha_phunw=0.40   # exponent for non-linear filtering (enter - for default: 0.40)
nfft_phunw=32      # filtering FFT window size, 2**N, 8 --> 512, (enter - for default: 32)
cc_win_phunw=5     # coherence parameter estimation window size odd, max: 15 (enter - for default: 5)
step_phunw="DEFAULT"  # processing step (enter - for default: nfft/8)
loff_phunw="DEFAULT"       # offset to starting line to process (enter - for default: 0)
nlines_phunw="DEFAULT"     # number of lines to process (enter - for default: to end of file)
wfrac_phunw=0.2    # minimum fraction of points required to be non-zero in the filter window (enter - for default: 0.200)

## cc_wave
bx_cc_wave=5 		# estimation window size in columns (enter - for default:5.0)
by_cc_wave=5		# estimation window size in lines (enter - for default:5.0)
wflg_cc_wave=0      # estimation window (enter - for default):
                    #  0: rectangular (default)
                    #  1: triangular
                    #  2: Gaussian
                    #  3: normalized vector sum with rectangular window
                    #     NOTE: This estimator does not use the MLI data, even when specified
#xmin_cc_wave=		# starting range pixel offset (enter - for default: 0)
#xmax_cc_wave=  	# last range pixel offset (enter - for default: width - 1)
#ymin_cc_wave=  	# starting azimuth row offset, relative to start (enter -  for default: 0)
#ymax_cc_wave=  	# last azimuth row offset, relative to start (enter - for default: nlines - 1)


## rascc_mask
#input parameters:
#cc         (input)interferometric correlation image (float)
#  pwr        (input)intensity image (float, enter - if not available)
#  width      number of samples/row

start_cc_phunw=1  	# starting line of coherence image (default: 1)
start_pwr_phunw=1 	# starting line of intensity image (default: 1)
nlines_phunw=0    	# number of lines to display (default=0: to end of file)
pixavr_phunw=1    	# number of pixels to average in range (default: 1)
pixavaz_phunw=1   	# number of pixels to average in azimuth (default: 1)
cc_thres_phunw=0.4  # coherence threshold for masking, pixels with cc < cc_thres are set to 0 (default: 0.0)
pwr_thres_phunw=0. 	# intensity threshold, pixels with relative intensity below pwr_thres are set to 0 (default: 0.)
cc_min_phunw=0.4    # minimum coherence value used for color display (default: 0.1)
cc_max_phunw=0.99   # maximum coherence value used for color display (default: 0.9)
scale_phunw=1.     	# intensity display scale factor (default: 1.)
exp_phunw=.35       # intensity display exponent (default: .35)
LR_phunw=1        	# left/right mirror image flag, (1: normal (default), -1: mirror image)
#rasf_phunw=      	# (output) image filename, extension determines the format, enter - for default: *.ras
            # *.bmp BMP format
            # *.ras Sun raster format
            # *.tif TIFF format

## rascc_mask_thinning
nmax_phunw=3           # number of sampling reduction runs (default: 3)
thresh_1_phunw=0.3     # first threshold (used for smallest scale sampling reduction)
thresh_2_phunw=0.5     # further thresholds
thresh_3_phunw=0.7     # threshold nmax (used for largest scale sampling reduction)

## mcf
tri_mode_mcf=0          # triangulation mode
                        #    0: filled triangular mesh (default)
                        #    1: Delaunay triangulation
roff_mcf="DEFAULT"              #  offset to starting range of section to unwrap (default: 0)
loff_mcf="DEFAULT"              #  offset to starting line of section to unwrap (default: 0)
nr_mcf="DEFAULT"             #  number of range samples of section to unwrap (default("DEFAULT"): width "DEFAULT" roff)
nlines_mcf="DEFAULT"          #  number of lines of section to unwrap (default("DEFAULT"): total number of lines "DEFAULT" loff)
npat_r_mcf=1            #  number of patches in range
npat_az_mcf=1           #  number of patches in azimuth
ovrlap_mcf="DEFAULT"            #  overlap between patches in pixels (overlap >= 7, default("DEFAULT"): 512)
r_init_mcf= 10694        #  phase reference point range offset (default("DEFAULT"): roff)
az_init_mcf= 141        #  phase reference point azimuth offset (default("DEFAULT"): loff)
init_flag_mcf=1         #  flag to set phase at reference point
                        #    0: use initial point phase value (default)
                        #    1: set phase to 0.0 at initial point

## rasrmg
start_unw_rasrmg=1 		#  starting line of unwrapped phase file (default: 1)
start_pwr_rasrmg=1 		#  starting line of intensity file (default: 1)
nlines_rasrmg=0    		#  number of lines to display (default 0: to end of file)
pixavr_rasrmg=1    		#  number of pixels to average in range (default: 1)
pixavaz_rasrmg=1   		#  number of pixels to average in azimuth (default: 1)
ph_scale_rasrmg=1.2 	#  phase display scale factor (default:  0.33333) NOTE: one color cycle: 2PI/ph_scale
scale_rasrmg=1      	#  pwr display scale factor (default: 1.0)
exp_rasrmg=.35      	#  pwr display exponent (default: .35)
ph_offset_rasrmg=0.0 	#  phase offset in radians subtracted from unw (default: 0.0)
LR_rasrmg=1          	#  left/right mirror image flag, (1: normal (default), -1: mirror image)
						#  (output) image filename, extension determines the format, enter - for default: *.ras
						#    *.bmp BMP format; *.ras Sun raster format; *.tif TIFF format
#cc       				#  (input) correlation threshold data file (float): e.g. correlation
start_cc_rasrmg=1 		#  starting line of threshold data file (default=1)
cc_min_rasrmg=.2   		#  pixels with cc values below cc_min are displayed using greyscale (default=.2)

########################### gamma_displ.sh ###################

# dispmap
mode_displ=0 # flag indicating displacement mode:
             # 0: along light of sight (LOS) [m] (+: towards sensor) (default)
             # 1: vertical displacement [m] (+: increasing height)
             # 2: horizontal displacement [m] (+: decreasing ground range)
sflg_displ=0 # sign flag:
             # 0: do not change sign (default)
             # 1: multiply displacement by -1

# rashgt
start_hgt_displ=1    # starting line of hgt (enter - for default: 1)
start_pwr_displ=1    # starting line of pwr (enter - for default: 1)
nlines_displ="DEFAULT"       # number of lines to display (enter - or 0 for default: to end of file)
pixavr_displ=1       # number of pixels to average in range (enter - for default: 1)
pixavaz_displ=1      # number of pixels to average in azimuth (enter - for default: 1)
m_cycle_displ=0.3    # meters per color cycle (enter - for default: 160.0)
scale_displ=1        # display scale factor (enter - for default: 1.0)
exp_displ=0.35       # display exponent (enter - for default: 0.35)
LR_displ=1           # left/right mirror image flag, (enter - for default: 1: normal (default), -1: mirror image)

############################ gamma_stacking.sh ####################
#usage: stacking <DIFF_tab> <width> <ph_rate> <sig_ph_rate> <sig_ph> <roff> <loff> [nr] [nl] [np_min] [tscale]

#roff_stack=1941       	# range pixel offset to center of the phase reference region
roff_stack=10694  # range pixel offset to center of the phase reference region
#loff_stack=641      	# line offset to center of the phase reference region
loff_stack=141 # line offset to center of the phase reference region
nr_stack=16          	# number of range pixels to average in the phase reference region (enter - for default:16)
nl_stack=16         	# number of lines average in the phase reference region (enter - for default: 16)
np_min_stack="DEFAULT"     		# min. number of phase values required to accept phase rate estimate (enter - for default: all files)

## rasrmg
start_unw_stack=1 	#  starting line of unwrapped phase file (default: 1)
start_pwr_stack=1 	#  starting line of intensity file (default: 1)
nlines_stack=0    	#  number of lines to display (default 0: to end of file)
pixavr_stack=1    	#  number of pixels to average in range (default: 1)
pixavaz_stack=1   	#  number of pixels to average in azimuth (default: 1)
ph_scale_stack=0.05 #  phase display scale factor (default:  0.33333) NOTE: one color cycle: 2PI/ph_scale
scale_stack=1       #  pwr display scale factor (default: 1.0)
exp_stack=.35       #  pwr display exponent (default: .35)
ph_offset_stack=0.0 #  phase offset in radians subtracted from unw (default: 0.0)
LR_stack=1          #  left/right mirror image flag, (1: normal (default), -1: mirror image)
