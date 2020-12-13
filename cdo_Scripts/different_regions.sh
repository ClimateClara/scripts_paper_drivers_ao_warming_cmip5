#!/bin/bash

#USE MASKS OF DIFFERENT REGIONS TO DEEPEN ANALYSIS

################# DECLARE THE PATHS #############################3
path1=/work/mh0033/m300411/DataEB/RAW_DATA
path2=/work/mh0033/m300411/DataEB/WORK_DATA
###################################################################

#resolution of ocean in the different models
res_ACCESS1_0='360x300'
res_ACCESS1_3='360x300'
res_CanESM2='128x64'
res_CCSM4='320x384'
res_CESM1_BGC='320x384'
res_CESM1_CAM5='320x384'
res_CMCC_CM='182x149'
res_CMCC_CMS='182x149'
res_CNRM_CM5='362x292'
res_CSIRO_Mk3_6_0='192x189'
res_FGOALS_g2='360x196'
res_GFDL_CM3='360x200'
res_GFDL_ESM2G='360x210'
res_GFDL_ESM2M='360x200'
res_GISS_E2_R='288x180'
res_HadGEM2_CC='360x216'
res_HadGEM2_ES='360x216'
res_IPSL_CM5A_LR='182x149'
res_IPSL_CM5A_MR='182x149'
res_IPSL_CM5B_LR='182x149'
res_MIROC5='256x224'
res_MPI_ESM_LR='256x220'
res_MPI_ESM_MR='802x404'
res_MRI_CGCM3='360x368'
res_NorESM1_M='320x384'
res_NorESM1_ME='320x384'

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do

mod2=`echo res_$mod | sed -e 's/-/_/g'`
echo $mod2

###### grid to work with
cdo -f nc -sellonlatbox,0,360,66,90 -random,r"${!mod2}" $path2/$mod/"$mod"_grid.nc


##### to look at temperature gradients
#Low latitudes
cdo -f nc -sellonlatbox,0,360,50,66 -random,r"${!mod2}"   $path2/$mod/"$mod"_grid_lowlat.nc
#North Atlantic
cdo -f nc -sellonlatbox,315,15,45,66 -random,r"${!mod2}" $path2/AOarea2/"$mod"_NAtl.nc
###########

######## to look at surface pressure
#Barents Sea 
cdo -f nc -sellonlatbox,17.5,60,66,80 -random,r"${!mod2}" $path2/AOarea2/"$mod"_BarentsSeaDomain.nc
###########

######## to look at temperature
#Fram Strait 
cdo -f nc -sellonlatbox,348.5,10.5,79.5,81.5 -random,r"${!mod2}" $path2/$mod/"$mod"_FS.nc


#Canadian Archipelago
cdo -f nc -sellonlatbox,231.5,300.5,70.5,82 -random,r"${!mod2}" $path2/$mod/"$mod"_CA.nc

#Barents Opening
cdo -f nc -sellonlatbox,16.5,19,70,76.5 -random,r"${!mod2}" $path2/$mod/"$mod"_BO.nc

#Bering Strait
cdo -f nc -sellonlatbox,189,194,65,66 -random,r"${!mod2}" $path2/$mod/"$mod"_BS.nc

#Denmark Strait
cdo -f nc -sellonlatbox,323,337.5,65.5,66.5 -random,r"${!mod2}" $path2/$mod/"$mod"_DS.nc

#Iceland - Faroe
cdo -f nc -sellonlatbox,346.4,352.6,62.2,64.9 -random,r"${!mod2}" $path2/$mod/"$mod"_IF.nc

#Faroe - Scotland
cdo -f nc -sellonlatbox,353.1,355,58.7,62 -random,r"${!mod2}" $path2/$mod/"$mod"_FaSc.nc

#Fram Strait right
cdo -f nc -sellonlatbox,0,15,77,80 -random,r"${!mod2}" $path2/$mod/"$mod"_FS_r.nc

#Fram Strait left
cdo -f nc -sellonlatbox,330,0,77,80 -random,r"${!mod2}" $path2/$mod/"$mod"_FS_l.nc

#Subpolar region
cdo -f nc -sellonlatbox,0,360,50,77 -random,r"${!mod2}" $path2/$mod/"$mod"_subpolar.nc

##############

###### for Arctic Amplification
#cdo -f nc -sellonlatbox,0,360,-30,30 -random,r"${!mod2}" $path2/AOarea2/"$mod"_Tropics.nc

############################

done




