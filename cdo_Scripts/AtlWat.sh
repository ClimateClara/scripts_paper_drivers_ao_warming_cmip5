#!/bin/bash

#CUT OUT THE UPPER 1000M TO LOOK AT ATLANTIC WATER

################# DECLARE THE PATHS #############################3
path1=/work/mh0033/m300411/DataEB/RAW_DATA
path2=/work/mh0033/m300411/DataEB/WORK_DATA
###################################################################

#indexes for levels englobing the upper 1000m of ocean in the different models
lev_ACCESS1_0='1/30'
lev_ACCESS1_3='1/30'
lev_CanESM2='1/27'
lev_CCSM4='1/40'
lev_CESM1_BGC='1/40'
lev_CESM1_CAM5='1/40'
lev_CMCC_CM='1/22'
lev_CMCC_CMS='1/22'
lev_CNRM_CM5='1/26'
lev_CSIRO_Mk3_6_0='1/20'
lev_FGOALS_g2='1/22'
lev_GFDL_CM3='1/35'
lev_GFDL_ESM2G='1/35'
lev_GFDL_ESM2M='1/35'
lev_GISS_E2_R='1/13'
lev_HadGEM2_CC='1/26'
lev_HadGEM2_ES='1/26'
lev_IPSL_CM5A_LR='1/22'
lev_IPSL_CM5A_MR='1/22'
lev_IPSL_CM5B_LR='1/22'
lev_MIROC5='1/30'
lev_MPI_ESM_LR='1/24'
lev_MPI_ESM_MR='1/24'
lev_MRI_CGCM3='1/31'
lev_NorESM1_M='1/37'
lev_NorESM1_ME='1/37'

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do

mod2=`echo lev_$mod | sed -e 's/-/_/g'`
echo $mod2

###### vertical mean over the upper 1000 m

for strait in {FS_l,FS_r}
do

echo 'Remap map' $mod $strait
cdo remapcon,$path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/"$mod"_"$strait".nc $path2/AOarea2/"$strait"_ocean_grid_"$mod".nc
echo 'cut out strait' $mod $strait
cdo ifthen $path2/AOarea2/"$strait"_ocean_grid_"$mod".nc -selvar,thetao $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_"$strait"_"$mod"_ensmean_186101-209912.nc

echo 'vertmean over upper 1000 m' $mod $strait
cdo vertmean -setctomiss,0 -sellevidx,"${!mod2}" $path2/$mod/hist_rcp2/thetao_"$strait"_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_vertmean_1000_"$strait"_"$mod"_ensmean_186101-209912.nc
echo 'fldmean' $mod $strait
cdo fldmean $path2/$mod/Method170516/thetao_vertmean_1000_"$strait"_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_fldmean_"$strait"_"$mod"_ensmean_186101-209912.nc
done
done




