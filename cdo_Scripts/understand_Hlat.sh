#!/bin/bash

#LOOK AT DIFFERENT VARIABLES TO UNDERSTAND CHANGES IN THE MERIDIONAL OCEANIC HEAT FLUX

################# DECLARE THE PATHS #############################
path1=/work/mh0033/m300411/DataEB/RAW_DATA
path2=/work/mh0033/m300411/DataEB/WORK_DATA
###################################################################

####### CHANGES IN ICE PRODUCTION

#substract September sea-ice volume (minimum) from next year March sea-ice volume (maximum sea ice) to have the ice production
for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,NorESM1-M,NorESM1-ME,MRI-CGCM3,MPI-ESM-LR,MPI-ESM-MR}
do
cdo sub -selmon,03 -selyear,1862/2099 $path2/$mod/Method170516/totsiv_Arctic_fldmean_"$mod"_ensmean_186101-209912.nc -selmon,09 -selyear,1861/2098 $path2/$mod/Method170516/totsiv_Arctic_fldmean_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/iceprod0903_Arctic_fldmean_"$mod"_ensmean_186101-209912.nc
done

###################################

########### COMPUTE GLOBAL MEAN AIR TEMPERATURE AND ARCTIC TEMPERATURE (OCEAN+LAND) TO USE FOR COMPUTATION OF ARCTIC AMPLIFICATION

#for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,NorESM1-M,NorESM1-ME,MRI-CGCM3,MPI-ESM-LR,MPI-ESM-MR}
#do
#cdo fldmean $path2/$mod/hist_rcp2/tas_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/AA_"$mod"_ensmean_186101-209912.nc
#done

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,NorESM1-M,NorESM1-ME,MRI-CGCM3,MPI-ESM-LR,MPI-ESM-MR}
do
cdo fldmean -sellonlatbox,0,360,-30,30 $path2/$mod/hist_rcp2/tas_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/AAtrop_"$mod"_ensmean_186101-209912.nc
done

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,NorESM1-M,NorESM1-ME,MRI-CGCM3,MPI-ESM-LR,MPI-ESM-MR}
do
cdo fldmean -sellonlatbox,0,360,66,90 $path2/$mod/hist_rcp2/tas_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/AAarc_"$mod"_ensmean_186101-209912.nc
done


##############################

############# COMPUTE MEAN AIR TEMPERATURE OVER OCEAN OF LOWER LATITUDES (50 TO 66 DEGREES NORTH)

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,NorESM1-M,NorESM1-ME,MRI-CGCM3,MPI-ESM-LR,MPI-ESM-MR}
do
cdo fldmean -sellonlatbox,0,360,50,66 -ifthen $path2/tos_mask/ocean_mask_agrid_"$mod".nc $path2/$mod/hist_rcp2/tas_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/tas_lowlat_ocean_"$mod"_ensmean_186101-209912.nc
done

########################################


############# COMPUTE SURFACE HEAT FLUXES OVER SUBPOLAR OCEAN (50 TO 77 DEGREES NORTH)

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,NorESM1-M,NorESM1-ME,MRI-CGCM3,MPI-ESM-LR,MPI-ESM-MR}
do
for j in {rlds,rsds,rlus,rsus,hfss,hfls}
do
echo $mod $j
cdo ifthen $path2/tos_mask/ocean_mask_agrid_"$mod".nc $path2/$mod/hist_rcp2/"$j"_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/"$j"_ls_"$mod"_ensmean_186101-209912.nc
cdo ifthen $path2/$mod/"$mod"_subpolar.nc -remapnn,$path2/$mod/"$mod"_subpolar.nc $path2/$mod/Method170516/"$j"_ls_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/"$j"_ls_subpolar_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/"$j"_ls_subpolar_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/"$j"_fldmean_ls_subpolar_"$mod"_ensmean_186101-209912.nc
done
done

########################################

################### COMPUTE MEAN OCEAN TEMPERATURE (OVER WHOLE DEPTH)

#FOR LOWER LATITUDES (50 TO 66 DEGREES NORTH)

#Step 1 : cut out low latitudes
#Step 2 : Weighted mean over depth
#Step 3 : Field mean

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2M,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MPI-ESM-LR,MPI-ESM-MR,NorESM1-M,NorESM1-ME}
do
cdo ifthen $path2/AOarea2/lowlat_ocean_grid_"$mod".nc $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_lowlat_"$mod"_ensmean_186101-209912.nc
cdo vertmean $path2/$mod/hist_rcp2/thetao_lowlat_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_vertmean_lowlat_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/thetao_vertmean_lowlat_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_fldmean_lowlat_"$mod"_ensmean_186101-209912.nc
done

#FOR ARCTIC (66 TO 90 DEGREES NORTH)

#Step 1 : cut out Arctic
#Step 2 : Weighted mean over depth
#Step 3 : Field mean

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2M,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MPI-ESM-LR,MPI-ESM-MR,NorESM1-M,NorESM1-ME}
do
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc
cdo vertmean $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_vertmean_Arctic_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/thetao_vertmean_Arctic_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_fldmean_Arctic_"$mod"_ensmean_186101-209912.nc
done

#####
#Special issues
for mod in {MIROC5,MRI-CGCM3}
do
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc -selvar,thetao $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc
cdo vertmean -setctomiss,0 $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_vertmean_Arctic_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/thetao_vertmean_Arctic_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_fldmean_Arctic_"$mod"_ensmean_186101-209912.nc

cdo ifthen $path2/AOarea2/lowlat_ocean_grid_"$mod".nc -selvar,thetao $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_lowlat_"$mod"_ensmean_186101-209912.nc
cdo vertmean -setctomiss,0 $path2/$mod/hist_rcp2/thetao_lowlat_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_vertmean_lowlat_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/thetao_vertmean_lowlat_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_fldmean_lowlat_"$mod"_ensmean_186101-209912.nc
done
###
mod=GISS-E2-R
cdo sellonlatbox,0,360,66,90 -selvar,thetao $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc
cdo vertmean $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_vertmean_Arctic_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/thetao_vertmean_Arctic_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_fldmean_Arctic_"$mod"_ensmean_186101-209912.nc

cdo sellonlatbox,0,360,50,66 -selvar,thetao $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_lowlat_"$mod"_ensmean_186101-209912.nc
cdo vertmean $path2/$mod/hist_rcp2/thetao_lowlat_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_vertmean_lowlat_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/thetao_vertmean_lowlat_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_fldmean_lowlat_"$mod"_ensmean_186101-209912.nc
####
mod=GFDL-ESM2G
cdo ifthen $path2/AOarea2/lowlat_ocean_grid_"$mod".nc -selgrid,2 $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_lowlat_"$mod"_ensmean_186101-209912.nc
cdo vertmean $path2/$mod/hist_rcp2/thetao_lowlat_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_vertmean_lowlat_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/thetao_vertmean_lowlat_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_fldmean_lowlat_"$mod"_ensmean_186101-209912.nc

cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc -selgrid,2 $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc
cdo vertmean $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_vertmean_Arctic_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/thetao_vertmean_Arctic_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_fldmean_Arctic_"$mod"_ensmean_186101-209912.nc

####################################################################################

############################ COMPUTE SURFACE PRESSURE OVER BARENTS SEA

#Step 1 : cut out ocean
#Step 2 : cut out Barents Sea
#Step 3 : Field mean

j=ps
for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
echo $mod
cdo ifthen $path2/AOarea2/Arctic_lsmask10_ocean_"$mod".nc -remapnn,$path2/AOarea2/Arctic_lsmask10_ocean_"$mod".nc $path2/$mod/hist_rcp2/"$j"_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/"$j"_Arctic_remap_"$mod"_ensmean_186101-209912.nc
cdo ifthen $path2/AOarea2/"$mod"_BarentsSeaDomain_ocean.nc $path2/$mod/Method170516/"$j"_Arctic_remap_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/"$j"_BarentsSea_remap_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/"$j"_BarentsSea_remap_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/"$j"_BarentsSea_remap_fldmean_"$mod"_ensmean_186101-209912.nc
done

##############################################################################

################################ MEAN OCEAN TEMPERATURE AT THE DIFFERENT GATEWAYS

#masks for the different regions created in different_regions.sh

#Step 1 : adapt the mask to the oceanic grid 

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
echo $mod
for strait in  {FS,CA,BO,BS,DS,IF,FaSc}
do
cdo remapcon,$path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/"$mod"_"$strait".nc $path2/AOarea2/"$strait"_ocean_grid_"$mod".nc
done
done

#####
#Special issue
mod=GFDL-ESM2G
cdo selgrid,2 $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912_grid2.nc
for line in  {FS,CA,BO,BS,DS,IF,FaSc}
do
cdo remapcon,$path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912_grid2.nc $path2/$mod/"$mod"_"$strait".nc $path2/AOarea2/"$strait"_ocean_grid_"$mod".nc
done
#####


#Step 2 : cut out the strait
#Step 3 : field mean of vertical mean of ocean temperature

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do

for strait in {FS,CA,BO,BS,DS,IF,FaSc}
do

echo $mod $strait

cdo ifthen $path2/AOarea2/"$strait"_ocean_grid_"$mod".nc -selvar,thetao $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_"$strait"_"$mod"_ensmean_186101-209912.nc
cdo vertmean -setctomiss,0 $path2/$mod/hist_rcp2/thetao_"$strait"_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_vertmean_"$strait"_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/thetao_vertmean_"$strait"_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_fldmean_"$strait"_"$mod"_ensmean_186101-209912.nc

done
done

#####
#Special issue
mod=GFDL-ESM2G
for strait in {FS,CA,BO,BS,DS,IF,FaSc}
do

echo $mod $strait

cdo ifthen $path2/AOarea2/"$strait"_ocean_grid_"$mod".nc -selvar,thetao $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912_grid2.nc $path2/$mod/hist_rcp2/thetao_"$strait"_"$mod"_ensmean_186101-209912.nc
cdo vertmean -setctomiss,0 $path2/$mod/hist_rcp2/thetao_"$strait"_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_vertmean_"$strait"_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/thetao_vertmean_"$strait"_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_fldmean_"$strait"_"$mod"_ensmean_186101-209912.nc

done
done

#### REVIEW
## dividing the Fram Strait
for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
echo $mod
for strait in  {FS_r,FS_l}
do
cdo remapcon,$path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/"$mod"_"$strait".nc $path2/AOarea2/"$strait"_ocean_grid_"$mod".nc
done
done

for mod in {GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,MIROC5}
do
for strait in  {FS_r,FS_l}
do
echo $mod
cdo selgrid,2 $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_grid2.nc 
cdo remapcon,$path2/$mod/hist_rcp2/thetao_grid2.nc $path2/$mod/"$mod"_"$strait".nc $path2/AOarea2/"$strait"_ocean_grid_"$mod".nc
done
done

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do

for strait in {FS_r,FS_l}
do

echo $mod $strait

cdo ifthen $path2/AOarea2/"$strait"_ocean_grid_"$mod".nc -selvar,thetao $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/thetao_"$strait"_"$mod"_ensmean_186101-209912.nc
cdo vertmean -sellevel,0/1000 -setctomiss,0 $path2/$mod/hist_rcp2/thetao_"$strait"_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_vertmean_"$strait"_"$mod"_ensmean_186101-209912.nc
cdo fldmean $path2/$mod/Method170516/thetao_vertmean_"$strait"_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/thetao_fldmean_"$strait"_"$mod"_ensmean_186101-209912.nc

done
done


