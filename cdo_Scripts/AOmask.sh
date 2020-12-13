#!/bin/bash

#DEFINE THE ARCTIC OCEAN MASK
#compute the grid areas, land sea masks...


################# DECLARE THE PATHS #############################3
path1=/work/mh0033/m300411/DataEB/RAW_DATA
path2=/work/mh0033/m300411/DataEB/WORK_DATA
###################################################################

####### CUT OUT OCEAN AREA WITH THE HELP OF SEA SURFACE TEMPERATURE #################

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,GISS-E2-R-CC,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M}
do
cdo seltimestep,1 -selvar,tos $path2/tos_mask/tos_Omon_"$mod"_historical_r1i1p1_*.nc $path2/tos_mask/tos_mask_"$mod".nc
done

mod=CSIRO-Mk3-6-0
cdo seltimestep,1 -sellevidx,1 -chname,thetao,tos $path2/$mod/hist_rcp/thetao_"$mod"_ensmean_186101-209912.nc $path2/tos_mask/tos_mask_"$mod".nc

#landsea mask for atmosphere
cdo remapcon,$path2/$mod/hist_rcp2/tas_"$mod"_ensmean_186101-209912.nc $path2/tos_mask/tos_mask_"$mod".nc $path2/tos_mask/ocean_mask_agrid_"$mod".nc

#####################################################################################


##################################################################

### ARCTIC MASK ON OCEAN GRID #######################################################
for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,GISS-E2-R-CC,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
#setup the sellonlatbox for ocean
cdo remapcon,$path2/tos_mask/tos_mask_"$mod".nc $path2/$mod/"$mod"_grid.nc $path2/AOarea2/Arctic_ocean_grid_"$mod".nc
#get gridarea for gridboxes in Arctic
cdo gridarea $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/AOarea2/Arctic_ocean_gridarea_"$mod".nc

#lanseamask with 1 over ocean and missing value otherwise
cdo ifthenc,1 -ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/tos_mask/tos_mask_"$mod".nc $path2/AOarea2/Arctic_lsmask10_ocean_"$mod".nc
#multiply gridarea by 1/0 mask and make a fldsum => Total area of Arctic Ocean
cdo fldsum -mul $path2/AOarea2/Arctic_lsmask10_ocean_"$mod".nc $path2/AOarea2/Arctic_ocean_gridarea_"$mod".nc $path2/AOarea2/Arctic_ocean_totalarea_"$mod".nc
done

#####
#Special issues
mod=GISS-E2-R
cdo fldsum -mul -ifthenc,1 -sellonlatbox,0,360,66,90 $path2/$mod/hist_rcp/heatcolumn_"$mod"_ensmean_186101-209912.nc -gridarea -sellonlatbox,0,360,66,90 $path2/$mod/hist_rcp/heatcolumn_"$mod"_ensmean_186101-209912.nc $path2/AOarea2/Arctic_ocean_totalarea_"$mod".nc
#####
########################################################################################################

######## ARCTIC MASK ON SEA ICE GRID ###################################################################

#The models where ocean and sea ice grid are the same
for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
cp $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/AOarea2/Arctic_seaice_grid_"$mod".nc
cp $path2/AOarea2/Arctic_ocean_gridarea_"$mod".nc $path2/AOarea2/Arctic_seaice_gridarea_"$mod".nc
cp $path2/AOarea2/Arctic_lsmask10_ocean_"$mod".nc $path2/AOarea2/Arctic_lsmask10_seaice_"$mod".nc
cp $path2/AOarea2/Arctic_ocean_totalarea_"$mod".nc $path2/AOarea2/Arctic_seaice_totalarea_"$mod".nc
done

#Other models
for mod in {CSIRO-Mk3-6-0,GISS-E2-R}
do
#get gridarea for gridboxes in Arctic
cdo gridarea -sellonlatbox,0,360,66,90 $path2/$mod/hist_rcp/sic_"$mod"_ensmean_186101-209912.nc $path2/AOarea2/Arctic_seaice_gridarea_"$mod".nc
#setup a landsea mask
cdo remapcon,$path2/AOarea2/Arctic_seaice_gridarea_"$mod".nc $path2/tos_mask/tos_mask_"$mod".nc  $path2/AOarea2/Arctic_lsmask_seaice_"$mod".nc

#lanseamask with 1 over ocean and missing value otherwise
cdo ifthenc,1 $path2/AOarea2/Arctic_lsmask_seaice_"$mod".nc $path2/AOarea2/Arctic_lsmask10_seaice_"$mod".nc
#multiply gridarea by 1/0 mask and make a fldsum => Total area of Arctic Ocean
cdo fldsum -mul $path2/AOarea2/Arctic_lsmask10_seaice_"$mod".nc $path2/AOarea2/Arctic_seaice_gridarea_"$mod".nc $path2/AOarea2/Arctic_seaice_totalarea_"$mod".nc
done

######################################################################################################


