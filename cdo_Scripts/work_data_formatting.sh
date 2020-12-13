#!/bin/bash

#WORK DATA FORMATTING
#format the data to work with it

################# DECLARE THE PATHS #############################3
path1=/work/mh0033/m300411/DataEB/RAW_DATA
path2=/work/mh0033/m300411/DataEB/WORK_DATA
###################################################################

############### CUT OUT PERIOD 1861 TO 2099 ##############################
for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
#for j in {hfls,hfss,rsds,rsus,rlds,rlus,tas,ps,mfo,tos,thetao,sic,sit,transifs}
for j in {rlut,rsut,rsdt} 
do
for i in {1..10}
do
echo $mod $j "Ensemble $i"
cdo selyear,1861/2099 $path1/$mod/"$j"_*_"$mod"_r"$i"i1p1_merged.nc $path2/$mod/Ensembles/"$j"_"$mod"_r"$i"_186101-209912.nc
done
done
done

#only for MPI-ESM-LR - meridional oceanic heat flux and snow depth on sea ice, model output
mod=MPI-ESM-LR
for i in {1..3}
do
for j in {hfx,hfy,snd}
echo $mod $j "Ensemble $i"
cdo selyear,1861/2099 $path1/$mod/"$j"_*_"$mod"_r"$i"i1p1_merged.nc $path2/$mod/Ensembles/"$j"_"$mod"_r"$i"_186101-209912.nc
done
done
############################################################################

############ BUILD ENSEMBLE MEAN SO THAT WE HAVE ONLY ONE FILE PER VARIABLE PER MODEL ##############

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
#for j in {hfls,hfss,rsds,rsus,rlds,rlus,tas,ps,mfo,tos,thetao,sic,sit,transifs}
for j in {rlut,rsut,rsdt} 
do
echo $mod $j "Ensemble mean"
cdo ensmean $path2/$mod/Ensembles/"$j"_"$mod"_r*_186101-209912.nc $path2/$mod/hist_rcp2/"$j"_"$mod"_ensmean_186101-209912.nc
done
done

#only for MPI-ESM-LR - meridional oceanic heat flux and snow depth on sea ice, model output
mod=MPI-ESM-LR
for j in {hfx,hfy,snd}
echo $mod $j "Ensemble mean"
cdo ensmean $path2/$mod/Ensembles/"$j"_"$mod"_r*_186101-209912.nc $path2/$mod/hist_rcp2/"$j"_"$mod"_ensmean_186101-209912.nc
done
###################################################################################

#Arctic ocean grids, grid areas and land sea masks are computed from AOmask.sh

########################################################################
######### ATMOSPHERE ##########
##########################################################################

#Step 1 : remap atmospheric variables to oceanic grid for comparison
#Step 2 : mean over the whole area

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
#for j in {hfls,hfss,rsds,rsus,rlds,rlus,tas}
for j in {rlut,rsut,rsdt} 
do
echo $mod $j
echo "Remap + Arctic $mod $j"
cdo ifthen $path2/AOarea2/Arctic_lsmask10_ocean_"$mod".nc -remapnn,$path2/AOarea2/Arctic_lsmask10_ocean_"$mod".nc $path2/$mod/hist_rcp2/"$j"_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/"$j"_Arctic_remap_"$mod"_ensmean_186101-209912.nc
echo "fldmean $mod $j"
cdo fldmean $path2/$mod/Method170516/"$j"_Arctic_remap_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/"$j"_Arctic_monfldmean_remap_"$mod"_ensmean_186101-209912.nc
done
done

#############################################################
####### SEA ICE #########
##############################################################

#Compute total sea-ice area and total sea-ice volume by selecting the domain, multiplying by gridarea and making the sum over the domain

for mod in {ACCESS1-0,ACCESS1-3,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,FGOALS-g2,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
echo $mod
cdo fldsum -mul -divc,100 -ifthen $path2/AOarea2/Arctic_lsmask10_seaice_"$mod".nc $path2/$mod/hist_rcp2/sic_"$mod"_ensmean_186101-209912.nc $path2/AOarea2/Arctic_seaice_gridarea_"$mod".nc $path2/$mod/Method170516/totsia_Arctic_fldmean_"$mod"_ensmean_186101-209912.nc
cdo fldsum -mul -ifthen $path2/AOarea2/Arctic_lsmask10_seaice_"$mod".nc $path2/$mod/hist_rcp2/sit_"$mod"_ensmean_186101-209912.nc $path2/AOarea2/Arctic_seaice_gridarea_"$mod".nc $path2/$mod/Method170516/totsiv_Arctic_fldmean_"$mod"_ensmean_186101-209912.nc
done

#same method for total snow volume for MPI-ESM-LR
mod=MPI-ESM-LR
echo $mod
cdo fldsum -mul -ifthen $path2/AOarea2/Arctic_lsmask10_seaice_"$mod".nc $path2/$mod/hist_rcp2/snd_"$mod"_ensmean_186101-209912.nc $path2/AOarea2/Arctic_ocean_gridarea_"$mod".nc $path2/$mod/Method170516/totsnv_Arctic_fldmean_"$mod"_ensmean_186101-209912.nc

#####
#Some models have some special issues :
for mod in {CanESM2,CSIRO-Mk3-6-0,GISS-E2-R}
do
echo $mod
cdo fldsum -mul -sellonlatbox,0,360,66,90 $path2/$mod/hist_rcp2/sit_"$mod"_ensmean_186101-209912.nc $path2/AOarea2/Arctic_seaice_gridarea_"$mod".nc $path2/$mod/Method170516/totsiv_Arctic_fldmean_"$mod"_ensmean_186101-209912.nc
cdo fldsum -mul -divc,100 -sellonlatbox,0,360,66,90 $path2/$mod/hist_rcp2/sic_"$mod"_ensmean_186101-209912.nc $path2/AOarea2/Arctic_seaice_gridarea_"$mod".nc $path2/$mod/Method170516/totsia_Arctic_fldmean_"$mod"_ensmean_186101-209912.nc
done

for mod in {GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M}
do
cdo fldsum -mul -ifthen $path2/AOarea2/Arctic_lsmask10_seaice_"$mod".nc -selgrid,2 $path2/$mod/hist_rcp2/sit_"$mod"_ensmean_186101-209912.nc $path2/AOarea2/Arctic_seaice_gridarea_"$mod".nc $path2/$mod/Method170516/totsiv_Arctic_fldmean_"$mod"_ensmean_186101-209912.nc
cdo fldsum -mul -divc,100 -ifthen $path2/AOarea2/Arctic_lsmask10_seaice_"$mod".nc -selgrid,2 $path2/$mod/hist_rcp2/sic_"$mod"_ensmean_186101-209912.nc $path2/AOarea2/Arctic_seaice_gridarea_"$mod".nc $path2/$mod/Method170516/totsia_Arctic_fldmean_"$mod"_ensmean_186101-209912.nc
done 
#####

##################################################################
####### OCEAN #############
####################################################################

#for the heat content
#Step 1 : make a weighted sum over product of temperature and volume of the cells (volcello can be found in CMIP5 archive)
#Step 2 : select Arctic
#Step 3 : sum over the domain of the product temperature x volume (the heat content is then calculated in the python file later)

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-ESM2G,HadGEM2-CC,HadGEM2-ES,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3}
do
echo $mod
cdo vertsum -mul $path1/volcello_files/volcello_fx_"$mod"_historical_r0i0p0.nc $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/heatcolumn_"$mod"_ensmean_186101-209912.nc
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/$mod/hist_rcp2/heatcolumn_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/heatcolumn_Arctic_"$mod"_ensmean_186101-209912.nc
cdo fldsum $path2/$mod/hist_rcp2/heatcolumn_Arctic_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/heatcolumn_Arctic_monfldsum_"$mod"_ensmean_186101-209912.nc
done

#####
#Special issues
mod=GFDL-CM3
echo $mod
cdo vertsum -mul $path1/volcello_files/volcello_fx_"$mod"_historical_r0i0p0.nc -selgrid,2 $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/heatcolumn_"$mod"_ensmean_186101-209912.nc
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/$mod/hist_rcp2/heatcolumn_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/heatcolumn_Arctic_"$mod"_ensmean_186101-209912.nc

mod=MIROC5
echo $mod
cdo vertsum -mul $path1/volcello_files/volcello_fx_"$mod"_historical_r0i0p0.nc -selvar,thetao $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/heatcolumn_"$mod"_ensmean_186101-209912.nc
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/$mod/hist_rcp2/heatcolumn_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/heatcolumn_Arctic_"$mod"_ensmean_186101-209912.nc

mod=GFDL-ESM2M
echo $mod
cdo vertsum -mul $path1/volcello_files/volcello_fx_"$mod"_piControl_r0i0p0.nc -selgrid,2 $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/heatcolumn_"$mod"_ensmean_186101-209912.nc
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/$mod/hist_rcp2/heatcolumn_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/heatcolumn_Arctic_"$mod"_ensmean_186101-209912.nc


#Models without or with erroneous volcello files, use depth and weighted mean over column
for mod in {CESM1-CAM5,GISS-E2-R,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,NorESM1-M,NorESM1-ME}
do
cdo mul $path1/deptho/deptho_fx_"$mod"_historical_r0i0p0.nc -vertmean -selvar,thetao -setctomiss,0  $path2/$mod/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/heatcolumn_"$mod"_ensmean_186101-209912.nc
cdo mul $path2/AOarea2/Arctic_ocean_gridarea_"$mod".nc -ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/$mod/hist_rcp2/heatcolumn_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/heatcolumn_Arctic_"$mod"_ensmean_186101-209912.nc
done
#####

#for the sea surface temperature
j=tos
for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,NorESM1-M,NorESM1-ME,MPI-ESM-LR,MPI-ESM-MR}
do
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/$mod/hist_rcp2/tos_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/tos_Arctic_"$mod"_ensmean_186101-209912.nc
cdo fldmean -setctomiss,0 $path2/$mod/hist_rcp2/"$j"_Arctic_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/"$j"_Arctic_monfldmean_"$mod"_ensmean_186101-209912.nc
done

#####
#Special issues
j=tos
for mod in {MIROC5,MRI-CGCM3}
do
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/$mod/hist_rcp2/tos_"$mod"_ensmean_186101-209912.nc $path2/$mod/hist_rcp2/tos_Arctic_"$mod"_ensmean_186101-209912.nc
cdo fldmean -setctomiss,0 $path2/$mod/hist_rcp2/"$j"_Arctic_"$mod"_ensmean_186101-209912.nc $path2/$mod/Method170516/"$j"_Arctic_monfldmean_"$mod"_ensmean_186101-209912.nc
done
######



