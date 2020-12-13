#!/bin/bash
#MPI rcp85

################# DECLARE THE PATHS #############################3
path0=/work/mh0033/m300411/DataEB/RAW_DATA/MPI-ESM-LR_rcp85
path1=/work/mh0033/m300411/DataEB/RAW_DATA
###################################################################

mod=MPI-ESM-LR
for i in {1..10}
do
for j in {hfls,hfss,rsds,rsus,rlds,rlus}
do
echo $mod $j "Ensemble $i"
cdo selyear,1861/2005 $path1/$mod/"$j"_Amon_"$mod"_r"$i"i1p1_merged.nc $path0/"$j"_Amon_"$mod"_historical_r"$i"i1.nc
cdo mergetime $path0/"$j"_Amon_"$mod"_*_r"$i"i*.nc $path1/"$mod"_rcp85/"$j"_Amon_"$mod"_r"$i"i1p1_merged.nc
done
done

mod=MPI-ESM-LR
for i in {1..10}
do
for j in {thetao,tos}
do
echo $mod $j "Ensemble $i"
cdo selyear,1861/2005 $path1/$mod/"$j"_Omon_"$mod"_r"$i"i1p1_merged.nc $path0/"$j"_Omon_"$mod"_historical_r"$i"i1.nc
cdo mergetime $path0/"$j"_Omon_"$mod"_*_r"$i"i*.nc $path1/"$mod"_rcp85/"$j"_Omon_"$mod"_r"$i"i1p1_merged.nc
done
done


mod=MPI-ESM-LR
for i in {1..10}
do
for j in {sic,sit}
do
echo $mod $j "Ensemble $i"
cdo selyear,1861/2005 $path1/$mod/"$j"_OImon_"$mod"_r"$i"i1p1_merged.nc $path0/"$j"_OImon_"$mod"_historical_r"$i"i1.nc
cdo mergetime $path0/"$j"_OImon_"$mod"_*_r"$i"i*.nc $path1/"$mod"_rcp85/"$j"_OImon_"$mod"_r"$i"i1p1_merged.nc
done



#only for MPI-ESM-LR - meridional oceanic heat flux and snow depth on sea ice, model output
mod=MPI-ESM-LR
for i in {1..3}
do
for j in {hfx,hfy}
do
echo $mod $j "Ensemble $i"
cdo selyear,1861/2005 $path1/$mod/"$j"_Omon_"$mod"_r"$i"i1p1_merged.nc $path0/"$j"_Omon_"$mod"_historical_r"$i"i1.nc
cdo mergetime $path0/"$j"_Omon_"$mod"_*_r"$i"i*.nc $path1/"$mod"_rcp85/"$j"_Omon_"$mod"_r"$i"i1p1_merged.nc
done
j=snd
echo $mod $j "Ensemble $i"
cdo selyear,1861/2005 $path1/$mod/"$j"_OImon_"$mod"_r"$i"i1p1_merged.nc $path0/"$j"_OImon_"$mod"_historical_r"$i"i1.nc
cdo mergetime $path0/"$j"_OImon_"$mod"_*_r"$i"i*.nc $path1/"$mod"_rcp85/"$j"_OImon_"$mod"_r"$i"i1p1_merged.nc
done


################# DECLARE THE PATHS #############################3
path1=/work/mh0033/m300411/DataEB/RAW_DATA
path2=/work/mh0033/m300411/DataEB/WORK_DATA
###################################################################

############### CUT OUT PERIOD 1861 TO 2099 ##############################
mod=MPI-ESM-LR
for j in {hfls,hfss,rsds,rsus,rlds,rlus,thetao,sit,hfx,hfy,snd}
do
for i in {1..3}
do
echo $mod $j "Ensemble $i"
cdo selyear,1861/2099 $path1/"$mod"_rcp85/"$j"_*_"$mod"_r"$i"i1p1_merged.nc $path2/"$mod"_rcp85/Ensembles/"$j"_"$mod"_r"$i"_186101-209912.nc
done
done

############################################################################

############ BUILD ENSEMBLE MEAN SO THAT WE HAVE ONLY ONE FILE PER VARIABLE PER MODEL ##############

nod=MPI-ESM-LR
for j in {hfls,hfss,rsds,rsus,rlds,rlus,thetao,sit,hfx,hfy,snd}
do
echo $mod $j "Ensemble mean"
cdo ensmean $path2/"$mod"_rcp85/Ensembles/"$j"_"$mod"_r*_186101-209912.nc $path2/"$mod"_rcp85/hist_rcp2/"$j"_"$mod"_ensmean_186101-209912.nc
done


###################################################################################

#Arctic ocean grids, grid areas and land sea masks are computed from AOmask.sh

########################################################################
######### ATMOSPHERE ##########
##########################################################################

#Step 1 : remap atmospheric variables to oceanic grid for comparison
#Step 2 : mean over the whole area

mod=MPI-ESM-LR
for j in {hfls,hfss,rsds,rsus,rlds,rlus}
do
echo $mod $j
echo "Remap + Arctic $mod $j"
cdo ifthen $path2/AOarea2/Arctic_lsmask10_ocean_"$mod".nc -remapnn,$path2/AOarea2/Arctic_lsmask10_ocean_"$mod".nc $path2/"$mod"_rcp85/hist_rcp2/"$j"_"$mod"_ensmean_186101-209912.nc $path2/"$mod"_rcp85/Method170516/"$j"_Arctic_remap_"$mod"_ensmean_186101-209912.nc
echo "fldmean $mod $j"
cdo fldmean $path2/"$mod"_rcp85/Method170516/"$j"_Arctic_remap_"$mod"_ensmean_186101-209912.nc $path2/"$mod"_rcp85/Method170516/"$j"_Arctic_monfldmean_remap_"$mod"_ensmean_186101-209912.nc
done


#############################################################
####### SEA ICE #########
##############################################################

#Compute total sea-ice area and total sea-ice volume by selecting the domain, multiplying by gridarea and making the sum over the domain

mod=MPI-ESM-LR
echo $mod
cdo fldsum -mul -ifthen $path2/AOarea2/Arctic_lsmask10_seaice_"$mod".nc $path2/"$mod"_rcp85/hist_rcp2/sit_"$mod"_ensmean_186101-209912.nc $path2/AOarea2/Arctic_seaice_gridarea_"$mod".nc $path2/"$mod"_rcp85/Method170516/totsiv_Arctic_fldmean_"$mod"_ensmean_186101-209912.nc

#same method for total snow volume for MPI-ESM-LR
mod=MPI-ESM-LR
echo $mod
cdo fldsum -mul -ifthen $path2/AOarea2/Arctic_lsmask10_seaice_"$mod".nc $path2/"$mod"_rcp85/hist_rcp2/snd_"$mod"_ensmean_186101-209912.nc $path2/AOarea2/Arctic_ocean_gridarea_"$mod".nc $path2/"$mod"_rcp85/Method170516/totsnv_Arctic_fldmean_"$mod"_ensmean_186101-209912.nc


##################################################################
####### OCEAN #############
####################################################################

#for the heat content
#Step 1 : make a weighted sum over product of temperature and volume of the cells (volcello can be found in CMIP5 archive)
#Step 2 : select Arctic
#Step 3 : sum over the domain of the product temperature x volume (the heat content is then calculated in the python file later)

mod=MPI-ESM-LR
echo $mod
cdo vertsum -mul $path1/volcello_files/volcello_fx_"$mod"_historical_r0i0p0.nc $path2/"$mod"_rcp85/hist_rcp2/thetao_"$mod"_ensmean_186101-209912.nc $path2/"$mod"_rcp85/hist_rcp2/heatcolumn_"$mod"_ensmean_186101-209912.nc
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/"$mod"_rcp85/hist_rcp2/heatcolumn_"$mod"_ensmean_186101-209912.nc $path2/"$mod"_rcp85/hist_rcp2/heatcolumn_Arctic_"$mod"_ensmean_186101-209912.nc
cdo fldsum $path2/"$mod"_rcp85/hist_rcp2/heatcolumn_Arctic_"$mod"_ensmean_186101-209912.nc $path2/"$mod"_rcp85/Method170516/heatcolumn_Arctic_monfldsum_"$mod"_ensmean_186101-209912.nc
done



