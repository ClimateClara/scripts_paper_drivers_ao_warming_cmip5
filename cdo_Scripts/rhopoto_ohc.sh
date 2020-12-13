#!/bin/bash       

################# DECLARE THE PATHS #############################3
path0=/work/mh0033/m300411/DataEB/RAW_DATA/density_files
path1=/work/mh0033/m300411/DataEB/RAW_DATA
###################################################################

for mod in {CESM1-CAM5,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M}
do
for i in {1..3}
do
for j in {rhopoto,ccc}
do
echo $mod $j "Ensemble $i"
cdo mergetime $path0/"$j"_Omon_"$mod"_*_r"$i"i*.nc $path1/$mod/"$j"_Omon_"$mod"_r"$i"i1p1_merged.nc
done
done
done

#GFDL-ESM2G, GFDL-ESM2M and CESM1-CAM5 r1 full
#rhopoto*cp*integral(thetao*volcello)

################# DECLARE THE PATHS #############################3
path1=/work/mh0033/m300411/DataEB/RAW_DATA
path2=/work/mh0033/m300411/DataEB/WORK_DATA
###################################################################
for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
cdo mul -gridarea $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc $path1/deptho/deptho_fx_"$mod"_historical_r0i0p0.nc $path1/deptho/depvol/depvol_fx_"$mod"_historical_r0i0p0.nc
done

mod=GISS-E2-R

for mod in {GFDL-CM3,GFDL-ESM2M}
do
cdo mul -gridarea -selgrid,2 $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc $path1/deptho/deptho_fx_"$mod"_historical_r0i0p0.nc $path1/deptho/depvol/depvol_fx_"$mod"_historical_r0i0p0.nc
done


mod=GFDL-ESM2G
echo "1"
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc $path2/density_stuff/$mod/thetao_Arctic_"$mod"_ensmean_186101-209912.nc
echo "2"
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc -selyear,1861/2099 $path1/$mod/rhopoto_Omon_"$mod"_r"$i"i1p1_merged.nc  $path2/density_stuff/$mod/rhopoto_Arctic_"$mod"_ensmean_186101-209912.nc
echo "3"
cdo mul -selgrid,2 $path2/density_stuff/$mod/rhopoto_Arctic_"$mod"_ensmean_186101-209912.nc $path2/density_stuff/$mod/thetao_Arctic_"$mod"_ensmean_186101-209912.nc  $path2/density_stuff/$mod/txrho_"$mod".nc
echo "4"
cdo mul -gridarea $path1/deptho/deptho_fx_GFDL-ESM2G_historical_r0i0p0.nc $path1/deptho/deptho_fx_GFDL-ESM2G_historical_r0i0p0.nc $path1/deptho/volcell_fx_GFDL-ESM2G_historical_r0i0p0.nc 
echo "5"
cdo mul -vertmean $path2/density_stuff/$mod/txrho_"$mod".nc $path2/density_stuff/$mod/volcell_"$mod".nc  $path2/density_stuff/$mod/txrho_allcolumns_pot_Arctic_"$mod"_ensmean_186101-209912.nc
echo "6"
cdo fldmean -mulc,3992.1 $path2/density_stuff/$mod/txrho_allcolumns_pot_Arctic_"$mod"_ensmean_186101-209912.nc $path2/density_stuff/$mod/ohc_tot_pot_Arctic_"$mod"_ensmean_186101-209912.nc

#test for NorESM1-M
mod=NorESM1-M
cdo mul -gridarea $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc  $path1/deptho/deptho_fx_"$mod"_historical_r0i0p0.nc $path2/density_stuff/$mod/volcell_"$mod".nc
cdo vertsum -mul $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc $path2/density_stuff/$mod/volcell_"$mod".nc $path2/density_stuff/$mod/compare_methods.nc

mod=GFDL-ESM2M
echo "1"
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc $path2/density_stuff/$mod/thetao_Arctic_"$mod"_ensmean_186101-209912.nc
echo "2"
cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc -selyear,1861/2099 $path1/$mod/rhopoto_Omon_"$mod"_r"$i"i1p1_merged.nc  $path2/density_stuff/$mod/rhopoto_Arctic_"$mod"_ensmean_186101-209912.nc
echo "3"
cdo vertmean -mul $path2/density_stuff/$mod/thetao_Arctic_"$mod"_ensmean_186101-209912.nc -selgrid,2 $path2/density_stuff/$mod/rhopoto_Arctic_"$mod"_ensmean_186101-209912.nc  $path2/density_stuff/$mod/txrho_column_pot_Arctic_"$mod"_ensmean_186101-209912.nc
echo "4"
cdo mulc,3992.1 -mul $path1/volcello_files/volcello_fx_GFDL-ESM2M_piControl_r0i0p0.nc $path2/density_stuff/$mod/txrho_column_pot_Arctic_"$mod"_ensmean_186101-209912.nc $path2/density_stuff/$mod/ohc_tot_pot_Arctic_"$mod"_ensmean_186101-209912.nc

#mod=GFDL-CM3
#echo "1"
#cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc $path2/$mod/hist_rcp2/thetao_Arctic_"$mod"_ensmean_186101-209912.nc $path2/density_stuff/$mod/thetao_Arctic_"$mod"_ensmean_186101-209912.nc
#echo "2"
#cdo ifthen $path2/AOarea2/Arctic_ocean_grid_"$mod".nc -selyear,1861/2099 $path1/$mod/rhopoto_Omon_"$mod"_r"$i"i1p1_merged.nc  $path2/density_stuff/$mod/rhopoto_Arctic_"$mod"_ensmean_186101-209912.nc
#echo "3"
#cdo vertsum -mul $path1/volcello_files/volcello_fx_"$mod"_historical_r0i0p0.nc -mul $path2/density_stuff/$mod/thetao_Arctic_"$mod"_ensmean_186101-209912.nc -selgrid,2 $path2/density_stuff/$mod/rhopoto_Arctic_"$mod"_ensmean_186101-209912.nc  $path2/density_stuff/$mod/txrho_column_pot_Arctic_"$mod"_ensmean_186101-209912.nc
#echo "4"
#cdo fldmean -mulc,3992.1 $path2/density_stuff/$mod/txrho_column_pot_Arctic_"$mod"_ensmean_186101-209912.nc $path2/density_stuff/$mod/ohc_tot_pot_Arctic_"$mod"_ensmean_186101-209912.nc


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
