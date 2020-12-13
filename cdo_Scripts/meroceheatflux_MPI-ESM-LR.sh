#!/bin/bash

#MERIDIONAL OCEANIC HEAT FLUX FROM MPI-ESM-LR, LOOKING ONLY AT 66N
#use uline_66.nc and vline_66.nc, and uline_b.nc and vline_b.nc

################# DECLARE THE PATHS #############################3
path1=/work/mh0033/m300411/DataEB/RAW_DATA
path2=/work/mh0033/m300411/DataEB/WORK_DATA
###################################################################

mod=MPI-ESM-LR

#masks for the two regions of the 66N line

cdo selvar,secmapv_atlantic_60n $path2/hf_ocean_test/HelmuthsDing/2nd/hel16130-PGG_mpiom_map_0000-01-01_0000-01-01.nc $path2/hf_ocean/vline_66.nc
cdo selvar,secmapu_atlantic_60n $path2/hf_ocean_test/HelmuthsDing/2nd/hel16130-PGG_mpiom_map_0000-01-01_0000-01-01.nc $path2/hf_ocean/uline_66.nc
cdo selvar,secmapv_bering_strait $path2/hf_ocean_test/HelmuthsDing/2nd/hel16130-PGG_mpiom_map_0000-01-01_0000-01-01.nc $path2/hf_ocean/vline_b.nc
cdo selvar,secmapu_bering_strait $path2/hf_ocean_test/HelmuthsDing/2nd/hel16130-PGG_mpiom_map_0000-01-01_0000-01-01.nc $path2/hf_ocean/uline_b.nc

#Step 1 : multiply ocean heat transport in x-direction and y-direction with the corresponding factor (-1 if southward, +1 if northward) using a mask for the 66N line (for more information, contact the author)
#Step 2 : add the x- and y- components to get the total over each strait

#Fram Strait region 66N
cdo mul $path2/hf_ocean/uline_66.nc $path2/$mod/hist_rcp2/hfx_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hfx_mulu_"$mod"_ensmean_186101-209912.nc
cdo mul $path2/hf_ocean/vline_66.nc $path2/$mod/hist_rcp2/hfy_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hfy_mulv_"$mod"_ensmean_186101-209912.nc
cdo add $path2/hf_ocean/hfx_mulu_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hfy_mulv_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hfxhfy_"$mod"_ensmean_186101-209912.nc 
cdo fldsum $path2/hf_ocean/hfxhfy_"$mod"_ensmean_186101-209912.nc  $path2/hf_ocean/hf_f_fldsum_"$mod"_ensmean_186101-209912.nc

#Bering Strait 66N 
cdo mul $path2/hf_ocean/uline_b.nc $path2/$mod/hist_rcp2/hfx_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hfx_mulu_b_"$mod"_ensmean_186101-209912.nc
cdo mul $path2/hf_ocean/vline_b.nc $path2/$mod/hist_rcp2/hfy_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hfy_mulv_b_"$mod"_ensmean_186101-209912.nc
cdo fldsum -chname,hfy,hfx $path2/hf_ocean/hfy_mulv_b_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hf_b_fldsum_"$mod"_ensmean_186101-209912.nc

#Step 3 : add both straits
cdo add $path2/hf_ocean/hf_b_fldsum_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hf_f_fldsum_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hf_tot_fldsum_mon_"$mod"_ensmean_186101-209912.nc

###############################for rcp85
#Fram Strait region 66N
cdo mul $path2/hf_ocean/uline_66.nc $path2/"$mod"_rcp85/hist_rcp2/hfx_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hfx_mulu_"$mod"_ensmean_186101-209912_rcp85.nc
cdo mul $path2/hf_ocean/vline_66.nc $path2/"$mod"_rcp85/hist_rcp2/hfy_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hfy_mulv_"$mod"_ensmean_186101-209912_rcp85.nc
cdo add $path2/hf_ocean/hfx_mulu_"$mod"_ensmean_186101-209912_rcp85.nc $path2/hf_ocean/hfy_mulv_"$mod"_ensmean_186101-209912_rcp85.nc $path2/hf_ocean/hfxhfy_"$mod"_ensmean_186101-209912_rcp85.nc 
cdo fldsum $path2/hf_ocean/hfxhfy_"$mod"_ensmean_186101-209912_rcp85.nc  $path2/hf_ocean/hf_f_fldsum_"$mod"_ensmean_186101-209912_rcp85.nc

#Bering Strait 66N 
cdo mul $path2/hf_ocean/uline_b.nc $path2/"$mod"_rcp85/hist_rcp2/hfx_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hfx_mulu_b_"$mod"_ensmean_186101-209912_rcp85.nc
cdo mul $path2/hf_ocean/vline_b.nc $path2/"$mod"_rcp85/hist_rcp2/hfy_"$mod"_ensmean_186101-209912.nc $path2/hf_ocean/hfy_mulv_b_"$mod"_ensmean_186101-209912_rcp85.nc
cdo fldsum -chname,hfy,hfx $path2/hf_ocean/hfy_mulv_b_"$mod"_ensmean_186101-209912_rcp85.nc $path2/hf_ocean/hf_b_fldsum_"$mod"_ensmean_186101-209912_rcp85.nc

#Step 3 : add both straits
cdo add $path2/hf_ocean/hf_b_fldsum_"$mod"_ensmean_186101-209912_rcp85.nc $path2/hf_ocean/hf_f_fldsum_"$mod"_ensmean_186101-209912_rcp85.nc $path2/hf_ocean/hf_tot_fldsum_mon_"$mod"_ensmean_186101-209912_rcp85.nc

