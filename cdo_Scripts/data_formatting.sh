#!/bin/bash

#RAW DATA FORMATTING
#merge all single files into one long file : historical and RCP4.5 files

################# DECLARE THE PATHS #############################3
path0=path_with_downloaded_data
path1=/work/mh0033/m300411/DataEB/RAW_DATA
###################################################################

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
for i in {1..10}
do
#for j in {hfls,hfss,rsds,rsus,rlds,rlus,tas,ps} #rlut,rsut,rsdt
for j in {rlut,rsut,rsdt} 
do
echo $mod $j "Ensemble $i"
SKIP_SAME_TIME=1 ; export SKIP_SAME_TIME
cdo mergetime $path0/"$j"_Amon_"$mod"_*_r"$i"i*.nc $path1/$mod/"$j"_Amon_"$mod"_r"$i"i1p1_merged.nc
done
done
done

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
for i in {1..10}
do
for j in {mfo,tos,thetao}
do
echo $mod $j "Ensemble $i"
cdo mergetime $path0/"$j"_Omon_"$mod"_*_r"$i"i*.nc $path1/$mod/"$j"_Omon_"$mod"_r"$i"i1p1_merged.nc
done
done
done

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
for i in {1..10}
do
for j in {sic,sit,transifs}
do
echo $mod $j "Ensemble $i"
cdo mergetime $path0/"$j"_OImon_"$mod"_*_r"$i"i*.nc $path1/$mod/"$j"_OImon_"$mod"_r"$i"i1p1_merged.nc
done
done
done

#only for MPI-ESM-LR - meridional oceanic heat flux and snow depth on sea ice, model output
mod=MPI-ESM-LR
for i in {1..3}
do
for j in {hfx,hfy}
echo $mod $j "Ensemble $i"
cdo mergetime $path0/"$j"_Omon_"$mod"_*_r"$i"i*.nc $path1/$mod/"$j"_Omon_"$mod"_r"$i"i1p1_merged.nc
done
j=snd
echo $mod $j "Ensemble $i"
cdo mergetime $path0/"$j"_OImon_"$mod"_*_r"$i"i*.nc $path1/$mod/"$j"_OImon_"$mod"_r"$i"i1p1_merged.nc
done


#for water vapor pressure

for mod in {ACCESS1-0,ACCESS1-3,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,IPSL-CM5A-MR,IPSL-CM5B-LR,MIROC5,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M,NorESM1-ME}
do
for i in {1..10}
do
echo $mod $j "Ensemble $i"
cdo mul $path1/$mod/ps_Amon_"$mod"_r"$i"i1p1_merged.nc $path1/$mod/prw_Amon_"$mod"_r"$i"i1p1_merged.nc $path1/$mod/psxprw_Amon_"$mod"_r"$i"i1p1_merged.nc
done
done

for mod in {GFDL-ESM2G,GFDL-ESM2M}
do
i=1
cdo mul -selvar,ps $path1/$mod/ps_Amon_"$mod"_r"$i"i1p1_merged.nc -selvar,prw $path1/$mod/prw_Amon_"$mod"_r"$i"i1p1_merged.nc $path1/$mod/psxprw_Amon_"$mod"_r"$i"i1p1_merged.nc
done
